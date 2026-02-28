"""
matrix.py — STAR count file parsing and count matrix construction
=================================================================

FILE MATCHING STRATEGY
-----------------------
When GDC returns a bulk ZIP from the /data endpoint, it unpacks to:

    raw_counts/
        <file_id_uuid>/         <- directory named with file_id UUID
            <filename>.tsv      <- the STAR counts file

We match each extracted TSV back to its sample metadata using this
priority order:

    1. UUID directory name  (most reliable — always present and unique)
       The parent directory of the TSV is the file_id UUID from our
       metadata table. We look this up in metadata_df['file_id'].

    2. Filename match       (fallback — works but filenames can collide
       across projects in edge cases)

    3. UNMATCHED prefix     (last resort — file is parsed but gets a
       temporary ID. These will be flagged in the output.)

COLUMN DETECTION
----------------
STAR-Counts files from the GDC have a header row followed by data rows.
The column containing "unstranded" counts is identified dynamically by
reading the header — never assumed to be at a fixed column index.

Current GDC STAR format:
    gene_id  gene_name  gene_type  unstranded  stranded_first  stranded_second
    (indices 0          1          2           3               4               5)

Older GDC STAR format:
    gene_id  unstranded  stranded_first  stranded_second
    (indices 0           1               2               3)

N_ SUMMARY ROWS
---------------
STAR counts files include 4 alignment summary rows that must be excluded:
    N_unmapped, N_multimapping, N_noFeature, N_ambiguous
These start with "N_" and are not gene-level counts.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from tcga_downloader.constants import (
    UNSTRANDED_SYNONYMS, STRANDED1_SYNONYMS, STRANDED2_SYNONYMS,
    is_valid_gdc_uuid,
)
from tcga_downloader.exceptions import ColumnDetectionError, GDCError


# =============================================================================
#   Column detection
# =============================================================================

def detect_column_map(filepath: Path) -> dict[str, int]:
    """
    Read the header of a STAR counts TSV and return column index mapping.

    Returns dict with keys 'unstranded' (required), 'stranded_first',
    'stranded_second' (optional). Values are zero-based column indices.

    Raises ColumnDetectionError if the unstranded column cannot be found.
    """
    try:
        lines = _read_lines(filepath, max_lines=20)
    except Exception as e:
        raise GDCError(
            f"Cannot open file for column detection: {filepath.name}",
            fix=(
                "The file may be corrupted. Try deleting the raw_counts\n"
                "folder and re-running to trigger a fresh download."
            ),
            step="column detection",
        ) from e

    header_line = None
    first_data_line = None
    for line in lines:
        line = line.strip()
        if not line or line.startswith("N_") or line.startswith("#"):
            continue
        if header_line is None:
            header_line = line
        elif first_data_line is None:
            first_data_line = line
            break

    if header_line is None:
        raise ColumnDetectionError(str(filepath), [])

    cols_lower = [c.strip().lower() for c in header_line.split("\t")]
    col_map: dict[str, int] = {}

    # Try named synonyms first
    for name in UNSTRANDED_SYNONYMS:
        if name in cols_lower:
            col_map["unstranded"] = cols_lower.index(name)
            break

    for name in STRANDED1_SYNONYMS:
        if name in cols_lower:
            col_map["stranded_first"] = cols_lower.index(name)
            break

    for name in STRANDED2_SYNONYMS:
        if name in cols_lower:
            col_map["stranded_second"] = cols_lower.index(name)
            break

    # Fallback: scan ALL candidate lines for the first integer-valued column.
    # This handles both unrecognised header names AND headerless files where
    # candidate_lines[0] is already a data row (e.g. "ENSG00000000003\t1500...").
    if "unstranded" not in col_map:
        for candidate in [header_line, first_data_line]:
            if candidate is None:
                continue
            for idx, val in enumerate(candidate.split("\t")[1:], start=1):
                try:
                    int(val.strip())
                    col_map["unstranded"] = idx
                    break
                except ValueError:
                    continue
            if "unstranded" in col_map:
                break

    if "unstranded" not in col_map:
        raise ColumnDetectionError(str(filepath), cols_lower)

    return col_map


def _read_lines(filepath: Path, max_lines: int | None = None) -> list[str]:
    if filepath.suffix == ".gz":
        with gzip.open(filepath, "rt", encoding="utf-8") as f:
            return [f.readline() for _ in range(max_lines)] if max_lines else f.readlines()
    with filepath.open("r", encoding="utf-8") as f:
        return [f.readline() for _ in range(max_lines)] if max_lines else f.readlines()


# =============================================================================
#   Single-file parser
# =============================================================================

def parse_star_file(filepath: Path, unstranded_col_idx: int) -> pd.Series | None:
    """
    Parse one STAR counts TSV. Returns Series (gene_id → count) or None.

    Excludes: N_ summary rows, header row, any row where count is not integer.
    """
    try:
        lines = _read_lines(filepath)
    except Exception as e:
        print(f"    ⚠️  Cannot read {filepath.name}: {e}")
        return None

    counts: dict[str, int] = {}
    header_seen = False

    for line in lines:
        line = line.strip()
        if not line or line.startswith("N_") or line.startswith("#"):
            continue
        parts = line.split("\t")
        gene_id = parts[0].strip()

        # Skip header row
        if gene_id.lower() in ("gene_id", "#gene_id", "geneid") and not header_seen:
            header_seen = True
            continue

        if len(parts) <= unstranded_col_idx:
            continue

        try:
            counts[gene_id] = int(parts[unstranded_col_idx].strip())
        except ValueError:
            if not header_seen:
                header_seen = True
            continue

    if not counts:
        print(f"    ⚠️  No counts parsed from {filepath.name}")
        return None

    return pd.Series(counts, dtype="int64", name=filepath.name)


# =============================================================================
#   Matrix builder
# =============================================================================

def build_count_matrix(
    counts_dir: Path,
    metadata_df: pd.DataFrame,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Build a gene × sample count matrix from extracted STAR TSV files.

    File-to-sample matching uses this priority:
      1. UUID directory name (parent dir of TSV = file_id from GDC ZIP)
      2. Filename match
      3. UNMATCHED_ fallback (file is parsed, flagged in output)

    Parameters
    ----------
    counts_dir : Path
        Directory containing extracted GDC ZIP contents.
    metadata_df : pd.DataFrame
        Must contain 'file_id', 'file_name', 'sample_submitter_id'.
    verbose : bool

    Returns
    -------
    pd.DataFrame
        Shape: (n_genes, n_samples). Index = gene_id.
    """
    tsv_files = list(counts_dir.rglob("*.tsv")) + list(counts_dir.rglob("*.tsv.gz"))

    if not tsv_files:
        raise GDCError(
            f"No TSV files found under {counts_dir}.",
            fix=(
                "The download may have failed or produced an empty archive.\n"
                "Delete the 'downloaded' checkpoint entry and re-run to retry."
            ),
            step="count matrix build",
        )

    if verbose:
        print(f"  🔍  Found {len(tsv_files)} TSV files in {counts_dir.name}/")

    # ── Detect column layout from first parseable file ────────────────────────
    col_map: dict[str, int] = {}
    reference_file: Path | None = None

    for fp in tsv_files:
        try:
            col_map = detect_column_map(fp)
            reference_file = fp
            break
        except (ColumnDetectionError, GDCError):
            continue

    if not col_map or "unstranded" not in col_map:
        raise GDCError(
            "Could not detect the unstranded count column in any downloaded file.",
            fix=(
                "Open one of the TSV files in a text editor to inspect its columns.\n"
                "Expected: gene_id  unstranded  stranded_first  stranded_second\n"
                "If the format differs, report at the project GitHub issues page."
            ),
            step="column detection",
        )

    unstranded_idx = col_map["unstranded"]
    if verbose:
        print(f"  ✅  Column map (from {reference_file.name}):")
        for k, v in col_map.items():
            print(f"        {k:<20s} → column index {v}")

    # ── Build TWO lookup tables: UUID→sample_id and filename→sample_id ────────
    # Primary:  file_id UUID   (= directory name in GDC ZIP extraction)
    # Fallback: file_name      (actual TSV filename)
    uuid_to_sample: dict[str, str] = {}
    name_to_sample: dict[str, str] = {}

    for _, row in metadata_df.iterrows():
        fid  = str(row.get("file_id",   "")).strip()
        fname = str(row.get("file_name", "")).strip()
        sid   = str(row.get("sample_submitter_id", "")).strip()
        if fid and sid:
            uuid_to_sample[fid] = sid
        if fname and sid:
            name_to_sample[fname] = sid

    # ── Parse all files ───────────────────────────────────────────────────────
    count_series: dict[str, pd.Series] = {}
    failed: list[str] = []
    unmatched: list[str] = []
    match_method_counts = {"uuid": 0, "filename": 0, "fallback": 0}

    iter_files = tqdm(tsv_files, desc="  Parsing files", ncols=70, unit="file") \
        if verbose else tsv_files

    for fp in iter_files:
        # Priority 1: match by parent directory name (= file_id UUID from GDC ZIP)
        parent_name = fp.parent.name
        sample_id = None
        if is_valid_gdc_uuid(parent_name):
            sample_id = uuid_to_sample.get(parent_name)
            if sample_id:
                match_method_counts["uuid"] += 1

        # Priority 2: match by filename
        if sample_id is None:
            sample_id = name_to_sample.get(fp.name)
            if sample_id:
                match_method_counts["filename"] += 1

        # Priority 3: fallback — use UUID dir name or stem as column label
        if sample_id is None:
            label = parent_name if is_valid_gdc_uuid(parent_name) else fp.stem[:16]
            sample_id = f"UNMATCHED_{label[:8]}"
            unmatched.append(fp.name)
            match_method_counts["fallback"] += 1

        series = parse_star_file(fp, unstranded_idx)
        if series is None:
            failed.append(fp.name)
            continue

        # Deduplicate sample IDs (multiple aliquots from same sample)
        base, n = sample_id, 1
        while sample_id in count_series:
            sample_id = f"{base}_dup{n}"
            n += 1
        count_series[sample_id] = series

    # ── Matching report ───────────────────────────────────────────────────────
    if verbose:
        print(f"\n  🔗  File matching report:")
        print(f"        Matched by UUID directory: {match_method_counts['uuid']}")
        print(f"        Matched by filename:        {match_method_counts['filename']}")
        print(f"        Unmatched (fallback ID):    {match_method_counts['fallback']}")
        if unmatched:
            print(f"\n  ⚠️  {len(unmatched)} unmatched file(s) — "
                  f"sample IDs begin with UNMATCHED_ in output:")
            for name in unmatched[:5]:
                print(f"       • {name}")
            if len(unmatched) > 5:
                print(f"       ... and {len(unmatched)-5} more")

    if failed:
        if verbose:
            print(f"\n  ⚠️  {len(failed)} file(s) could not be parsed:")
            for name in failed[:10]:
                print(f"       • {name}")

    if not count_series:
        raise GDCError(
            "No count files could be successfully parsed.",
            fix=(
                "Try deleting the raw_counts directory and re-running.\n"
                "If the problem persists, open one TSV file to inspect its format."
            ),
            step="count matrix build",
        )

    if verbose:
        print(f"\n  🔧  Assembling matrix ({len(count_series)} samples)...")

    matrix = pd.DataFrame(count_series)
    matrix.index.name = "gene_id"

    if verbose:
        _report_matrix_qc(matrix)

    return matrix


def _report_matrix_qc(matrix: pd.DataFrame) -> None:
    n_genes, n_samples = matrix.shape
    zero_samples = (matrix.sum(axis=0) == 0).sum()
    zero_genes   = (matrix.sum(axis=1) == 0).sum()
    median_lib   = matrix.sum(axis=0).median()

    print(f"\n  📊  Matrix QC:")
    print(f"      Genes:                          {n_genes}")
    print(f"      Samples:                        {n_samples}")
    print(f"      Median library size (counts):  {median_lib:,.0f}")
    print(f"      Genes with zero counts in all samples: {zero_genes}")
    if zero_samples:
        print(f"      ⚠️  Samples with ALL-ZERO counts: {zero_samples}")
        print("         (Possible failed library preps — review before analysis)")
    else:
        print("      All samples have non-zero counts ✅")
