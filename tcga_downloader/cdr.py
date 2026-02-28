"""
cdr.py — TCGA PanCanAtlas Clinical Data Resource (CDR) integration
==================================================================

SOURCE
------
Liu J et al. "An Integrated TCGA Pan-Cancer Clinical Data Resource to
Drive High-Quality Survival Outcome Analytics." Cell. 2018.
https://doi.org/10.1016/j.cell.2018.02.052

The CDR Excel file is hosted on the GDC with a stable UUID.
It covers all 33 TCGA cancer types with harmonised clinical annotations,
four curated survival endpoints, and molecular subtype calls.

SCOPE GATE
----------
CDR integration is only available for TCGA projects.
is_tcga_project() must be called before any other function here.
Non-TCGA projects raise a clear ProgramNotSupportedError immediately.

CDR COVERAGE REALITY
--------------------
Not all patients in a TCGA GDC download will have a CDR entry.
Three situations are handled explicitly:

  1. Case matched CDR, all key fields populated  → cdr_matched=True
  2. Case matched CDR, some fields missing        → cdr_matched=True,
                                                    per-field NaN
  3. Case NOT in CDR (post-2018 or excluded)      → cdr_matched=False,
                                                    all CDR cols = NaN

The caller is never silently given empty data — every unmatched case
is flagged, reported, and written to a separate file.

OUTPUT COLUMNS ADDED
--------------------
All CDR columns are prefixed with 'cdr_' in the output to distinguish
them from GDC API clinical fields that cover the same concept.
Three audit flag columns are also added:

  cdr_matched            bool  — did this case appear in the CDR at all
  cdr_subtype_available  bool  — is Subtype_Selected (or Subtype_mRNA) populated
  cdr_survival_complete  bool  — are OS, OS.time, PFI, PFI.time all present
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tcga_downloader.constants import (
    CDR_FILE_UUID,
    CDR_CACHE_FILENAME,
    CDR_JOIN_KEY,
    CDR_SURVIVAL_COLS,
    CDR_SUBTYPE_COL,
    CDR_SUBTYPE_COL_LEGACY,
    CDR_COLUMNS_TO_KEEP,
    CDR_SUPPORTED_PROGRAMS,
    SUPPORTED_PROGRAMS,
)
from tcga_downloader.exceptions import GDCError


# =============================================================================
#   Program gate
# =============================================================================

class ProgramNotSupportedError(GDCError):
    """Raised when the user attempts to use a non-TCGA project."""
    def __init__(self, program_name: str, project_id: str):
        super().__init__(
            message=(
                f"Project '{project_id}' belongs to program '{program_name}', not TCGA.\n"
                f"This tool currently supports TCGA projects only."
            ),
            fix=(
                "Browse supported TCGA projects at:\n"
                "  https://portal.gdc.cancer.gov/projects\n"
                "Filter by program = TCGA. Examples: TCGA-BRCA, TCGA-GBM, TCGA-LUAD."
            ),
            step="program validation",
        )


def is_tcga_project(program_name: str, project_id: str = "") -> bool:
    """
    Return True if the program is TCGA. Raise ProgramNotSupportedError if not.

    Parameters
    ----------
    program_name : str
        The 'program.name' field from the GDC /projects API response.
    project_id : str
        Project ID for the error message (optional, improves message quality).

    Raises
    ------
    ProgramNotSupportedError
        If program_name is not in SUPPORTED_PROGRAMS.
    """
    if program_name in SUPPORTED_PROGRAMS:
        return True
    raise ProgramNotSupportedError(program_name, project_id)


def is_cdr_supported(program_name: str) -> bool:
    """Return True if CDR integration is available for this program."""
    return program_name in CDR_SUPPORTED_PROGRAMS


# =============================================================================
#   CDR download and caching
# =============================================================================

def download_cdr_file(client, cache_dir: Path, verbose: bool = True) -> Path:
    """
    Download the CDR Excel file from the GDC and cache it locally.

    Uses the stable GDC UUID for the CDR file. If the file is already
    cached in cache_dir, the download is skipped entirely.

    Parameters
    ----------
    client : GDCClient
        Authenticated (or open-access) GDC client instance.
    cache_dir : Path
        Directory where the CDR Excel file will be cached.
    verbose : bool

    Returns
    -------
    Path
        Path to the cached CDR Excel file.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / CDR_CACHE_FILENAME

    if dest.exists() and dest.stat().st_size > 100_000:
        if verbose:
            print(f"  ✅  CDR file already cached: {dest.name}")
        return dest

    if verbose:
        print(f"  📥  Downloading PanCanAtlas CDR file (one-time, ~2 MB)...")

    # CDR file is open-access — no token required
    n_bytes = client.stream_download([CDR_FILE_UUID], str(dest))

    if dest.stat().st_size < 100_000:
        raise GDCError(
            "CDR file downloaded but appears too small — may be corrupted.",
            fix=(
                f"Delete {dest} and re-run to retry.\n"
                "If the problem persists, the GDC UUID may have changed.\n"
                "Check: https://portal.gdc.cancer.gov/files/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
            ),
            step="CDR download",
        )

    if verbose:
        print(f"  ✅  CDR downloaded: {n_bytes / (1024**2):.1f} MB → {dest.name}")

    return dest


# =============================================================================
#   CDR parsing
# =============================================================================

def parse_cdr(filepath: Path, verbose: bool = True) -> pd.DataFrame:
    """
    Parse the CDR Excel file into a clean DataFrame.

    Handles known CDR quirks:
      - '#N/A' strings → NaN
      - '[Not Available]', '[Not Applicable]', '[Unknown]' → NaN
      - Strips whitespace from string columns
      - Selects only the columns defined in CDR_COLUMNS_TO_KEEP
      - Prefixes all columns (except join key) with 'cdr_'

    Parameters
    ----------
    filepath : Path
        Path to the CDR Excel file.

    Returns
    -------
    pd.DataFrame
        One row per TCGA case. Join key column is 'bcr_patient_barcode'
        (unprefixed). All other columns are prefixed 'cdr_'.
    """
    try:
        raw = pd.read_excel(
            filepath,
            sheet_name=0,
            dtype=str,           # read everything as string first
            na_values=[
                "#N/A", "N/A", "NA", "",
                "[Not Available]", "[Not Applicable]",
                "[Unknown]", "[Discrepancy]", "NaN",
            ],
            keep_default_na=True,
        )
    except Exception as e:
        raise GDCError(
            f"Cannot parse CDR Excel file: {filepath.name}",
            fix=(
                f"Delete {filepath} and re-run to download a fresh copy.\n"
                "If this persists, check that openpyxl is installed:\n"
                "  pip install openpyxl"
            ),
            step="CDR parsing",
        ) from e

    # ── Keep only the columns we need ────────────────────────────────────────
    available = [c for c in CDR_COLUMNS_TO_KEEP if c in raw.columns]
    missing   = [c for c in CDR_COLUMNS_TO_KEEP if c not in raw.columns]

    if CDR_JOIN_KEY not in raw.columns:
        raise GDCError(
            f"CDR file is missing the join key column '{CDR_JOIN_KEY}'.",
            fix=(
                "The CDR file may be a different version than expected.\n"
                "Delete the cached file and re-run to download a fresh copy."
            ),
            step="CDR parsing",
        )

    if missing and verbose:
        print(f"  ℹ️   {len(missing)} CDR columns not found in file "
              f"(may be in a different CDR version): {missing[:5]}")

    df = raw[available].copy()

    # ── Clean string values ───────────────────────────────────────────────────
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip()

    # ── Convert numeric-looking columns to float ──────────────────────────────
    numeric_hints = {*CDR_SURVIVAL_COLS,
                     "age_at_initial_pathologic_diagnosis",
                     "birth_days_to", "last_contact_days_to",
                     "death_days_to", "new_tumor_event_dx_days_to",
                     "initial_pathologic_dx_year"}
    for col in df.columns:
        if col in numeric_hints:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # ── Prefix all columns except the join key ────────────────────────────────
    rename = {c: f"cdr_{c}" for c in df.columns if c != CDR_JOIN_KEY}
    df = df.rename(columns=rename)

    if verbose:
        print(f"  ✅  CDR parsed: {len(df)} cases, {len(df.columns)} columns")

    return df


# =============================================================================
#   CDR join
# =============================================================================

def join_cdr(
    metadata_df: pd.DataFrame,
    cdr_df: pd.DataFrame,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Left-join CDR annotations onto the sample metadata DataFrame.

    Every sample row is preserved. Samples whose case_submitter_id does
    not appear in the CDR get NaN for all CDR columns.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Must contain 'case_submitter_id' column.
    cdr_df : pd.DataFrame
        Output of parse_cdr() — must contain 'bcr_patient_barcode' column.

    Returns
    -------
    pd.DataFrame
        metadata_df with CDR columns appended. Does not modify inputs.
    """
    if "case_submitter_id" not in metadata_df.columns:
        raise GDCError(
            "metadata_df is missing 'case_submitter_id' column.",
            fix="This is an internal error — please report it.",
            step="CDR join",
        )

    merged = metadata_df.merge(
        cdr_df,
        left_on="case_submitter_id",
        right_on=CDR_JOIN_KEY,
        how="left",
        suffixes=("", "_cdr_dup"),
    )

    # Drop the redundant join key column from CDR side
    if CDR_JOIN_KEY in merged.columns:
        merged = merged.drop(columns=[CDR_JOIN_KEY])

    # Drop any accidental duplicate columns from the merge
    dup_cols = [c for c in merged.columns if c.endswith("_cdr_dup")]
    if dup_cols:
        merged = merged.drop(columns=dup_cols)

    return merged


# =============================================================================
#   Audit flags
# =============================================================================

def audit_cdr_join(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add three boolean audit flag columns to the merged DataFrame.

    Columns added (inserted near the front, after sample identifiers):
      cdr_matched            — True if case appeared in CDR at all
      cdr_subtype_available  — True if Subtype_Selected or Subtype_mRNA is populated
      cdr_survival_complete  — True if OS, OS.time, PFI, PFI.time all present

    The minimum survival completeness threshold is OS + PFI because:
      - OS (overall survival) is the most universally used endpoint
      - PFI (progression-free interval) is the most informative for
        event-driven analyses
      - DSS and DFI have lower coverage and are not required for the flag

    Parameters
    ----------
    merged_df : pd.DataFrame
        Output of join_cdr().

    Returns
    -------
    pd.DataFrame
        Same DataFrame with three new columns added. Modifies in place
        and returns for chaining.
    """
    # cdr_matched: True if ANY cdr_ column (other than flags) is non-null
    cdr_data_cols = [c for c in merged_df.columns
                     if c.startswith("cdr_") and c not in
                     ("cdr_matched", "cdr_subtype_available", "cdr_survival_complete")]

    if cdr_data_cols:
        merged_df["cdr_matched"] = merged_df[cdr_data_cols].notna().any(axis=1)
    else:
        merged_df["cdr_matched"] = False

    # cdr_subtype_available
    # Check primary column (Subtype_Selected in GDC-hosted CDR) first,
    # then fall back to legacy column (Subtype_mRNA in original paper supplement)
    subtype_col = f"cdr_{CDR_SUBTYPE_COL}"
    subtype_col_legacy = f"cdr_{CDR_SUBTYPE_COL_LEGACY}"
    if subtype_col in merged_df.columns:
        merged_df["cdr_subtype_available"] = merged_df[subtype_col].notna()
    elif subtype_col_legacy in merged_df.columns:
        merged_df["cdr_subtype_available"] = merged_df[subtype_col_legacy].notna()
    else:
        merged_df["cdr_subtype_available"] = False

    # cdr_survival_complete: OS + OS.time + PFI + PFI.time all non-null
    min_survival_cols = ["cdr_OS", "cdr_OS.time", "cdr_PFI", "cdr_PFI.time"]
    present_survival = [c for c in min_survival_cols if c in merged_df.columns]
    if present_survival:
        merged_df["cdr_survival_complete"] = (
            merged_df[present_survival].notna().all(axis=1)
        )
    else:
        merged_df["cdr_survival_complete"] = False

    return merged_df


# =============================================================================
#   Coverage report
# =============================================================================

def generate_coverage_report(
    merged_df: pd.DataFrame,
    project_id: str,
    output_dir: Path,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compute and save a per-field CDR coverage report.

    For every CDR column, reports how many samples have a non-null value
    and what percentage of total samples that represents.

    Also reports the three audit flag summaries and lists unmatched
    case barcodes to a separate text file if any exist.

    Parameters
    ----------
    merged_df : pd.DataFrame
        Output of audit_cdr_join().
    project_id : str
        Used for output filenames.
    output_dir : Path
    verbose : bool

    Returns
    -------
    pd.DataFrame
        Coverage report table (also saved to disk).
    """
    n_total = len(merged_df)
    safe_pid = project_id.replace("/", "_")

    cdr_cols = [c for c in merged_df.columns if c.startswith("cdr_")
                and c not in ("cdr_matched", "cdr_subtype_available",
                               "cdr_survival_complete")]

    rows = []
    for col in cdr_cols:
        n_present = merged_df[col].notna().sum()
        rows.append({
            "field":       col,
            "n_present":   int(n_present),
            "n_total":     n_total,
            "pct_present": round(100 * n_present / n_total, 1) if n_total else 0.0,
        })

    report_df = pd.DataFrame(rows).sort_values("pct_present", ascending=False)

    # ── Save coverage report ──────────────────────────────────────────────────
    report_path = output_dir / f"{safe_pid}_CDR_coverage_report.tsv"
    report_df.to_csv(report_path, sep="\t", index=False)

    # ── Save unmatched cases ──────────────────────────────────────────────────
    unmatched_path = None
    if "cdr_matched" in merged_df.columns:
        unmatched = merged_df.loc[
            ~merged_df["cdr_matched"], "case_submitter_id"
        ].dropna().unique()

        if len(unmatched) > 0:
            unmatched_path = output_dir / f"{safe_pid}_CDR_unmatched_cases.txt"
            unmatched_path.write_text(
                "\n".join(sorted(unmatched)) + "\n", encoding="utf-8"
            )

    # ── Print summary ─────────────────────────────────────────────────────────
    if verbose:
        n_matched   = int(merged_df["cdr_matched"].sum()) if "cdr_matched" in merged_df.columns else 0
        n_subtype   = int(merged_df["cdr_subtype_available"].sum()) if "cdr_subtype_available" in merged_df.columns else 0
        n_survival  = int(merged_df["cdr_survival_complete"].sum()) if "cdr_survival_complete" in merged_df.columns else 0
        n_unmatched = n_total - n_matched

        print(f"\n  📊  CDR coverage report for {project_id} ({n_total} samples):")
        print(f"      cdr_matched:           {n_matched:>5} / {n_total}  "
              f"({100*n_matched/n_total:.1f}%)")
        print(f"      cdr_subtype_available: {n_subtype:>5} / {n_total}  "
              f"({100*n_subtype/n_total:.1f}%)")
        print(f"      cdr_survival_complete: {n_survival:>5} / {n_total}  "
              f"({100*n_survival/n_total:.1f}%)")

        if n_unmatched > 0:
            print(f"\n      ⚠️  {n_unmatched} case(s) not matched to CDR.")
            print("         Likely causes:")
            print("           • Sample added to GDC after the 2018 TCGA data freeze")
            print("           • Case excluded from CDR due to QC failure")
            print(f"         Unmatched barcodes saved to: "
                  f"{unmatched_path.name if unmatched_path else 'N/A'}")

        print(f"\n      Top CDR fields by coverage:")
        for _, row in report_df.head(10).iterrows():
            bar = "█" * int(row["pct_present"] / 5)
            print(f"        {row['field']:<35} {row['n_present']:>5}  "
                  f"({row['pct_present']:5.1f}%)  {bar}")

        print(f"\n      Full report saved to: {report_path.name}")

    return report_df


# =============================================================================
#   Complete-cases split
# =============================================================================

def split_complete_cases(
    merged_df: pd.DataFrame,
    project_id: str,
    output_dir: Path,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split merged DataFrame into complete-CDR and incomplete-CDR subsets.

    'Complete' is defined as cdr_survival_complete == True.
    If the project has subtype data available for >50% of matched cases,
    cdr_subtype_available is also required for the complete subset.

    Both subsets are written to disk.

    Parameters
    ----------
    merged_df : pd.DataFrame
        Output of audit_cdr_join().
    project_id : str
    output_dir : Path
    verbose : bool

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        (complete_df, incomplete_df)
    """
    safe_pid = project_id.replace("/", "_")

    survival_complete = merged_df.get("cdr_survival_complete",
                                      pd.Series([False] * len(merged_df)))

    # Decide whether to require subtype in the complete-cases definition
    if "cdr_subtype_available" in merged_df.columns and "cdr_matched" in merged_df.columns:
        n_matched  = merged_df["cdr_matched"].sum()
        n_subtype  = merged_df["cdr_subtype_available"].sum()
        subtype_coverage = (n_subtype / n_matched) if n_matched > 0 else 0.0
        require_subtype  = subtype_coverage >= 0.5
    else:
        require_subtype = False

    if require_subtype:
        subtype_ok = merged_df.get("cdr_subtype_available",
                                    pd.Series([False] * len(merged_df)))
        complete_mask = survival_complete & subtype_ok
        completeness_criteria = "cdr_survival_complete AND cdr_subtype_available"
    else:
        complete_mask = survival_complete
        completeness_criteria = "cdr_survival_complete only (subtype coverage <50%)"

    complete_df   = merged_df[complete_mask].copy()
    incomplete_df = merged_df[~complete_mask].copy()

    # ── Write files ───────────────────────────────────────────────────────────
    complete_path = output_dir / f"{safe_pid}_FULL_merged_CDR_complete_cases.tsv"
    complete_df.to_csv(complete_path, sep="\t", index=False)

    if verbose:
        print(f"\n  ✂️   Complete-cases split ({completeness_criteria}):")
        print(f"      Complete:   {len(complete_df):>5} samples → {complete_path.name}")
        print(f"      Incomplete: {len(incomplete_df):>5} samples "
              f"(kept in ALL and TUMOR/NORMAL files with NaN for missing CDR fields)")

    return complete_df, incomplete_df


# =============================================================================
#   Top-level orchestrator
# =============================================================================

def run_cdr_pipeline(
    client,
    metadata_df: pd.DataFrame,
    project_id: str,
    program_name: str,
    output_dir: Path,
    verbose: bool = True,
) -> dict:
    """
    Run the full CDR pipeline for a TCGA project.

    Steps:
      1. Gate check — abort immediately if not TCGA
      2. Download CDR file (or use cache)
      3. Parse CDR Excel
      4. Join to metadata
      5. Add audit flags
      6. Generate coverage report
      7. Split complete vs incomplete cases
      8. Save CDR-only annotation file

    Parameters
    ----------
    client : GDCClient
    metadata_df : pd.DataFrame
        Must contain 'case_submitter_id'.
    project_id : str
    program_name : str
        From project info — used for TCGA gate check.
    output_dir : Path
    verbose : bool

    Returns
    -------
    dict
        Keys: 'merged_df', 'complete_df', 'incomplete_df',
              'coverage_report_df', 'cdr_df', 'output_paths'
    """
    # ── Step 1: Gate check ────────────────────────────────────────────────────
    is_tcga_project(program_name, project_id)

    safe_pid   = project_id.replace("/", "_")
    cache_dir  = output_dir / "cdr_cache"

    # ── Step 2: Download ──────────────────────────────────────────────────────
    cdr_path = download_cdr_file(client, cache_dir, verbose=verbose)

    # ── Step 3: Parse ─────────────────────────────────────────────────────────
    cdr_df = parse_cdr(cdr_path, verbose=verbose)

    # ── Step 4: Join ──────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  🔗  Joining CDR to {len(metadata_df)} samples...")
    merged_df = join_cdr(metadata_df, cdr_df, verbose=verbose)

    # ── Step 5: Audit flags ───────────────────────────────────────────────────
    merged_df = audit_cdr_join(merged_df)

    # ── Step 6: Coverage report ───────────────────────────────────────────────
    coverage_df = generate_coverage_report(
        merged_df, project_id, output_dir, verbose=verbose
    )

    # ── Step 7: Complete-cases split ──────────────────────────────────────────
    complete_df, incomplete_df = split_complete_cases(
        merged_df, project_id, output_dir, verbose=verbose
    )

    # ── Step 8: Save CDR-only annotation file ─────────────────────────────────
    cdr_only_cols = (["case_submitter_id", "sample_submitter_id",
                      "cdr_matched", "cdr_subtype_available", "cdr_survival_complete"]
                     + [c for c in merged_df.columns if c.startswith("cdr_")
                        and c not in ("cdr_matched","cdr_subtype_available",
                                      "cdr_survival_complete")])
    cdr_only_cols = [c for c in cdr_only_cols if c in merged_df.columns]
    cdr_annot_path = output_dir / f"{safe_pid}_CDR_annotations.tsv"
    merged_df[cdr_only_cols].to_csv(cdr_annot_path, sep="\t", index=False)

    if verbose:
        print(f"  💾  CDR annotations saved: {cdr_annot_path.name}")

    output_paths = {
        "cdr_annotations":  cdr_annot_path,
        "coverage_report":  output_dir / f"{safe_pid}_CDR_coverage_report.tsv",
        "complete_cases":   output_dir / f"{safe_pid}_FULL_merged_CDR_complete_cases.tsv",
    }
    unmatched_path = output_dir / f"{safe_pid}_CDR_unmatched_cases.txt"
    if unmatched_path.exists():
        output_paths["unmatched_cases"] = unmatched_path

    return {
        "merged_df":       merged_df,
        "complete_df":     complete_df,
        "incomplete_df":   incomplete_df,
        "coverage_report_df": coverage_df,
        "cdr_df":          cdr_df,
        "output_paths":    output_paths,
    }
