"""
merge.py
========
Merges the count matrix, file metadata, and clinical data into
final output TSV files.

Output files produced
---------------------
- {project}_STAR_unstranded_merged_ALL.tsv    — all samples
- {project}_STAR_unstranded_TUMOR_ONLY.tsv    — tumor samples only
- {project}_STAR_unstranded_NORMAL_ONLY.tsv   — normal samples only
- {project}_sample_metadata_clinical.tsv      — metadata + clinical, no counts
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

from tcga_downloader.constants import TUMOR_CODES, NORMAL_CODES, SAMPLE_TYPE_CODES
from tcga_downloader.exceptions import GDCError


# Column ordering: these appear first (in this order) in every output file
PRIORITY_COLUMNS = [
    "sample_submitter_id",
    "case_submitter_id",
    "sample_id",
    "case_id",
    "file_id",
    "file_name",
    "sample_category",       # Tumor / Normal / Unknown
    "sample_type_label",     # Primary Tumor / Solid Tissue Normal / ...
    "sample_type_code",      # 01 / 11 / ...
    "gdc_sample_type",       # GDC's own label (may differ from barcode-derived)
    "tissue_type",
    # clinical demographics
    "gender",
    "age_at_index",
    "age_at_diagnosis",
    "race",
    "ethnicity",
    "vital_status",
    "days_to_death",
    "days_to_last_follow_up",
    # tumour characteristics
    "primary_diagnosis",
    "tumor_stage",
    "tumor_grade",
    "ajcc_pathologic_stage",
    "ajcc_pathologic_t",
    "ajcc_pathologic_n",
    "ajcc_pathologic_m",
    "morphology",
    "tissue_or_organ_of_origin",
    "site_of_resection_or_biopsy",
    "prior_malignancy",
    "synchronous_malignancy",
    # exposures
    "cigarettes_per_day",
    "years_smoked",
    "pack_years_smoked",
    "alcohol_history",
    "bmi",
    "height",
    "weight",
]


def flatten_clinical(case_hits: list[dict]) -> pd.DataFrame:
    """
    Flatten nested GDC case JSON records into a flat DataFrame.

    GDC returns nested sub-documents for demographic, diagnoses, exposures.
    We take the first element of each list-type field (most common case).

    Parameters
    ----------
    case_hits : list[dict]
        Raw 'hits' list from the GDC /cases API response.

    Returns
    -------
    pd.DataFrame
        One row per case, one column per clinical variable.
    """
    if not case_hits:
        return pd.DataFrame()

    rows = []
    for case in case_hits:
        row: dict = {
            "case_id":           case.get("case_id", ""),
            "case_submitter_id": case.get("submitter_id", ""),
        }

        # Demographic (single dict, not a list)
        demo = case.get("demographic") or {}
        for k in ("gender", "age_at_index", "race", "ethnicity",
                  "vital_status", "days_to_death", "year_of_birth"):
            row[k] = demo.get(k, "")

        # Diagnoses (list — take first element)
        diag = (case.get("diagnoses") or [{}])[0]
        for k in ("age_at_diagnosis", "days_to_last_follow_up",
                  "primary_diagnosis", "tumor_stage", "tumor_grade",
                  "morphology", "tissue_or_organ_of_origin",
                  "site_of_resection_or_biopsy", "prior_malignancy",
                  "synchronous_malignancy", "ajcc_pathologic_stage",
                  "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m"):
            row[k] = diag.get(k, "")

        # Exposures (list — take first element)
        expo = (case.get("exposures") or [{}])[0]
        for k in ("cigarettes_per_day", "years_smoked", "pack_years_smoked",
                  "alcohol_history", "bmi", "height", "weight"):
            row[k] = expo.get(k, "")

        rows.append(row)

    df = pd.DataFrame(rows)
    # Safety net: deduplicate any columns that appear more than once
    df = df.loc[:, ~df.columns.duplicated()]
    return df


def flatten_file_metadata(file_hits: list[dict]) -> pd.DataFrame:
    """
    Flatten GDC file-level metadata into a DataFrame.

    Parameters
    ----------
    file_hits : list[dict]
        Raw 'hits' list from the GDC /files API response.

    Returns
    -------
    pd.DataFrame
        One row per file, with file_id and sample identifiers as columns.
    """
    from tcga_downloader.sample import annotate_metadata

    rows = []
    for fhit in file_hits:
        case  = (fhit.get("cases") or [{}])[0]
        sample = (case.get("samples") or [{}])[0]
        rows.append({
            "file_id":              fhit.get("file_id", ""),
            "file_name":            fhit.get("file_name", ""),
            "file_size_bytes":      fhit.get("file_size", 0),
            "case_id":              case.get("case_id", ""),
            "case_submitter_id":    case.get("submitter_id", ""),
            "sample_id":            sample.get("sample_id", ""),
            "sample_submitter_id":  sample.get("submitter_id", ""),
            "gdc_sample_type":      sample.get("sample_type", ""),
            "tissue_type":          sample.get("tissue_type", ""),
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        annotate_metadata(df)   # adds sample_category, sample_type_label, sample_type_code
    return df


def merge_outputs(
    count_matrix: pd.DataFrame,
    metadata_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    project_id: str,
    output_dir: Path,
    verbose: bool = True,
) -> dict[str, Path]:
    """
    Produce all final output TSV files.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Gene × sample matrix (from matrix.build_count_matrix).
    metadata_df : pd.DataFrame
        File-level metadata (from flatten_file_metadata).
    clinical_df : pd.DataFrame
        Clinical data (from flatten_clinical).
    project_id : str
    output_dir : Path
    verbose : bool

    Returns
    -------
    dict[str, Path]
        Keys: 'all', 'tumor', 'normal', 'metadata'.
        Values: Path objects of written files (only keys with actual output
        are included).
    """
    if verbose:
        print("  🔗  Joining counts + metadata + clinical data...")

    # Transpose: samples become rows
    counts_T = (
        count_matrix.T
        .reset_index()
        .rename(columns={"index": "sample_submitter_id"})
    )

    # Left-join metadata
    merged = counts_T.merge(metadata_df, on="sample_submitter_id", how="left")
    # Deduplicate any columns that appear twice after the metadata join
    merged = merged.loc[:, ~merged.columns.duplicated()]

    # Left-join clinical — avoid duplicating already-present columns
    if not clinical_df.empty:
        # Deduplicate clinical_df columns first — GDC JSON sometimes produces
        # duplicate column names (e.g. case_submitter_id appearing twice)
        clinical_df = clinical_df.loc[:, ~clinical_df.columns.duplicated()]

        existing = set(merged.columns)
        clin_cols = ["case_submitter_id"] + [
            c for c in clinical_df.columns if c not in existing and c != "case_submitter_id"
        ]
        merged = merged.merge(clinical_df[clin_cols], on="case_submitter_id", how="left")

    # Re-order: priority metadata/clinical first, then genes
    gene_cols = list(count_matrix.index)  # all gene IDs, in original order
    present_priority = [c for c in PRIORITY_COLUMNS if c in merged.columns]
    other_meta = [c for c in merged.columns
                  if c not in present_priority and c not in gene_cols]
    ordered = merged[present_priority + other_meta + gene_cols]

    # ── Validation ────────────────────────────────────────────────────────────
    if verbose:
        _print_validation_report(ordered, clinical_df)

    # ── Write files ───────────────────────────────────────────────────────────
    safe_pid = project_id.replace("/", "_")
    output_paths: dict[str, Path] = {}

    def _save(df: pd.DataFrame, suffix: str) -> Path:
        path = output_dir / f"{safe_pid}_{suffix}.tsv"
        df.to_csv(path, sep="\t", index=False)
        size_mb = path.stat().st_size / (1024 ** 2)
        if verbose:
            print(f"  💾  {path.name}  ({len(df)} samples, {size_mb:.1f} MB)")
        return path

    output_paths["all"]      = _save(ordered, "STAR_unstranded_merged_ALL")

    if "sample_category" in ordered.columns:
        tumor_df  = ordered[ordered["sample_category"] == "Tumor"]
        normal_df = ordered[ordered["sample_category"] == "Normal"]
        if not tumor_df.empty:
            output_paths["tumor"]  = _save(tumor_df,  "STAR_unstranded_TUMOR_ONLY")
        if not normal_df.empty:
            output_paths["normal"] = _save(normal_df, "STAR_unstranded_NORMAL_ONLY")

    meta_only = ordered[present_priority + other_meta]
    output_paths["metadata"] = _save(meta_only, "sample_metadata_clinical")

    return output_paths


def _print_validation_report(merged: pd.DataFrame, clinical_df: pd.DataFrame) -> None:
    """Print a plain-English data integrity report before saving."""
    print("\n  🔎  Validation report:")

    n = len(merged)
    print(f"      Total samples: {n}")

    if "sample_category" in merged.columns:
        cats = merged["sample_category"].value_counts()
        print(f"      Tumor samples:   {cats.get('Tumor', 0)}")
        print(f"      Normal samples:  {cats.get('Normal', 0)}")
        print(f"      Unknown samples: {cats.get('Unknown', 0)}")

    if not clinical_df.empty and "case_submitter_id" in merged.columns:
        unmatched = merged["case_submitter_id"].isna().sum()
        if unmatched:
            print(f"      ⚠️  {unmatched} samples have no matched clinical data")
            print("          (This can happen if clinical records are not yet available in GDC)")
        else:
            print("      All samples matched to clinical data ✅")

    if "sample_submitter_id" in merged.columns:
        dupes = merged["sample_submitter_id"].duplicated().sum()
        if dupes:
            print(f"      ⚠️  {dupes} duplicate sample_submitter_id values")
            print("          (Multiple aliquots from the same sample — check before DE analysis)")


def save_full_merged_with_cdr(
    count_matrix,
    metadata_with_cdr: "pd.DataFrame",
    project_id: str,
    output_dir: "Path",
    verbose: bool = True,
) -> dict:
    """
    Produce the final FULL merged files that include CDR annotations.

    Called after run_cdr_pipeline() has added CDR columns to metadata_df.
    Produces:
        {project}_FULL_merged_with_CDR.tsv   — all samples, counts + GDC + CDR

    The existing TUMOR_ONLY / NORMAL_ONLY / ALL files (counts + GDC clinical
    only) are produced separately by merge_outputs() and are not duplicated here.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Gene x sample count matrix.
    metadata_with_cdr : pd.DataFrame
        Sample metadata already joined with CDR annotations (from cdr.py).
    project_id : str
    output_dir : Path
    verbose : bool

    Returns
    -------
    dict
        Keys: 'full_with_cdr'. Values: Path objects.
    """
    import pandas as pd
    safe_pid = project_id.replace("/", "_")

    # Transpose counts so samples are rows
    counts_T = (
        count_matrix.T
        .reset_index()
        .rename(columns={"index": "sample_submitter_id"})
    )

    merged = counts_T.merge(metadata_with_cdr, on="sample_submitter_id", how="left")

    # Column ordering: audit flags first, then priority metadata, then cdr_ cols, then genes
    gene_cols  = list(count_matrix.index)
    audit_cols = ["cdr_matched", "cdr_subtype_available", "cdr_survival_complete"]
    cdr_cols   = [c for c in merged.columns
                  if c.startswith("cdr_") and c not in audit_cols]
    meta_cols  = [c for c in PRIORITY_COLUMNS if c in merged.columns]
    other_meta = [c for c in merged.columns
                  if c not in audit_cols + cdr_cols + meta_cols + gene_cols]

    ordered_cols = [c for c in audit_cols if c in merged.columns]
    ordered_cols += meta_cols
    ordered_cols += [c for c in other_meta if c not in ordered_cols]
    ordered_cols += cdr_cols
    ordered_cols += [g for g in gene_cols if g in merged.columns]

    final = merged[[c for c in ordered_cols if c in merged.columns]]

    path = output_dir / f"{safe_pid}_FULL_merged_with_CDR.tsv"
    final.to_csv(path, sep="\t", index=False)
    size_mb = path.stat().st_size / (1024 ** 2)

    if verbose:
        print(f"  💾  {path.name}  ({len(final)} samples, {size_mb:.1f} MB)")

    return {"full_with_cdr": path}
