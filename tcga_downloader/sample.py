"""
sample.py
=========
TCGA barcode parsing and tumor / normal classification.

TCGA barcode anatomy
--------------------
TCGA - TSS - Participant - SampleType - Vial - Portion - Analyte - Plate - Center
 0       1        2             3          4       5         6         7       8

Field 3 (index 3, zero-based) holds a 3-character code whose first 2 digits
identify the sample type (tumor vs. normal) and whose third character is the
vial letter.

Examples:
  TCGA-BH-A0B3-01A   →  code "01" = Primary Tumor
  TCGA-BH-A0B3-11A   →  code "11" = Solid Tissue Normal
  TCGA-06-0645-01A-11R-A10J-07  (full aliquot barcode, same logic applies)
"""

from __future__ import annotations
from tcga_downloader.constants import SAMPLE_TYPE_CODES, TUMOR_CODES, NORMAL_CODES


def extract_sample_code(barcode: str) -> str:
    """
    Extract the two-digit sample type code from a TCGA barcode string.

    Parameters
    ----------
    barcode : str
        Any TCGA barcode or submitter ID containing at least 4 hyphen-
        separated fields (e.g. 'TCGA-BH-A0B3-01A').

    Returns
    -------
    str
        Two-digit code, e.g. '01', '11'. Returns '' if unparseable.

    Examples
    --------
    >>> extract_sample_code('TCGA-BH-A0B3-01A')
    '01'
    >>> extract_sample_code('TCGA-BH-A0B3-11A')
    '11'
    >>> extract_sample_code('not-a-tcga-barcode-xx')
    ''
    """
    if not barcode:
        return ""
    parts = barcode.split("-")
    if len(parts) < 4:
        return ""
    field = parts[3]          # e.g. '01A'
    if len(field) < 2:
        return ""
    code = field[:2]
    if not code.isdigit():
        return ""
    return code


def classify_sample(barcode: str) -> str:
    """
    Return 'Tumor', 'Normal', or 'Unknown' for a TCGA sample barcode.

    Parameters
    ----------
    barcode : str
        TCGA sample or aliquot barcode.

    Returns
    -------
    str
        One of 'Tumor', 'Normal', 'Unknown'.

    Examples
    --------
    >>> classify_sample('TCGA-BH-A0B3-01A')
    'Tumor'
    >>> classify_sample('TCGA-BH-A0B3-11A')
    'Normal'
    >>> classify_sample('UNKNOWN-XX')
    'Unknown'
    """
    code = extract_sample_code(barcode)
    if code in TUMOR_CODES:
        return "Tumor"
    if code in NORMAL_CODES:
        return "Normal"
    return "Unknown"


def get_sample_type_label(barcode: str) -> str:
    """
    Return the human-readable GDC sample type label for a TCGA barcode.

    Parameters
    ----------
    barcode : str

    Returns
    -------
    str
        E.g. 'Primary Tumor', 'Solid Tissue Normal', or 'Unknown (code: XX)'.

    Examples
    --------
    >>> get_sample_type_label('TCGA-BH-A0B3-01A')
    'Primary Tumor'
    >>> get_sample_type_label('TCGA-BH-A0B3-11A')
    'Solid Tissue Normal'
    """
    code = extract_sample_code(barcode)
    if not code:
        return "Unknown (unrecognised barcode)"
    return SAMPLE_TYPE_CODES.get(code, f"Unknown (code: {code!r})")


def annotate_metadata(metadata_df) -> None:
    """
    Add tumor/normal annotation columns to a metadata DataFrame in-place.

    Adds three columns:
      - sample_type_code   : two-digit TCGA code, e.g. '01'
      - sample_type_label  : human-readable name, e.g. 'Primary Tumor'
      - sample_category    : 'Tumor', 'Normal', or 'Unknown'

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Must contain a 'sample_submitter_id' column.
    """
    if "sample_submitter_id" not in metadata_df.columns:
        raise ValueError("metadata_df must contain a 'sample_submitter_id' column")

    metadata_df["sample_type_code"]  = metadata_df["sample_submitter_id"].map(extract_sample_code)
    metadata_df["sample_type_label"] = metadata_df["sample_submitter_id"].map(get_sample_type_label)
    metadata_df["sample_category"]   = metadata_df["sample_submitter_id"].map(classify_sample)


def summarise_sample_types(metadata_df) -> dict:
    """
    Return a summary dict of sample category counts.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Must have a 'sample_category' column (added by annotate_metadata).

    Returns
    -------
    dict
        E.g. {'Tumor': 1095, 'Normal': 113, 'Unknown': 0}
    """
    if "sample_category" not in metadata_df.columns:
        return {}
    counts = metadata_df["sample_category"].value_counts().to_dict()
    return {
        "Tumor":   counts.get("Tumor",   0),
        "Normal":  counts.get("Normal",  0),
        "Unknown": counts.get("Unknown", 0),
    }
