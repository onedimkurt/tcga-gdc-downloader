"""
tcga_downloader
===============
User-friendly TCGA RNA-seq count downloader from the NCI GDC portal.

Public API
----------
GDCClient               — authenticate and query the GDC REST API
build_count_matrix      — parse downloaded STAR-Counts TSVs into a DataFrame
merge_outputs           — join counts + GDC clinical data, write output files
save_full_merged_with_cdr — join counts + GDC clinical + CDR annotations
run_cdr_pipeline        — full PanCanAtlas CDR integration pipeline
classify_sample         — classify a TCGA barcode as Tumor or Normal
get_sample_type_label   — human-readable label for a TCGA sample type code
"""

from tcga_downloader.client import GDCClient
from tcga_downloader.matrix import build_count_matrix
from tcga_downloader.merge  import merge_outputs, save_full_merged_with_cdr
from tcga_downloader.sample import classify_sample, get_sample_type_label
from tcga_downloader.cdr    import run_cdr_pipeline

__version__ = "2.1.1"

__all__ = [
    "GDCClient",
    "build_count_matrix",
    "merge_outputs",
    "save_full_merged_with_cdr",
    "run_cdr_pipeline",
    "classify_sample",
    "get_sample_type_label",
]


def _gui_path():
    from pathlib import Path
    return Path(__file__).parent / "app.py"
