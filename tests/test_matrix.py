"""
tests/test_matrix.py
====================
Unit tests for STAR count file parsing and count matrix construction.

Tests cover:
- Current GDC column format (6 columns)
- Older GDC column format (4 columns)
- Compressed (.tsv.gz) files
- Files with no header
- Files with only N_ rows (edge case)
- Mixed parseable / unparseable files
- All-zero sample detection
- Correct gene count values
"""

import gzip
import textwrap
from pathlib import Path

import pandas as pd
import pytest

from tcga_downloader.matrix import (
    detect_column_map,
    parse_star_file,
    build_count_matrix,
    _report_matrix_qc,
)
from tcga_downloader.exceptions import ColumnDetectionError, GDCError


# =============================================================================
#   Helpers
# =============================================================================

def write_star_file(path: Path, content: str, compress: bool = False) -> Path:
    if compress:
        path = path.with_suffix(".tsv.gz")
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write(content)
    else:
        path.write_text(content, encoding="utf-8")
    return path


CURRENT_FORMAT = textwrap.dedent("""\
    gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second
    N_unmapped\t\t\t100\t50\t50
    N_multimapping\t\t\t200\t100\t100
    N_noFeature\t\t\t300\t150\t150
    N_ambiguous\t\t\t50\t25\t25
    ENSG00000000003\tTSPAN6\tprotein_coding\t1500\t750\t750
    ENSG00000000005\tTNMD\tprotein_coding\t0\t0\t0
    ENSG00000000419\tDPM1\tprotein_coding\t2300\t1150\t1150
""")

OLD_FORMAT = textwrap.dedent("""\
    gene_id\tunstranded\tstranded_first\tstranded_second
    N_unmapped\t100\t50\t50
    ENSG00000000003\t1500\t750\t750
    ENSG00000000419\t2300\t1150\t1150
""")

TWO_COLUMN_FORMAT = textwrap.dedent("""\
    gene_id\tcount
    N_unmapped\t100
    ENSG00000000003\t1500
    ENSG00000000419\t2300
""")


# =============================================================================
#   detect_column_map
# =============================================================================

class TestDetectColumnMap:

    def test_current_format_detects_unstranded(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT)
        col_map = detect_column_map(fp)
        assert "unstranded" in col_map
        assert col_map["unstranded"] == 3

    def test_current_format_detects_stranded_columns(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT)
        col_map = detect_column_map(fp)
        assert col_map.get("stranded_first")  == 4
        assert col_map.get("stranded_second") == 5

    def test_old_format_detects_unstranded(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", OLD_FORMAT)
        col_map = detect_column_map(fp)
        assert col_map["unstranded"] == 1

    def test_two_column_format(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", TWO_COLUMN_FORMAT)
        col_map = detect_column_map(fp)
        assert "unstranded" in col_map

    def test_gzipped_file(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT, compress=True)
        col_map = detect_column_map(fp)
        assert col_map["unstranded"] == 3

    def test_unknown_columns_raise_error(self, tmp_path):
        bad_content = "gene_id\tfoo\tbar\tbaz\n"
        fp = write_star_file(tmp_path / "test.tsv", bad_content)
        with pytest.raises(ColumnDetectionError):
            detect_column_map(fp)

    def test_empty_file_raises_error(self, tmp_path):
        fp = tmp_path / "empty.tsv"
        fp.write_text("", encoding="utf-8")
        with pytest.raises((ColumnDetectionError, GDCError)):
            detect_column_map(fp)


# =============================================================================
#   parse_star_file
# =============================================================================

class TestParseStarFile:

    def test_current_format_correct_values(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT)
        s = parse_star_file(fp, unstranded_col_idx=3)
        assert s is not None
        assert s["ENSG00000000003"] == 1500
        assert s["ENSG00000000419"] == 2300

    def test_zero_count_gene_included(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT)
        s = parse_star_file(fp, unstranded_col_idx=3)
        assert s["ENSG00000000005"] == 0

    def test_n_rows_excluded(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT)
        s = parse_star_file(fp, unstranded_col_idx=3)
        assert s is not None
        for key in s.index:
            assert not key.startswith("N_"), f"N_ row '{key}' should have been excluded"

    def test_header_row_excluded(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT)
        s = parse_star_file(fp, unstranded_col_idx=3)
        assert "gene_id" not in s.index

    def test_old_format_correct_values(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", OLD_FORMAT)
        s = parse_star_file(fp, unstranded_col_idx=1)
        assert s["ENSG00000000003"] == 1500

    def test_gzipped_file_parses_correctly(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT, compress=True)
        s = parse_star_file(fp, unstranded_col_idx=3)
        assert s is not None
        assert s["ENSG00000000003"] == 1500

    def test_returns_none_for_empty_file(self, tmp_path):
        fp = tmp_path / "empty.tsv"
        fp.write_text("", encoding="utf-8")
        s = parse_star_file(fp, unstranded_col_idx=3)
        assert s is None

    def test_returns_none_for_only_n_rows(self, tmp_path):
        content = "N_unmapped\t100\t50\t50\nN_noFeature\t200\t100\t100\n"
        fp = write_star_file(tmp_path / "test.tsv", content)
        s = parse_star_file(fp, unstranded_col_idx=1)
        assert s is None

    def test_dtype_is_int64(self, tmp_path):
        fp = write_star_file(tmp_path / "test.tsv", CURRENT_FORMAT)
        s = parse_star_file(fp, unstranded_col_idx=3)
        assert s.dtype == "int64"


# =============================================================================
#   build_count_matrix
# =============================================================================

class TestBuildCountMatrix:

    @pytest.fixture
    def two_sample_dir(self, tmp_path, fake_metadata_df):
        """Create a raw_counts directory with two properly-formatted STAR files."""
        counts_dir = tmp_path / "raw_counts"
        for i, row in fake_metadata_df.iterrows():
            d = counts_dir / f"uuid-{i:03d}"
            d.mkdir(parents=True, exist_ok=True)
            lines = ["gene_id\tunstranded\tstranded_first\tstranded_second"]
            lines.append("N_unmapped\t50\t25\t25")
            for j, g in enumerate(["ENSG00000000003", "ENSG00000000419", "ENSG00000000457"]):
                count = (j + 1) * 100 + i * 500
                lines.append(f"{g}\t{count}\t{count//2}\t{count//2}")
            (d / row["file_name"]).write_text("\n".join(lines), encoding="utf-8")
        return counts_dir

    def test_matrix_shape(self, two_sample_dir, fake_metadata_df):
        matrix = build_count_matrix(two_sample_dir, fake_metadata_df, verbose=False)
        assert matrix.shape == (3, 2)  # 3 genes × 2 samples

    def test_gene_ids_as_index(self, two_sample_dir, fake_metadata_df):
        matrix = build_count_matrix(two_sample_dir, fake_metadata_df, verbose=False)
        assert "ENSG00000000003" in matrix.index
        assert "N_unmapped" not in matrix.index

    def test_sample_ids_as_columns(self, two_sample_dir, fake_metadata_df):
        matrix = build_count_matrix(two_sample_dir, fake_metadata_df, verbose=False)
        # Columns should be sample_submitter_id values from metadata
        expected_samples = set(fake_metadata_df["sample_submitter_id"])
        assert set(matrix.columns) == expected_samples

    def test_count_values_are_integers(self, two_sample_dir, fake_metadata_df):
        matrix = build_count_matrix(two_sample_dir, fake_metadata_df, verbose=False)
        for col in matrix.columns:
            assert matrix[col].dtype in ("int64", "int32")

    def test_raises_if_no_tsv_files(self, tmp_path, fake_metadata_df):
        empty_dir = tmp_path / "empty_counts"
        empty_dir.mkdir()
        with pytest.raises(GDCError, match="No TSV files"):
            build_count_matrix(empty_dir, fake_metadata_df, verbose=False)

    def test_unmatched_files_get_fallback_id(self, tmp_path):
        """Files not in metadata should get UNMATCHED_ prefix, not raise an error."""
        counts_dir = tmp_path / "raw_counts"
        d = counts_dir / "unknown-uuid"
        d.mkdir(parents=True, exist_ok=True)
        content = "gene_id\tunstranded\n" + "N_unmapped\t50\n" + "ENSG00000000003\t100\n"
        (d / "unexpected_file.tsv").write_text(content, encoding="utf-8")

        meta = pd.DataFrame({
            "file_name":           ["different_name.tsv"],
            "sample_submitter_id": ["TCGA-AA-0001-01A"],
        })
        matrix = build_count_matrix(counts_dir, meta, verbose=False)
        # Should not raise; column name starts with UNMATCHED_
        assert any(c.startswith("UNMATCHED_") for c in matrix.columns)

    def test_duplicate_sample_ids_deduplicated(self, tmp_path):
        """Two files mapped to the same sample_submitter_id should get _dup suffix."""
        counts_dir = tmp_path / "raw_counts"
        content = "gene_id\tunstranded\n" + "ENSG00000000003\t100\n"
        for i in range(2):
            d = counts_dir / f"uuid-{i}"
            d.mkdir(parents=True, exist_ok=True)
            (d / f"sample_{i}.tsv").write_text(content, encoding="utf-8")

        meta = pd.DataFrame({
            "file_name":           [f"sample_{i}.tsv" for i in range(2)],
            "sample_submitter_id": ["TCGA-AA-0001-01A", "TCGA-AA-0001-01A"],  # same ID
        })
        matrix = build_count_matrix(counts_dir, meta, verbose=False)
        # Should have 2 columns, second with _dup suffix
        assert matrix.shape[1] == 2
        assert "TCGA-AA-0001-01A" in matrix.columns
        assert any("_dup" in c for c in matrix.columns)
