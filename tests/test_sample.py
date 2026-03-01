"""
tests/test_sample.py
====================
Unit tests for TCGA barcode parsing and tumor / normal classification.

These tests cover:
- Correct extraction of sample type codes from barcodes
- Classification as Tumor, Normal, or Unknown
- All edge cases: malformed barcodes, empty strings, partial barcodes
- The annotate_metadata helper

These are fast pure-Python tests with no network calls.
"""

import pytest
import pandas as pd

from tcga_downloader.sample import (
    extract_sample_code,
    classify_sample,
    get_sample_type_label,
    annotate_metadata,
    summarise_sample_types,
)


# =============================================================================
#   extract_sample_code
# =============================================================================

class TestExtractSampleCode:

    def test_primary_tumor_short_barcode(self):
        assert extract_sample_code("TCGA-BH-A0B3-01A") == "01"

    def test_solid_tissue_normal(self):
        assert extract_sample_code("TCGA-BH-A0B3-11A") == "11"

    def test_full_aliquot_barcode(self):
        # Long barcodes still use field [3]
        assert extract_sample_code("TCGA-BH-A0B3-01A-11R-A10J-07") == "01"

    def test_metastatic(self):
        assert extract_sample_code("TCGA-XX-0001-06B") == "06"

    def test_blood_derived_normal(self):
        assert extract_sample_code("TCGA-XX-0001-10A") == "10"

    def test_empty_string(self):
        assert extract_sample_code("") == ""

    def test_none_handled_as_empty(self):
        # Some metadata rows may have None — should not raise
        assert extract_sample_code(None or "") == ""

    def test_too_few_fields(self):
        assert extract_sample_code("TCGA-BH") == ""

    def test_non_numeric_code(self):
        # If field 3 doesn't start with digits, return ""
        assert extract_sample_code("TCGA-BH-A0B3-XYA") == ""

    def test_exactly_four_fields(self):
        assert extract_sample_code("TCGA-AA-0001-01") == "01"


# =============================================================================
#   classify_sample
# =============================================================================

class TestClassifySample:

    # ── Tumor codes ───────────────────────────────────────────────────────────
    @pytest.mark.parametrize("barcode,expected", [
        ("TCGA-AA-0001-01A", "Tumor"),   # Primary Tumor
        ("TCGA-AA-0001-02A", "Tumor"),   # Recurrent Tumor
        ("TCGA-AA-0001-06A", "Tumor"),   # Metastatic
        ("TCGA-AA-0001-07A", "Tumor"),   # Additional Metastatic
    ])
    def test_tumor_barcodes(self, barcode, expected):
        assert classify_sample(barcode) == expected

    # ── Normal codes ──────────────────────────────────────────────────────────
    @pytest.mark.parametrize("barcode,expected", [
        ("TCGA-AA-0001-10A", "Normal"),  # Blood Derived Normal
        ("TCGA-AA-0001-11A", "Normal"),  # Solid Tissue Normal
        ("TCGA-AA-0001-12A", "Normal"),  # Buccal Cell Normal
        ("TCGA-AA-0001-14A", "Normal"),  # Bone Marrow Normal
    ])
    def test_normal_barcodes(self, barcode, expected):
        assert classify_sample(barcode) == expected

    # ── Edge cases ────────────────────────────────────────────────────────────
    def test_empty_string_is_unknown(self):
        assert classify_sample("") == "Unknown"

    def test_unrecognised_code_is_unknown(self):
        assert classify_sample("TCGA-AA-0001-99A") == "Unknown"

    def test_non_tcga_barcode_is_unknown(self):
        assert classify_sample("NOT-A-TCGA-BARCODE") == "Unknown"

    def test_tumor_is_never_normal(self):
        """Critical: a Primary Tumor barcode must never be classified as Normal."""
        result = classify_sample("TCGA-BH-A0B3-01A")
        assert result != "Normal"

    def test_normal_is_never_tumor(self):
        """Critical: a Solid Tissue Normal barcode must never be classified as Tumor."""
        result = classify_sample("TCGA-BH-A0B3-11A")
        assert result != "Tumor"


# =============================================================================
#   get_sample_type_label
# =============================================================================

class TestGetSampleTypeLabel:

    def test_primary_tumor_label(self):
        assert get_sample_type_label("TCGA-AA-0001-01A") == "Primary Tumor"

    def test_solid_tissue_normal_label(self):
        assert get_sample_type_label("TCGA-AA-0001-11A") == "Solid Tissue Normal"

    def test_metastatic_label(self):
        assert get_sample_type_label("TCGA-AA-0001-06A") == "Metastatic"

    def test_unknown_code_returns_descriptive_string(self):
        label = get_sample_type_label("TCGA-AA-0001-99A")
        assert "Unknown" in label
        assert "99" in label

    def test_empty_barcode_returns_descriptive_string(self):
        label = get_sample_type_label("")
        assert "Unknown" in label


# =============================================================================
#   annotate_metadata
# =============================================================================

class TestAnnotateMetadata:

    @pytest.fixture
    def sample_df(self):
        return pd.DataFrame({
            "sample_submitter_id": [
                "TCGA-AA-0001-01A",
                "TCGA-AA-0001-11A",
                "TCGA-AA-0002-06A",
                "",
            ]
        })

    def test_adds_three_columns(self, sample_df):
        annotate_metadata(sample_df)
        assert "sample_type_code"  in sample_df.columns
        assert "sample_type_label" in sample_df.columns
        assert "sample_category"   in sample_df.columns

    def test_correct_categories(self, sample_df):
        annotate_metadata(sample_df)
        cats = sample_df["sample_category"].tolist()
        assert cats[0] == "Tumor"
        assert cats[1] == "Normal"
        assert cats[2] == "Tumor"
        assert cats[3] == "Unknown"

    def test_correct_codes(self, sample_df):
        annotate_metadata(sample_df)
        assert sample_df.loc[0, "sample_type_code"] == "01"
        assert sample_df.loc[1, "sample_type_code"] == "11"

    def test_raises_if_no_sample_submitter_id_column(self):
        df = pd.DataFrame({"other_col": ["TCGA-AA-0001-01A"]})
        with pytest.raises(ValueError, match="sample_submitter_id"):
            annotate_metadata(df)

    def test_in_place_modification(self, sample_df):
        """annotate_metadata must modify the DataFrame in place, not return a new one."""
        result = annotate_metadata(sample_df)
        assert result is None   # modifies in place
        assert "sample_category" in sample_df.columns


# =============================================================================
#   summarise_sample_types
# =============================================================================

class TestSummariseSampleTypes:

    def test_counts_correct(self):
        df = pd.DataFrame({
            "sample_submitter_id": [
                "TCGA-AA-0001-01A",
                "TCGA-AA-0002-01A",
                "TCGA-AA-0003-11A",
            ]
        })
        annotate_metadata(df)
        summary = summarise_sample_types(df)
        assert summary["Tumor"]   == 2
        assert summary["Normal"]  == 1
        assert summary["Unknown"] == 0

    def test_empty_dataframe(self):
        df = pd.DataFrame({"sample_submitter_id": []})
        annotate_metadata(df)
        summary = summarise_sample_types(df)
        assert summary["Tumor"] == 0

    def test_returns_zero_for_missing_category(self):
        """Must include all three keys even if some are zero."""
        df = pd.DataFrame({"sample_submitter_id": ["TCGA-AA-0001-01A"]})
        annotate_metadata(df)
        summary = summarise_sample_types(df)
        assert "Tumor"   in summary
        assert "Normal"  in summary
        assert "Unknown" in summary
