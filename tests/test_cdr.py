"""
tests/test_cdr.py
=================
Unit tests for the CDR integration module.

All tests are pure-unit — no network calls, no real Excel files.
The CDR Excel structure is faked using in-memory DataFrames written
to temporary files, so tests run offline and fast.

Tests cover:
  - TCGA program gate (accepts TCGA, rejects others)
  - CDR parsing: column selection, NA handling, cdr_ prefix
  - CDR join: left join preserves all samples, correct key used
  - Audit flags: cdr_matched, cdr_subtype_available, cdr_survival_complete
  - Coverage report: per-field counts, unmatched cases file
  - Complete-cases split: correct rows selected, files written
  - Post-2018 samples: unmatched rows handled correctly, not silently dropped
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from tcga_downloader.cdr import (
    is_tcga_project,
    is_cdr_supported,
    parse_cdr,
    join_cdr,
    audit_cdr_join,
    generate_coverage_report,
    split_complete_cases,
    ProgramNotSupportedError,
)
from tcga_downloader.constants import CDR_JOIN_KEY


# =============================================================================
#   Helpers — fake CDR data
# =============================================================================

def _make_fake_cdr_excel(path: Path, extra_rows: list[dict] | None = None) -> Path:
    """
    Write a minimal fake CDR Excel file to path and return path.

    Contains 3 real TCGA cases (matching metadata fixture) plus
    optional extra rows for edge case testing.
    """
    rows = [
        {
            "bcr_patient_barcode":                  "TCGA-AA-0001",
            "type":                                  "BRCA",
            "age_at_initial_pathologic_diagnosis":   52,
            "gender":                                "FEMALE",
            "race":                                  "WHITE",
            "vital_status":                          "Alive",
            "tumor_status":                          "TUMOR FREE",
            "OS":    0, "OS.time":  730,
            "DSS":   0, "DSS.time": 730,
            "DFI":   0, "DFI.time": 680,
            "PFI":   0, "PFI.time": 730,
            "Subtype_mRNA":  "LumA",
            "Subtype_miRNA": "miRNA cluster 1",
            "Subtype_protein": "#N/A",
            "Subtype_CNV_low": "[Not Available]",
            "Subtype_CNV_high": "",
            "ajcc_pathologic_tumor_stage": "Stage IIA",
            "histological_type": "Infiltrating Ductal Carcinoma",
            "histological_grade": "G2",
            "initial_pathologic_dx_year": 2009,
            "death_days_to": "[Not Applicable]",
            "last_contact_days_to": 730,
            "cause_of_death": "[Not Available]",
        },
        {
            "bcr_patient_barcode":                  "TCGA-AA-0002",
            "type":                                  "BRCA",
            "age_at_initial_pathologic_diagnosis":   61,
            "gender":                                "FEMALE",
            "race":                                  "WHITE",
            "vital_status":                          "Dead",
            "tumor_status":                          "WITH TUMOR",
            "OS":    1, "OS.time":  420,
            "DSS":   1, "DSS.time": 420,
            "DFI":   "#N/A", "DFI.time": "#N/A",
            "PFI":   1, "PFI.time": 390,
            "Subtype_mRNA":  "Basal",
            "Subtype_miRNA": "[Not Available]",
            "Subtype_protein": "[Not Available]",
            "Subtype_CNV_low": "[Not Available]",
            "Subtype_CNV_high": "[Not Available]",
            "ajcc_pathologic_tumor_stage": "Stage III",
            "histological_type": "Infiltrating Ductal Carcinoma",
            "histological_grade": "G3",
            "initial_pathologic_dx_year": 2011,
            "death_days_to": 420,
            "last_contact_days_to": "[Not Applicable]",
            "cause_of_death": "Cancer Related",
        },
        # Third case exists in CDR but OS.time is missing → survival incomplete
        {
            "bcr_patient_barcode":                  "TCGA-AA-0003",
            "type":                                  "BRCA",
            "age_at_initial_pathologic_diagnosis":   45,
            "gender":                                "FEMALE",
            "race":                                  "BLACK OR AFRICAN AMERICAN",
            "vital_status":                          "Alive",
            "tumor_status":                          "TUMOR FREE",
            "OS":    0, "OS.time":  "[Not Available]",
            "DSS":   0, "DSS.time": "[Not Available]",
            "DFI":   0, "DFI.time": "[Not Available]",
            "PFI":   0, "PFI.time": "[Not Available]",
            "Subtype_mRNA":  "[Not Available]",
            "Subtype_miRNA": "[Not Available]",
            "Subtype_protein": "[Not Available]",
            "Subtype_CNV_low": "[Not Available]",
            "Subtype_CNV_high": "[Not Available]",
            "ajcc_pathologic_tumor_stage": "[Not Available]",
            "histological_type": "Infiltrating Ductal Carcinoma",
            "histological_grade": "[Not Available]",
            "initial_pathologic_dx_year": "[Not Available]",
            "death_days_to": "[Not Applicable]",
            "last_contact_days_to": "[Not Available]",
            "cause_of_death": "[Not Available]",
        },
    ]
    if extra_rows:
        rows.extend(extra_rows)

    df = pd.DataFrame(rows)
    df.to_excel(path, index=False, sheet_name="TCGA-CDR")
    return path


@pytest.fixture
def fake_cdr_excel(tmp_path):
    return _make_fake_cdr_excel(tmp_path / "fake_cdr.xlsx")


@pytest.fixture
def fake_metadata_df():
    """4 samples: 2 from cases in CDR, 1 with survival gaps, 1 post-2018 (not in CDR)."""
    return pd.DataFrame({
        "sample_submitter_id": [
            "TCGA-AA-0001-01A",   # in CDR, complete
            "TCGA-AA-0002-01A",   # in CDR, complete
            "TCGA-AA-0003-01A",   # in CDR, survival incomplete
            "TCGA-AA-9999-01A",   # NOT in CDR — post-2018 or excluded
        ],
        "case_submitter_id": [
            "TCGA-AA-0001",
            "TCGA-AA-0002",
            "TCGA-AA-0003",
            "TCGA-AA-9999",       # no CDR entry
        ],
        "sample_category":   ["Tumor", "Tumor", "Tumor", "Tumor"],
        "file_id":           ["uuid-001", "uuid-002", "uuid-003", "uuid-004"],
    })


# =============================================================================
#   Program gate
# =============================================================================

class TestIsTcgaProject:

    def test_tcga_returns_true(self):
        assert is_tcga_project("TCGA", "TCGA-BRCA") is True

    def test_cptac_raises(self):
        with pytest.raises(ProgramNotSupportedError) as exc_info:
            is_tcga_project("CPTAC", "CPTAC-3")
        assert "CPTAC" in str(exc_info.value)
        assert "TCGA" in exc_info.value.fix

    def test_target_raises(self):
        with pytest.raises(ProgramNotSupportedError):
            is_tcga_project("TARGET", "TARGET-AML")

    def test_empty_program_raises(self):
        with pytest.raises(ProgramNotSupportedError):
            is_tcga_project("", "UNKNOWN-PROJ")

    def test_error_message_contains_project_id(self):
        with pytest.raises(ProgramNotSupportedError) as exc_info:
            is_tcga_project("CPTAC", "CPTAC-3")
        assert "CPTAC-3" in exc_info.value.message

    def test_error_fix_mentions_gdc_portal(self):
        with pytest.raises(ProgramNotSupportedError) as exc_info:
            is_tcga_project("CGCI", "CGCI-BLGSP")
        assert "portal.gdc.cancer.gov" in exc_info.value.fix

    def test_is_cdr_supported_tcga(self):
        assert is_cdr_supported("TCGA") is True

    def test_is_cdr_supported_other(self):
        assert is_cdr_supported("CPTAC") is False


# =============================================================================
#   CDR parsing
# =============================================================================

class TestParseCdr:

    def test_returns_dataframe(self, fake_cdr_excel):
        df = parse_cdr(fake_cdr_excel, verbose=False)
        assert isinstance(df, pd.DataFrame)

    def test_has_join_key_column(self, fake_cdr_excel):
        df = parse_cdr(fake_cdr_excel, verbose=False)
        assert CDR_JOIN_KEY in df.columns

    def test_columns_prefixed_with_cdr(self, fake_cdr_excel):
        df = parse_cdr(fake_cdr_excel, verbose=False)
        non_join_cols = [c for c in df.columns if c != CDR_JOIN_KEY]
        for col in non_join_cols:
            assert col.startswith("cdr_"), f"Column '{col}' missing 'cdr_' prefix"

    def test_hash_na_converted_to_nan(self, fake_cdr_excel):
        df = parse_cdr(fake_cdr_excel, verbose=False)
        # Subtype_protein was '#N/A' for first case
        if "cdr_Subtype_protein" in df.columns:
            first_row = df[df[CDR_JOIN_KEY] == "TCGA-AA-0001"].iloc[0]
            assert pd.isna(first_row["cdr_Subtype_protein"])

    def test_not_available_converted_to_nan(self, fake_cdr_excel):
        df = parse_cdr(fake_cdr_excel, verbose=False)
        if "cdr_Subtype_CNV_low" in df.columns:
            first_row = df[df[CDR_JOIN_KEY] == "TCGA-AA-0001"].iloc[0]
            assert pd.isna(first_row["cdr_Subtype_CNV_low"])

    def test_survival_columns_are_numeric(self, fake_cdr_excel):
        df = parse_cdr(fake_cdr_excel, verbose=False)
        for col in ("cdr_OS", "cdr_OS.time", "cdr_PFI", "cdr_PFI.time"):
            if col in df.columns:
                # Should be numeric (int/float), not string
                assert pd.api.types.is_numeric_dtype(df[col]) or df[col].isna().all(), \
                    f"{col} should be numeric"

    def test_correct_number_of_rows(self, fake_cdr_excel):
        df = parse_cdr(fake_cdr_excel, verbose=False)
        assert len(df) == 3

    def test_missing_join_key_raises(self, tmp_path):
        bad_df = pd.DataFrame({"wrong_column": ["TCGA-AA-0001"]})
        bad_path = tmp_path / "bad_cdr.xlsx"
        bad_df.to_excel(bad_path, index=False)
        from tcga_downloader.exceptions import GDCError
        with pytest.raises(GDCError, match="join key"):
            parse_cdr(bad_path, verbose=False)


# =============================================================================
#   CDR join
# =============================================================================

class TestJoinCdr:

    @pytest.fixture
    def parsed_cdr(self, fake_cdr_excel):
        return parse_cdr(fake_cdr_excel, verbose=False)

    def test_all_metadata_rows_preserved(self, fake_metadata_df, parsed_cdr):
        """Left join must keep ALL samples, including post-2018 unmatched ones."""
        merged = join_cdr(fake_metadata_df, parsed_cdr, verbose=False)
        assert len(merged) == len(fake_metadata_df)

    def test_post_2018_sample_not_dropped(self, fake_metadata_df, parsed_cdr):
        """TCGA-AA-9999 has no CDR entry — must be in output with NaN CDR columns."""
        merged = join_cdr(fake_metadata_df, parsed_cdr, verbose=False)
        assert "TCGA-AA-9999-01A" in merged["sample_submitter_id"].values

    def test_post_2018_sample_has_nan_cdr_columns(self, fake_metadata_df, parsed_cdr):
        merged = join_cdr(fake_metadata_df, parsed_cdr, verbose=False)
        row = merged[merged["case_submitter_id"] == "TCGA-AA-9999"].iloc[0]
        cdr_cols = [c for c in merged.columns if c.startswith("cdr_")]
        for col in cdr_cols:
            assert pd.isna(row[col]), f"{col} should be NaN for post-2018 sample"

    def test_matched_sample_has_cdr_data(self, fake_metadata_df, parsed_cdr):
        merged = join_cdr(fake_metadata_df, parsed_cdr, verbose=False)
        row = merged[merged["case_submitter_id"] == "TCGA-AA-0001"].iloc[0]
        if "cdr_Subtype_mRNA" in merged.columns:
            assert row["cdr_Subtype_mRNA"] == "LumA"

    def test_no_duplicate_rows_created(self, fake_metadata_df, parsed_cdr):
        merged = join_cdr(fake_metadata_df, parsed_cdr, verbose=False)
        assert len(merged) == len(fake_metadata_df)

    def test_bcr_patient_barcode_column_dropped(self, fake_metadata_df, parsed_cdr):
        """The CDR join key should not appear as a duplicate column in output."""
        merged = join_cdr(fake_metadata_df, parsed_cdr, verbose=False)
        assert CDR_JOIN_KEY not in merged.columns

    def test_raises_without_case_submitter_id(self, parsed_cdr):
        from tcga_downloader.exceptions import GDCError
        bad_df = pd.DataFrame({"sample_submitter_id": ["TCGA-AA-0001-01A"]})
        with pytest.raises(GDCError, match="case_submitter_id"):
            join_cdr(bad_df, parsed_cdr, verbose=False)


# =============================================================================
#   Audit flags
# =============================================================================

class TestAuditCdrJoin:

    @pytest.fixture
    def merged_df(self, fake_metadata_df, fake_cdr_excel):
        cdr_df = parse_cdr(fake_cdr_excel, verbose=False)
        joined = join_cdr(fake_metadata_df, cdr_df, verbose=False)
        return audit_cdr_join(joined)

    def test_cdr_matched_column_added(self, merged_df):
        assert "cdr_matched" in merged_df.columns

    def test_cdr_subtype_available_column_added(self, merged_df):
        assert "cdr_subtype_available" in merged_df.columns

    def test_cdr_survival_complete_column_added(self, merged_df):
        assert "cdr_survival_complete" in merged_df.columns

    def test_post_2018_sample_not_matched(self, merged_df):
        row = merged_df[merged_df["case_submitter_id"] == "TCGA-AA-9999"].iloc[0]
        assert row["cdr_matched"] is False or row["cdr_matched"] == False

    def test_complete_sample_is_matched(self, merged_df):
        row = merged_df[merged_df["case_submitter_id"] == "TCGA-AA-0001"].iloc[0]
        assert row["cdr_matched"] is True or row["cdr_matched"] == True

    def test_subtype_available_when_populated(self, merged_df):
        row = merged_df[merged_df["case_submitter_id"] == "TCGA-AA-0001"].iloc[0]
        # TCGA-AA-0001 has Subtype_mRNA = LumA
        assert row["cdr_subtype_available"] == True

    def test_subtype_not_available_when_missing(self, merged_df):
        # TCGA-AA-0003 has Subtype_mRNA = [Not Available] → NaN after parsing
        row = merged_df[merged_df["case_submitter_id"] == "TCGA-AA-0003"].iloc[0]
        assert row["cdr_subtype_available"] == False

    def test_survival_complete_when_all_present(self, merged_df):
        # TCGA-AA-0001 has OS, OS.time, PFI, PFI.time all present
        row = merged_df[merged_df["case_submitter_id"] == "TCGA-AA-0001"].iloc[0]
        assert row["cdr_survival_complete"] == True

    def test_survival_incomplete_when_os_time_missing(self, merged_df):
        # TCGA-AA-0003 has OS.time = NaN
        row = merged_df[merged_df["case_submitter_id"] == "TCGA-AA-0003"].iloc[0]
        assert row["cdr_survival_complete"] == False

    def test_post_2018_survival_incomplete(self, merged_df):
        # TCGA-AA-9999 has no CDR data at all
        row = merged_df[merged_df["case_submitter_id"] == "TCGA-AA-9999"].iloc[0]
        assert row["cdr_survival_complete"] == False


# =============================================================================
#   Coverage report
# =============================================================================

class TestGenerateCoverageReport:

    @pytest.fixture
    def audited_df(self, fake_metadata_df, fake_cdr_excel, tmp_path):
        cdr_df = parse_cdr(fake_cdr_excel, verbose=False)
        joined = join_cdr(fake_metadata_df, cdr_df, verbose=False)
        return audit_cdr_join(joined)

    def test_coverage_report_file_created(self, audited_df, tmp_path):
        generate_coverage_report(audited_df, "TCGA-TEST", tmp_path, verbose=False)
        assert (tmp_path / "TCGA-TEST_CDR_coverage_report.tsv").exists()

    def test_coverage_report_has_expected_columns(self, audited_df, tmp_path):
        report = generate_coverage_report(
            audited_df, "TCGA-TEST", tmp_path, verbose=False
        )
        assert "field" in report.columns
        assert "n_present" in report.columns
        assert "pct_present" in report.columns

    def test_unmatched_file_created_when_unmatched_exist(self, audited_df, tmp_path):
        generate_coverage_report(audited_df, "TCGA-TEST", tmp_path, verbose=False)
        # TCGA-AA-9999 is unmatched → file should be created
        assert (tmp_path / "TCGA-TEST_CDR_unmatched_cases.txt").exists()

    def test_unmatched_file_contains_correct_barcode(self, audited_df, tmp_path):
        generate_coverage_report(audited_df, "TCGA-TEST", tmp_path, verbose=False)
        content = (tmp_path / "TCGA-TEST_CDR_unmatched_cases.txt").read_text()
        assert "TCGA-AA-9999" in content

    def test_unmatched_file_not_created_when_all_matched(self, tmp_path, fake_cdr_excel):
        """If every case matches CDR, the unmatched file should not be created."""
        all_matched_meta = pd.DataFrame({
            "sample_submitter_id": ["TCGA-AA-0001-01A", "TCGA-AA-0002-01A"],
            "case_submitter_id":   ["TCGA-AA-0001",     "TCGA-AA-0002"],
            "sample_category":     ["Tumor",             "Tumor"],
            "file_id":             ["uuid-001",          "uuid-002"],
        })
        cdr_df = parse_cdr(fake_cdr_excel, verbose=False)
        joined = join_cdr(all_matched_meta, cdr_df, verbose=False)
        audited = audit_cdr_join(joined)
        generate_coverage_report(audited, "TCGA-TEST", tmp_path, verbose=False)
        assert not (tmp_path / "TCGA-TEST_CDR_unmatched_cases.txt").exists()

    def test_pct_present_between_0_and_100(self, audited_df, tmp_path):
        report = generate_coverage_report(
            audited_df, "TCGA-TEST", tmp_path, verbose=False
        )
        assert (report["pct_present"] >= 0).all()
        assert (report["pct_present"] <= 100).all()


# =============================================================================
#   Complete-cases split
# =============================================================================

class TestSplitCompleteCases:

    @pytest.fixture
    def audited_df(self, fake_metadata_df, fake_cdr_excel, tmp_path):
        cdr_df = parse_cdr(fake_cdr_excel, verbose=False)
        joined = join_cdr(fake_metadata_df, cdr_df, verbose=False)
        return audit_cdr_join(joined)

    def test_complete_plus_incomplete_equals_total(self, audited_df, tmp_path):
        complete, incomplete = split_complete_cases(
            audited_df, "TCGA-TEST", tmp_path, verbose=False
        )
        assert len(complete) + len(incomplete) == len(audited_df)

    def test_complete_cases_file_written(self, audited_df, tmp_path):
        split_complete_cases(audited_df, "TCGA-TEST", tmp_path, verbose=False)
        assert (tmp_path / "TCGA-TEST_FULL_merged_CDR_complete_cases.tsv").exists()

    def test_complete_cases_all_have_survival(self, audited_df, tmp_path):
        complete, _ = split_complete_cases(
            audited_df, "TCGA-TEST", tmp_path, verbose=False
        )
        if len(complete) > 0:
            assert complete["cdr_survival_complete"].all()

    def test_incomplete_cases_have_missing_survival(self, audited_df, tmp_path):
        _, incomplete = split_complete_cases(
            audited_df, "TCGA-TEST", tmp_path, verbose=False
        )
        # Post-2018 and survival-incomplete cases should all be in incomplete
        assert "TCGA-AA-9999-01A" in incomplete["sample_submitter_id"].values
        assert "TCGA-AA-0003-01A" in incomplete["sample_submitter_id"].values

    def test_post_2018_never_in_complete_cases(self, audited_df, tmp_path):
        """Critical: post-2018 unmatched samples must never appear in complete cases."""
        complete, _ = split_complete_cases(
            audited_df, "TCGA-TEST", tmp_path, verbose=False
        )
        assert "TCGA-AA-9999-01A" not in complete["sample_submitter_id"].values
