"""
tests/test_merge.py
===================
Unit tests for the merge layer: clinical flattening, metadata flattening,
and the final merge_outputs function.

Tests cover:
- Clinical JSON flattening (nested diagnoses, demographics, exposures)
- Metadata flattening with tumor/normal annotation
- Tumor and normal sample separation in output files
- Column ordering in merged output
- Handling of missing clinical data (unmatched cases)
- Output file creation and content
"""

import pandas as pd
import pytest
from pathlib import Path

from tcga_downloader.merge import (
    flatten_clinical,
    flatten_file_metadata,
    merge_outputs,
    PRIORITY_COLUMNS,
)
from tests.conftest import FAKE_FILES_RESPONSE, FAKE_CASES_RESPONSE


# =============================================================================
#   flatten_clinical
# =============================================================================

class TestFlattenClinical:

    def test_returns_dataframe(self):
        df = flatten_clinical(FAKE_CASES_RESPONSE["data"]["hits"])
        assert isinstance(df, pd.DataFrame)

    def test_correct_number_of_rows(self):
        df = flatten_clinical(FAKE_CASES_RESPONSE["data"]["hits"])
        assert len(df) == 2

    def test_case_submitter_id_present(self):
        df = flatten_clinical(FAKE_CASES_RESPONSE["data"]["hits"])
        assert "case_submitter_id" in df.columns
        assert "TCGA-AA-0001" in df["case_submitter_id"].values

    def test_demographic_fields_extracted(self):
        df = flatten_clinical(FAKE_CASES_RESPONSE["data"]["hits"])
        row = df[df["case_submitter_id"] == "TCGA-AA-0001"].iloc[0]
        assert row["gender"] == "female"
        assert int(row["age_at_index"]) == 52
        assert row["vital_status"] == "Alive"

    def test_diagnosis_fields_extracted(self):
        df = flatten_clinical(FAKE_CASES_RESPONSE["data"]["hits"])
        row = df[df["case_submitter_id"] == "TCGA-AA-0001"].iloc[0]
        assert row["tumor_stage"] == "stage ii"
        assert row["ajcc_pathologic_t"] == "T2"

    def test_exposure_fields_extracted(self):
        df = flatten_clinical(FAKE_CASES_RESPONSE["data"]["hits"])
        row = df[df["case_submitter_id"] == "TCGA-AA-0001"].iloc[0]
        assert row["alcohol_history"] == "No"

    def test_missing_exposures_handled(self):
        """Case with empty exposures list should not raise."""
        df = flatten_clinical(FAKE_CASES_RESPONSE["data"]["hits"])
        row = df[df["case_submitter_id"] == "TCGA-AA-0002"].iloc[0]
        # Should be empty string, not raise KeyError
        assert row["alcohol_history"] == "" or pd.isna(row["alcohol_history"])

    def test_empty_input_returns_empty_df(self):
        df = flatten_clinical([])
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0


# =============================================================================
#   flatten_file_metadata
# =============================================================================

class TestFlattenFileMetadata:

    def test_returns_dataframe(self):
        df = flatten_file_metadata(FAKE_FILES_RESPONSE["data"]["hits"])
        assert isinstance(df, pd.DataFrame)

    def test_correct_number_of_rows(self):
        df = flatten_file_metadata(FAKE_FILES_RESPONSE["data"]["hits"])
        assert len(df) == 2

    def test_sample_category_added(self):
        """annotate_metadata must have been called — sample_category must exist."""
        df = flatten_file_metadata(FAKE_FILES_RESPONSE["data"]["hits"])
        assert "sample_category" in df.columns

    def test_tumor_classified_correctly(self):
        df = flatten_file_metadata(FAKE_FILES_RESPONSE["data"]["hits"])
        tumor_row = df[df["sample_submitter_id"] == "TCGA-AA-0001-01A"].iloc[0]
        assert tumor_row["sample_category"] == "Tumor"
        assert tumor_row["sample_type_label"] == "Primary Tumor"

    def test_normal_classified_correctly(self):
        df = flatten_file_metadata(FAKE_FILES_RESPONSE["data"]["hits"])
        normal_row = df[df["sample_submitter_id"] == "TCGA-AA-0002-11A"].iloc[0]
        assert normal_row["sample_category"] == "Normal"
        assert normal_row["sample_type_label"] == "Solid Tissue Normal"

    def test_tumor_and_normal_never_swapped(self):
        """Regression test: tumour row must not be labelled Normal, and vice versa."""
        df = flatten_file_metadata(FAKE_FILES_RESPONSE["data"]["hits"])
        tumor_cats  = df[df["sample_type_code"] == "01"]["sample_category"].unique()
        normal_cats = df[df["sample_type_code"] == "11"]["sample_category"].unique()
        assert "Normal" not in tumor_cats
        assert "Tumor"  not in normal_cats


# =============================================================================
#   merge_outputs
# =============================================================================

class TestMergeOutputs:

    def test_all_file_created(self, tmp_path, fake_count_matrix, fake_metadata_df,
                               fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        assert "all" in paths
        assert paths["all"].exists()

    def test_tumor_file_created(self, tmp_path, fake_count_matrix, fake_metadata_df,
                                 fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        assert "tumor" in paths
        assert paths["tumor"].exists()

    def test_normal_file_created(self, tmp_path, fake_count_matrix, fake_metadata_df,
                                  fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        assert "normal" in paths
        assert paths["normal"].exists()

    def test_metadata_only_file_created(self, tmp_path, fake_count_matrix,
                                         fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        assert "metadata" in paths
        assert paths["metadata"].exists()

    def test_tumor_only_file_has_only_tumor_samples(self, tmp_path, fake_count_matrix,
                                                      fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        tumor_df = pd.read_csv(paths["tumor"], sep="\t")
        assert set(tumor_df["sample_category"].unique()) == {"Tumor"}
        assert "Normal" not in tumor_df["sample_category"].values

    def test_normal_only_file_has_only_normal_samples(self, tmp_path, fake_count_matrix,
                                                        fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        normal_df = pd.read_csv(paths["normal"], sep="\t")
        assert set(normal_df["sample_category"].unique()) == {"Normal"}
        assert "Tumor" not in normal_df["sample_category"].values

    def test_sample_category_column_in_all_file(self, tmp_path, fake_count_matrix,
                                                  fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        df = pd.read_csv(paths["all"], sep="\t")
        assert "sample_category" in df.columns

    def test_gene_columns_present_in_all_file(self, tmp_path, fake_count_matrix,
                                               fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        df = pd.read_csv(paths["all"], sep="\t")
        gene_ids = list(fake_count_matrix.index)
        for gene in gene_ids:
            assert gene in df.columns, f"Gene {gene} missing from output"

    def test_clinical_columns_present(self, tmp_path, fake_count_matrix,
                                       fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        df = pd.read_csv(paths["all"], sep="\t")
        assert "gender"          in df.columns
        assert "vital_status"    in df.columns
        assert "tumor_stage"     in df.columns

    def test_priority_columns_come_before_genes(self, tmp_path, fake_count_matrix,
                                                  fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        df = pd.read_csv(paths["all"], sep="\t")
        all_cols = list(df.columns)
        gene_positions = [all_cols.index(g) for g in fake_count_matrix.index if g in all_cols]
        meta_positions = [all_cols.index(c) for c in PRIORITY_COLUMNS if c in all_cols]
        if gene_positions and meta_positions:
            assert max(meta_positions) < min(gene_positions), (
                "Metadata columns must appear before gene count columns"
            )

    def test_no_gene_columns_in_metadata_file(self, tmp_path, fake_count_matrix,
                                               fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        meta_df = pd.read_csv(paths["metadata"], sep="\t")
        for gene in fake_count_matrix.index:
            assert gene not in meta_df.columns

    def test_works_without_clinical_data(self, tmp_path, fake_count_matrix, fake_metadata_df):
        """merge_outputs must not fail if clinical_df is empty."""
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, pd.DataFrame(),
            "TCGA-TEST", tmp_path, verbose=False,
        )
        assert paths["all"].exists()

    def test_row_count_matches_sample_count(self, tmp_path, fake_count_matrix,
                                             fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        df = pd.read_csv(paths["all"], sep="\t")
        assert len(df) == fake_count_matrix.shape[1]

    def test_output_is_tab_separated(self, tmp_path, fake_count_matrix,
                                      fake_metadata_df, fake_clinical_df):
        paths = merge_outputs(
            fake_count_matrix, fake_metadata_df, fake_clinical_df,
            "TCGA-TEST", tmp_path, verbose=False,
        )
        first_line = paths["all"].read_text(encoding="utf-8").splitlines()[0]
        assert "\t" in first_line
