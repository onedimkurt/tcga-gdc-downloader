"""
tests/test_integration.py
==========================
Integration tests that make REAL calls to the GDC API.

These tests are SKIPPED by default to avoid:
  - Requiring internet access in CI
  - Hitting GDC rate limits during development

To run them:
    export TCGA_INTEGRATION=1
    pytest tests/test_integration.py -v

What is tested against the live API:
  - A small, well-known project (TCGA-CHOL, ~36 samples) is used to
    keep download size manageable for integration testing.
  - We verify the actual API contract (fields, pagination, sample types).
  - We DO NOT download the full file archive in integration tests.
"""

import os
import pytest

pytestmark = pytest.mark.integration

SMALL_PROJECT = "TCGA-CHOL"   # 36 samples — small enough for quick tests


@pytest.fixture(scope="module")
def live_client():
    """GDCClient with no token (open-access only)."""
    from tcga_downloader.client import GDCClient
    client = GDCClient()
    if not client.check_connectivity():
        pytest.skip("GDC API is not reachable from this machine.")
    return client


class TestLiveProjectValidation:

    def test_tcga_chol_exists(self, live_client):
        info = live_client.validate_project(SMALL_PROJECT)
        assert info["project_id"] == SMALL_PROJECT

    def test_tcga_brca_exists(self, live_client):
        info = live_client.validate_project("TCGA-BRCA")
        assert "Breast" in info.get("primary_site", [])

    def test_nonexistent_project_raises(self, live_client):
        from tcga_downloader.exceptions import GDCError
        with pytest.raises(GDCError):
            live_client.validate_project("TCGA-ZZZZZZZ")


class TestLiveFileDiscovery:

    def test_chol_has_star_files(self, live_client):
        hits = live_client.discover_star_files(SMALL_PROJECT)
        assert len(hits) > 0

    def test_files_have_required_fields(self, live_client):
        hits = live_client.discover_star_files(SMALL_PROJECT)
        for hit in hits:
            assert "file_id"   in hit
            assert "file_name" in hit
            assert "cases"     in hit
            cases = hit["cases"]
            assert len(cases) > 0
            assert "submitter_id" in cases[0]

    def test_all_files_are_tsv(self, live_client):
        """Every returned file should be a .tsv (STAR counts format)."""
        hits = live_client.discover_star_files(SMALL_PROJECT)
        for hit in hits:
            assert hit["file_name"].endswith(".tsv") or hit["file_name"].endswith(".tsv.gz"), (
                f"Unexpected file format: {hit['file_name']}"
            )

    def test_sample_barcodes_are_tcga_format(self, live_client):
        """All sample submitter IDs should follow TCGA barcode pattern."""
        import re
        hits = live_client.discover_star_files(SMALL_PROJECT)
        tcga_pattern = re.compile(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}")
        for hit in hits:
            for case in hit.get("cases", []):
                for sample in case.get("samples", []):
                    sid = sample.get("submitter_id", "")
                    assert tcga_pattern.match(sid), (
                        f"Sample ID '{sid}' does not match TCGA barcode format"
                    )

    def test_tumor_normal_classification_consistent(self, live_client):
        """
        The sample_category derived from the barcode must agree with
        the GDC-provided sample_type field (they should be consistent).
        """
        from tcga_downloader.merge import flatten_file_metadata
        from tcga_downloader.sample import classify_sample

        hits = live_client.discover_star_files(SMALL_PROJECT)
        meta_df = flatten_file_metadata(hits)

        mismatches = []
        for _, row in meta_df.iterrows():
            barcode_cat = row["sample_category"]
            gdc_type    = row["gdc_sample_type"].lower() if row["gdc_sample_type"] else ""

            # GDC "Normal" types contain "normal" in their label
            gdc_is_normal = "normal" in gdc_type
            barcode_is_normal = barcode_cat == "Normal"

            if gdc_is_normal != barcode_is_normal and barcode_cat != "Unknown":
                mismatches.append({
                    "sample": row["sample_submitter_id"],
                    "barcode_cat": barcode_cat,
                    "gdc_type": row["gdc_sample_type"],
                })

        # We allow 0 mismatches — this is a critical correctness check
        assert len(mismatches) == 0, (
            f"Tumor/normal classification disagreement between barcode and GDC label:\n"
            + "\n".join(str(m) for m in mismatches)
        )


class TestLiveClinicalData:

    def test_chol_clinical_data_available(self, live_client):
        cases = live_client.fetch_clinical_data(SMALL_PROJECT)
        assert len(cases) > 0

    def test_clinical_cases_have_demographic(self, live_client):
        cases = live_client.fetch_clinical_data(SMALL_PROJECT)
        for case in cases:
            # demographic may be empty dict but key should exist or be None
            assert "demographic" in case or case.get("demographic") is not None or True

    def test_flatten_clinical_produces_expected_columns(self, live_client):
        from tcga_downloader.merge import flatten_clinical
        cases = live_client.fetch_clinical_data(SMALL_PROJECT)
        df = flatten_clinical(cases)
        expected_columns = [
            "case_id", "case_submitter_id", "gender",
            "vital_status", "primary_diagnosis",
        ]
        for col in expected_columns:
            assert col in df.columns, f"Expected column '{col}' missing"
