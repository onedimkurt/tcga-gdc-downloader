"""
conftest.py
===========
Shared pytest fixtures used across all test modules.

Real GDC API calls are never made in the unit test suite.
Instead, the `responses` library intercepts HTTP requests
and returns controlled fake responses.

For integration tests (tests/test_integration.py), set the
environment variable:
    TCGA_INTEGRATION=1
to allow real network calls. These tests are skipped by default.
"""

from __future__ import annotations

import json
import os
import textwrap
from pathlib import Path

import pandas as pd
import pytest
import responses as responses_lib


# =============================================================================
#   Minimal fake GDC API responses
# =============================================================================

FAKE_PROJECT_ID = "TCGA-TEST"

FAKE_PROJECT_RESPONSE = {
    "data": {
        "hits": [{
            "project_id":    FAKE_PROJECT_ID,
            "name":          "Test Project",
            "primary_site":  ["Breast"],
            "disease_type":  ["Adenomas and Adenocarcinomas"],
        }],
        "pagination": {"total": 1, "count": 1, "page": 1, "pages": 1, "size": 1, "from": 0}
    }
}

FAKE_FILES_RESPONSE = {
    "data": {
        "hits": [
            {
                "file_id":   "550e8400-e29b-41d4-a716-446655440001",
                "file_name": "sample_001.star_counts.tsv",
                "file_size": 500_000,
                "cases": [{
                    "case_id":     "case-uuid-001",
                    "submitter_id": "TCGA-AA-0001",
                    "samples": [{
                        "sample_id":     "sample-uuid-001",
                        "submitter_id":  "TCGA-AA-0001-01A",  # Primary Tumor
                        "sample_type":   "Primary Tumor",
                        "tissue_type":   "Tumor",
                    }]
                }]
            },
            {
                "file_id":   "550e8400-e29b-41d4-a716-446655440002",
                "file_name": "sample_002.star_counts.tsv",
                "file_size": 490_000,
                "cases": [{
                    "case_id":     "case-uuid-002",
                    "submitter_id": "TCGA-AA-0002",
                    "samples": [{
                        "sample_id":     "sample-uuid-002",
                        "submitter_id":  "TCGA-AA-0002-11A",  # Solid Tissue Normal
                        "sample_type":   "Solid Tissue Normal",
                        "tissue_type":   "Normal",
                    }]
                }]
            },
        ],
        "pagination": {"total": 2, "count": 2, "page": 1, "pages": 1, "size": 2, "from": 0}
    }
}

FAKE_CASES_RESPONSE = {
    "data": {
        "hits": [
            {
                "case_id":     "case-uuid-001",
                "submitter_id": "TCGA-AA-0001",
                "demographic": {
                    "gender": "female", "age_at_index": 52, "race": "white",
                    "ethnicity": "not hispanic or latino", "vital_status": "Alive",
                    "days_to_death": None, "year_of_birth": 1960,
                },
                "diagnoses": [{
                    "age_at_diagnosis": 18980, "days_to_last_follow_up": 730,
                    "primary_diagnosis": "Infiltrating duct carcinoma, NOS",
                    "tumor_stage": "stage ii", "tumor_grade": "G2",
                    "ajcc_pathologic_stage": "Stage IIA",
                    "ajcc_pathologic_t": "T2", "ajcc_pathologic_n": "N0",
                    "ajcc_pathologic_m": "M0",
                }],
                "exposures": [{"alcohol_history": "No", "bmi": 24.5}],
            },
            {
                "case_id":     "case-uuid-002",
                "submitter_id": "TCGA-AA-0002",
                "demographic": {
                    "gender": "female", "age_at_index": 61, "race": "white",
                    "ethnicity": "not hispanic or latino", "vital_status": "Dead",
                    "days_to_death": 420, "year_of_birth": 1951,
                },
                "diagnoses": [{
                    "age_at_diagnosis": 22265, "days_to_last_follow_up": 420,
                    "primary_diagnosis": "Infiltrating duct carcinoma, NOS",
                    "tumor_stage": "stage iii", "tumor_grade": "G3",
                }],
                "exposures": [],
            },
        ],
        "pagination": {"total": 2}
    }
}

# A minimal STAR counts file with current GDC column layout
STAR_COUNTS_CONTENT = textwrap.dedent("""\
    gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second
    N_unmapped\t\t\t100\t50\t50
    N_multimapping\t\t\t200\t100\t100
    N_noFeature\t\t\t300\t150\t150
    N_ambiguous\t\t\t50\t25\t25
    ENSG00000000003\tTSPAN6\tprotein_coding\t1500\t750\t750
    ENSG00000000005\tTNMD\tprotein_coding\t0\t0\t0
    ENSG00000000419\tDPM1\tprotein_coding\t2300\t1150\t1150
    ENSG00000000457\tSCYL3\tprotein_coding\t800\t400\t400
    ENSG00000000460\tC1orf112\tprotein_coding\t1200\t600\t600
""")

# An older-format STAR file (no gene_name/gene_type columns)
STAR_COUNTS_OLD_FORMAT = textwrap.dedent("""\
    gene_id\tunstranded\tstranded_first\tstranded_second
    N_unmapped\t100\t50\t50
    ENSG00000000003\t1500\t750\t750
    ENSG00000000419\t2300\t1150\t1150
""")


# =============================================================================
#   Fixtures
# =============================================================================

@pytest.fixture
def tmp_output(tmp_path: Path) -> Path:
    """A clean temporary directory for each test."""
    return tmp_path


@pytest.fixture
def fake_star_file(tmp_path: Path) -> Path:
    """Write a fake STAR counts TSV and return its path."""
    p = tmp_path / "sample_001.star_counts.tsv"
    p.write_text(STAR_COUNTS_CONTENT, encoding="utf-8")
    return p


@pytest.fixture
def fake_star_file_old(tmp_path: Path) -> Path:
    """Write an older-format STAR counts TSV and return its path."""
    p = tmp_path / "sample_old.star_counts.tsv"
    p.write_text(STAR_COUNTS_OLD_FORMAT, encoding="utf-8")
    return p


@pytest.fixture
def fake_metadata_df() -> pd.DataFrame:
    """Minimal metadata DataFrame matching FAKE_FILES_RESPONSE."""
    from tcga_downloader.merge import flatten_file_metadata
    return flatten_file_metadata(FAKE_FILES_RESPONSE["data"]["hits"])


@pytest.fixture
def fake_clinical_df() -> pd.DataFrame:
    """Minimal clinical DataFrame matching FAKE_CASES_RESPONSE."""
    from tcga_downloader.merge import flatten_clinical
    return flatten_clinical(FAKE_CASES_RESPONSE["data"]["hits"])


@pytest.fixture
def fake_count_matrix(tmp_path: Path, fake_metadata_df: pd.DataFrame) -> pd.DataFrame:
    """A small 5-gene × 2-sample count matrix for merge tests."""
    # Create two TSV files in a UUID-subfolder structure (as GDC provides)
    genes = ["ENSG00000000003", "ENSG00000000005", "ENSG00000000419",
             "ENSG00000000457", "ENSG00000000460"]

    counts_dir = tmp_path / "raw_counts"
    for i, row in fake_metadata_df.iterrows():
        uuid_dir = counts_dir / f"uuid-{i:03d}"
        uuid_dir.mkdir(parents=True, exist_ok=True)
        lines = ["gene_id\tunstranded\tstranded_first\tstranded_second"]
        lines += [f"N_unmapped\t50\t25\t25"]
        for j, g in enumerate(genes):
            count = (j + 1) * 100 * (i + 1)
            lines.append(f"{g}\t{count}\t{count//2}\t{count//2}")
        (uuid_dir / row["file_name"]).write_text("\n".join(lines), encoding="utf-8")

    from tcga_downloader.matrix import build_count_matrix
    return build_count_matrix(counts_dir, fake_metadata_df, verbose=False)


@pytest.fixture
def mocked_gdc_api():
    """
    Activate the `responses` HTTP mock library for the duration of a test.
    Tests using this fixture will raise ConnectionError for any unregistered URL.
    """
    with responses_lib.RequestsMock(assert_all_requests_are_fired=False) as rsps:
        yield rsps


@pytest.fixture
def gdc_client_no_token():
    """A GDCClient with no authentication token."""
    from tcga_downloader.client import GDCClient
    return GDCClient()


# ── Integration test guard ────────────────────────────────────────────────────

def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: tests that make real GDC API calls (skipped unless TCGA_INTEGRATION=1)"
    )


def pytest_runtest_setup(item):
    if "integration" in item.keywords:
        if not os.environ.get("TCGA_INTEGRATION"):
            pytest.skip("Set TCGA_INTEGRATION=1 to run integration tests")
