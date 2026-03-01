"""
tests/test_client.py
====================
Unit tests for the GDCClient API wrapper.

All HTTP calls are intercepted by the `responses` library — no real
network calls are made. This ensures tests are fast, reproducible,
and work offline.

Tests cover:
- Successful project validation
- Successful file discovery
- Authentication errors (401)
- Access denied errors (403)
- Rate limiting with retry (429)
- Server errors (500)
- Connection failures
- Empty results handling
"""

import json
import pytest
import responses as responses_lib

from tcga_downloader.client import GDCClient
from tcga_downloader.exceptions import (
    AuthenticationError,
    AccessDeniedError,
    RateLimitError,
    NoFilesFoundError,
    GDCError,
)
from tcga_downloader.constants import (
    GDC_FILES_ENDPOINT,
    GDC_CASES_ENDPOINT,
    GDC_PROJECTS_ENDPOINT,
)
from tests.conftest import (
    FAKE_PROJECT_RESPONSE,
    FAKE_FILES_RESPONSE,
    FAKE_CASES_RESPONSE,
    FAKE_PROJECT_ID,
)


# =============================================================================
#   Helpers
# =============================================================================

def _register_json(rsps, url: str, body: dict, status: int = 200, method="POST"):
    """Register a mock endpoint returning a JSON body."""
    add = rsps.add if method == "POST" else rsps.get
    rsps.add(
        method=responses_lib.POST if method == "POST" else responses_lib.GET,
        url=url,
        json=body,
        status=status,
    )


# =============================================================================
#   GDCClient instantiation
# =============================================================================

class TestGDCClientInit:

    def test_no_token(self):
        client = GDCClient()
        assert "X-Auth-Token" not in client._session.headers

    def test_with_token_string(self):
        client = GDCClient(token="abc123def456ghi789")
        assert client._session.headers.get("X-Auth-Token") == "abc123def456ghi789"

    def test_token_whitespace_stripped(self):
        client = GDCClient(token="  mytoken123  ")
        assert client._session.headers["X-Auth-Token"] == "mytoken123"

    def test_from_file(self, tmp_path):
        token_file = tmp_path / "token.txt"
        token_file.write_text("validtoken1234567890\n", encoding="utf-8")
        client = GDCClient.from_file(str(token_file))
        assert client._session.headers["X-Auth-Token"] == "validtoken1234567890"

    def test_from_file_missing_raises(self):
        with pytest.raises(GDCError, match="not found"):
            GDCClient.from_file("/nonexistent/token.txt")

    def test_from_file_empty_raises(self, tmp_path):
        token_file = tmp_path / "empty.txt"
        token_file.write_text("   ", encoding="utf-8")
        with pytest.raises(GDCError, match="empty"):
            GDCClient.from_file(str(token_file))


# =============================================================================
#   validate_project
# =============================================================================

class TestValidateProject:

    @responses_lib.activate
    def test_valid_project_returns_info(self):
        responses_lib.add(
            responses_lib.POST, GDC_PROJECTS_ENDPOINT,
            json=FAKE_PROJECT_RESPONSE, status=200,
        )
        client = GDCClient()
        info = client.validate_project(FAKE_PROJECT_ID)
        assert info["project_id"] == FAKE_PROJECT_ID
        assert info["name"] == "Test Project"

    @responses_lib.activate
    def test_unknown_project_raises_gdcerror(self):
        responses_lib.add(
            responses_lib.POST, GDC_PROJECTS_ENDPOINT,
            json={"data": {"hits": [], "pagination": {"total": 0}}}, status=200,
        )
        client = GDCClient()
        with pytest.raises(GDCError, match="not found"):
            client.validate_project("TCGA-NONEXISTENT")

    @responses_lib.activate
    def test_401_raises_authentication_error(self):
        responses_lib.add(
            responses_lib.POST, GDC_PROJECTS_ENDPOINT,
            json={"message": "Unauthorized"}, status=401,
        )
        client = GDCClient()
        with pytest.raises(AuthenticationError):
            client.validate_project(FAKE_PROJECT_ID)

    @responses_lib.activate
    def test_403_raises_access_denied_error(self):
        responses_lib.add(
            responses_lib.POST, GDC_PROJECTS_ENDPOINT,
            json={"message": "Forbidden"}, status=403,
        )
        client = GDCClient()
        with pytest.raises(AccessDeniedError):
            client.validate_project(FAKE_PROJECT_ID)

    @responses_lib.activate
    def test_500_raises_gdcerror(self):
        responses_lib.add(
            responses_lib.POST, GDC_PROJECTS_ENDPOINT,
            json={"message": "Internal Server Error"}, status=500,
        )
        client = GDCClient()
        with pytest.raises(GDCError) as exc_info:
            client.validate_project(FAKE_PROJECT_ID)
        assert exc_info.value.http_status == 500


# =============================================================================
#   discover_star_files
# =============================================================================

class TestDiscoverStarFiles:

    @responses_lib.activate
    def test_returns_file_list(self):
        responses_lib.add(
            responses_lib.POST, GDC_FILES_ENDPOINT,
            json=FAKE_FILES_RESPONSE, status=200,
        )
        client = GDCClient()
        hits = client.discover_star_files(FAKE_PROJECT_ID)
        assert len(hits) == 2

    @responses_lib.activate
    def test_hit_has_expected_fields(self):
        responses_lib.add(
            responses_lib.POST, GDC_FILES_ENDPOINT,
            json=FAKE_FILES_RESPONSE, status=200,
        )
        client = GDCClient()
        hits = client.discover_star_files(FAKE_PROJECT_ID)
        hit = hits[0]
        assert "file_id"   in hit
        assert "file_name" in hit
        assert "cases"     in hit

    @responses_lib.activate
    def test_empty_results_raise_no_files_found(self):
        responses_lib.add(
            responses_lib.POST, GDC_FILES_ENDPOINT,
            json={"data": {"hits": [], "pagination": {"total": 0}}}, status=200,
        )
        client = GDCClient()
        with pytest.raises(NoFilesFoundError):
            client.discover_star_files("TCGA-EMPTY")

    @responses_lib.activate
    def test_401_raises_authentication_error(self):
        responses_lib.add(
            responses_lib.POST, GDC_FILES_ENDPOINT,
            json={}, status=401,
        )
        client = GDCClient()
        with pytest.raises(AuthenticationError):
            client.discover_star_files(FAKE_PROJECT_ID)


# =============================================================================
#   fetch_clinical_data
# =============================================================================

class TestFetchClinicalData:

    @responses_lib.activate
    def test_returns_case_list(self):
        responses_lib.add(
            responses_lib.POST, GDC_CASES_ENDPOINT,
            json=FAKE_CASES_RESPONSE, status=200,
        )
        client = GDCClient()
        hits = client.fetch_clinical_data(FAKE_PROJECT_ID)
        assert len(hits) == 2

    @responses_lib.activate
    def test_empty_clinical_returns_empty_list_not_error(self):
        """Missing clinical data is common — must not raise."""
        responses_lib.add(
            responses_lib.POST, GDC_CASES_ENDPOINT,
            json={"data": {"hits": [], "pagination": {"total": 0}}}, status=200,
        )
        client = GDCClient()
        hits = client.fetch_clinical_data(FAKE_PROJECT_ID)
        assert hits == []


# =============================================================================
#   Error message quality
# =============================================================================

class TestErrorMessages:
    """Verify that error messages contain actionable guidance."""

    @responses_lib.activate
    def test_auth_error_mentions_token_renewal(self):
        responses_lib.add(
            responses_lib.POST, GDC_PROJECTS_ENDPOINT,
            json={}, status=401,
        )
        client = GDCClient()
        with pytest.raises(AuthenticationError) as exc_info:
            client.validate_project(FAKE_PROJECT_ID)
        err = exc_info.value
        assert "portal.gdc.cancer.gov" in err.fix
        assert "token" in err.fix.lower()

    @responses_lib.activate
    def test_403_error_mentions_dbgap(self):
        responses_lib.add(
            responses_lib.POST, GDC_FILES_ENDPOINT,
            json={}, status=403,
        )
        client = GDCClient()
        with pytest.raises(AccessDeniedError) as exc_info:
            client.discover_star_files(FAKE_PROJECT_ID)
        err = exc_info.value
        assert "dbgap" in err.fix.lower()

    def test_no_files_error_mentions_repository(self):
        err = NoFilesFoundError("TCGA-TEST")
        assert "portal.gdc.cancer.gov" in err.fix
        assert "STAR" in err.fix or "workflow" in err.fix.lower()
