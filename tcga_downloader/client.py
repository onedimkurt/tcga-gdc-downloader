"""
client.py — GDC REST API wrapper
=================================

DOWNLOAD METHOD USED
--------------------
This client uses the GDC REST API /data endpoint to download files by UUID:

    POST https://api.gdc.cancer.gov/data
    Content-Type: application/json
    Body: {"ids": ["file_id_uuid_1", "file_id_uuid_2", ...]}

The response is a streaming ZIP archive. Each file inside is placed in a
subdirectory named after its file_id UUID:
    <file_id_uuid>/<filename>.tsv

This is NOT the manifest-based gdc-client approach. We never:
  - Create a manifest.txt file
  - Call gdc-client
  - Pass a plain list of UUIDs as a manifest (that would fail with
    "Invalid manifest" because -m expects a multi-column TSV, not UUIDs)

WHY THIS IS THE CORRECT APPROACH
---------------------------------
The GDC docs explicitly document two download methods:
  1. REST API /data (what we use) — accepts {"ids": [...]} JSON body
  2. gdc-client -m manifest.txt   — requires a full manifest TSV file

For programmatic downloads without the gdc-client tool installed,
method 1 is the correct and recommended approach.
See: https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/

UUID VALIDATION
---------------
All file_ids passed to stream_download() are validated against the
GDC UUID format before the request is sent. This prevents:
  - Accidentally passing case_ids or sample_ids instead of file_ids
  - Passing submitter barcodes (TCGA-XX-XXXX) instead of UUIDs
  - Passing empty strings or None values
"""

from __future__ import annotations
import time
import requests

from tcga_downloader.constants import (
    GDC_FILES_ENDPOINT, GDC_CASES_ENDPOINT, GDC_DATA_ENDPOINT,
    GDC_PROJECTS_ENDPOINT, GDC_STATUS_ENDPOINT,
    FILE_FIELDS, CLINICAL_FIELDS, STAR_COUNTS_FILTERS,
    REQUEST_TIMEOUT_SHORT, REQUEST_TIMEOUT_LONG,
    MAX_RETRIES, RETRY_WAIT_SECONDS, CHUNK_SIZE_BYTES,
    GDC_PAGE_SIZE, MAX_TOTAL_FILES,
    is_valid_gdc_uuid,
)
from tcga_downloader.exceptions import (
    AuthenticationError, AccessDeniedError, RateLimitError,
    NoFilesFoundError, GDCError,
)
import tcga_downloader.exceptions as exc_mod


class GDCClient:
    """
    Wrapper around the GDC REST API.

    Usage
    -----
    client = GDCClient()                          # open-access, no token
    client = GDCClient(token="your-token-string") # with authentication
    client = GDCClient.from_file("/path/to/gdc-user-token.txt")
    """

    def __init__(self, token: str | None = None):
        self._session = requests.Session()
        self._session.headers.update({"Content-Type": "application/json"})
        if token:
            self._session.headers.update({"X-Auth-Token": token.strip()})

    # ── Constructors ──────────────────────────────────────────────────────────

    @classmethod
    def from_file(cls, token_path: str) -> "GDCClient":
        """Load a GDC authentication token from a file."""
        from pathlib import Path
        p = Path(token_path)
        if not p.exists():
            raise GDCError(
                f"Token file not found: {token_path}",
                fix=(
                    "1. Check the path is correct.\n"
                    "2. Download a token at https://portal.gdc.cancer.gov\n"
                    "   (Login → click username → Download Token)"
                ),
                step="authentication",
            )
        token = p.read_text(encoding="utf-8").strip()
        if len(token) < 20:
            raise GDCError(
                "Token file appears empty or too short.",
                fix="Download a fresh token at https://portal.gdc.cancer.gov",
                step="authentication",
            )
        return cls(token=token)

    # ── Low-level HTTP ────────────────────────────────────────────────────────

    def _post(self, endpoint: str, payload: dict, context: str,
              stream: bool = False, timeout: int = REQUEST_TIMEOUT_SHORT) -> requests.Response:
        """POST with retry logic. Translates all errors to structured GDCError."""
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                resp = self._session.post(endpoint, json=payload,
                                          stream=stream, timeout=timeout)
                resp.raise_for_status()
                return resp

            except requests.exceptions.Timeout:
                if attempt < MAX_RETRIES:
                    print(f"  ⚠️  Timeout (attempt {attempt}/{MAX_RETRIES}), "
                          f"retrying in {RETRY_WAIT_SECONDS}s...")
                    time.sleep(RETRY_WAIT_SECONDS)
                else:
                    raise exc_mod.TimeoutError(context)

            except requests.exceptions.ConnectionError:
                raise exc_mod.ConnectionError(context)

            except requests.exceptions.HTTPError:
                status = resp.status_code
                body = ""
                try:
                    body = resp.text[:300]
                except Exception:
                    pass
                if status == 401: raise AuthenticationError(body)
                if status == 403: raise AccessDeniedError(body)
                if status == 429:
                    if attempt < MAX_RETRIES:
                        print(f"  ⚠️  Rate limited. Waiting 60s before retry...")
                        time.sleep(60)
                        continue
                    raise RateLimitError()
                if status == 500:
                    raise GDCError(
                        f"GDC server error (HTTP 500) during: {context}",
                        fix=(
                            "This is a server-side error at the GDC.\n"
                            "1. Wait a few minutes and re-run.\n"
                            "2. If it persists: https://gdc.cancer.gov/support"
                        ),
                        step=context, http_status=500,
                    )
                raise GDCError(f"HTTP {status} during {context}: {body}",
                               fix="Re-run. If repeated, check https://status.gdc.cancer.gov",
                               step=context, http_status=status)

        raise GDCError(f"Failed after {MAX_RETRIES} attempts: {context}", step=context)

    # ── Public API ────────────────────────────────────────────────────────────

    def check_connectivity(self) -> bool:
        """Return True if the GDC API is reachable."""
        try:
            return self._session.get(GDC_STATUS_ENDPOINT, timeout=10).status_code == 200
        except Exception:
            return False

    def validate_project(self, project_id: str) -> dict:
        """Verify a project ID exists. Returns project info dict or raises GDCError."""
        resp = self._post(
            GDC_PROJECTS_ENDPOINT,
            {"filters": {"op": "=", "content": {"field": "project_id", "value": project_id}},
             "fields": "project_id,name,primary_site,disease_type,summary,program.name",
             "size": 1},
            context="project validation",
        )
        hits = resp.json().get("data", {}).get("hits", [])
        if not hits:
            raise GDCError(
                f"Project '{project_id}' not found in the GDC.",
                fix=(
                    "1. Check spelling — IDs are case-sensitive: TCGA-BRCA\n"
                    "2. Browse all projects: https://portal.gdc.cancer.gov/projects"
                ),
                step="project validation",
            )
        return hits[0]

    def discover_star_files(self, project_id: str) -> list[dict]:
        """
        Return metadata for all open-access STAR-Counts TSV files in the project.

        Uses the STAR_COUNTS_FILTERS from constants.py to ensure ONLY
        STAR unstranded count files are returned (not FPKM, miRNA, WGS, etc.).

        Paginates automatically if the project has more than GDC_PAGE_SIZE files.

        Each returned dict contains:
            file_id   : UUID used for /data download (NOT a manifest UUID list)
            file_name : actual TSV filename
            file_size : bytes
            cases     : list with case_id, submitter_id, samples[...]

        Raises NoFilesFoundError if zero files match.
        """
        project_filter = {"op": "=", "content": {
            "field": "cases.project.project_id", "value": project_id
        }}
        all_filters = {"op": "and", "content": [project_filter] + STAR_COUNTS_FILTERS}

        payload = {
            "filters": all_filters,
            "fields": FILE_FIELDS,
            "size": GDC_PAGE_SIZE,
            "from": 0,
            "format": "json",
        }

        all_hits: list[dict] = []
        page = 0

        while True:
            payload["from"] = page * GDC_PAGE_SIZE
            resp = self._post(GDC_FILES_ENDPOINT, payload, context="file discovery")
            data = resp.json().get("data", {})
            hits = data.get("hits", [])
            total = data.get("pagination", {}).get("total", len(hits))

            all_hits.extend(hits)

            if len(all_hits) >= total:
                break
            if len(all_hits) >= MAX_TOTAL_FILES:
                print(f"  ⚠️  Safety limit reached ({MAX_TOTAL_FILES} files). "
                      f"Proceeding with {len(all_hits)} files.")
                break

            page += 1
            print(f"  📄  Retrieved {len(all_hits)}/{total} file records...")

        if not all_hits:
            raise NoFilesFoundError(project_id)

        # ── UUID validation: every file_id must be a valid GDC UUID ──────────
        invalid = [h.get("file_id", "") for h in all_hits
                   if not is_valid_gdc_uuid(h.get("file_id", ""))]
        if invalid:
            raise GDCError(
                f"GDC returned {len(invalid)} records with invalid file_id format:\n"
                + "\n".join(f"  '{v}'" for v in invalid[:5]),
                fix=(
                    "This is unexpected and likely indicates a GDC API change.\n"
                    "Please report this at the project GitHub issues page,\n"
                    "including the project ID and the invalid IDs shown above."
                ),
                step="file discovery UUID validation",
            )

        return all_hits

    def fetch_clinical_data(self, project_id: str) -> list[dict]:
        """Return raw clinical records for all cases in the project."""
        resp = self._post(
            GDC_CASES_ENDPOINT,
            {"filters": {"op": "=", "content": {"field": "project.project_id",
                                                  "value": project_id}},
             "fields": CLINICAL_FIELDS,
             "size": GDC_PAGE_SIZE,
             "format": "json"},
            context="clinical data fetch",
        )
        return resp.json().get("data", {}).get("hits", [])

    def stream_download(self, file_ids: list[str], dest_path: str) -> int:
        """
        Download files by UUID using the GDC REST API /data endpoint.

        METHOD: POST /data with {"ids": [file_id_uuid, ...]}
        This is the correct UUID-direct download method.
        It is NOT the gdc-client manifest method.

        The IDs must be file_id UUIDs (the 'file_id' field from /files queries).
        Do NOT pass:
          - case_id UUIDs
          - sample_id UUIDs
          - TCGA submitter barcodes (TCGA-XX-XXXX)
          - A plain text file of UUIDs (that's the gdc-client -m format)

        Parameters
        ----------
        file_ids : list[str]
            List of file_id UUIDs from discover_star_files().
        dest_path : str
            Where to write the downloaded ZIP archive.

        Returns
        -------
        int
            Number of bytes written.
        """
        if not file_ids:
            raise GDCError("No file IDs provided to download.", step="download")

        # ── Validate every UUID before sending to avoid silent failures ───────
        bad = [fid for fid in file_ids if not is_valid_gdc_uuid(fid)]
        if bad:
            raise GDCError(
                f"{len(bad)} invalid file_id value(s) detected before download:\n"
                + "\n".join(f"  '{v}'" for v in bad[:5]),
                fix=(
                    "Only file_id UUIDs (from the GDC /files API) should be passed\n"
                    "to stream_download(). Do not pass case_ids, sample_ids, or\n"
                    "TCGA submitter barcodes. Format: xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
                ),
                step="download UUID validation",
            )

        from pathlib import Path
        resp = self._post(
            GDC_DATA_ENDPOINT,
            {"ids": file_ids},
            context="bulk download",
            stream=True,
            timeout=REQUEST_TIMEOUT_LONG,
        )

        dest = Path(dest_path)
        total_size = int(resp.headers.get("content-length", 0))
        written = 0

        try:
            from tqdm import tqdm
            pbar = tqdm(total=total_size or None, unit="B", unit_scale=True,
                        unit_divisor=1024, desc="  Downloading", ncols=70)
        except ImportError:
            pbar = _DummyProgress()

        max_attempts = MAX_RETRIES + 1
        last_error = None

        for attempt in range(1, max_attempts + 1):
            try:
                # Re-request on retry attempts
                if attempt > 1:
                    import time
                    time.sleep(5 * attempt)
                    resp = self._post(
                        GDC_DATA_ENDPOINT,
                        {"ids": file_ids},
                        context="bulk download",
                        stream=True,
                        timeout=REQUEST_TIMEOUT_LONG,
                    )
                    written = 0
                    if hasattr(pbar, "reset"):
                        pbar.reset()

                with open(dest, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=CHUNK_SIZE_BYTES):
                        if chunk:
                            f.write(chunk)
                            written += len(chunk)
                            if hasattr(pbar, "update"):
                                pbar.update(len(chunk))
                break  # success — exit retry loop

            except (OSError, IOError) as e:
                raise exc_mod.DiskSpaceError(str(dest), str(e))
            except Exception as e:
                last_error = e
                if attempt < max_attempts:
                    print(f"\n  ⚠️  Download interrupted (attempt {attempt}/{max_attempts}), retrying...")
                    if dest.exists():
                        dest.unlink()
                else:
                    raise GDCError(
                        f"Download failed after {max_attempts} attempts: {e}",
                        fix="Check your internet connection and try again.",
                        step="download",
                    )
        if hasattr(pbar, "close"):
            pbar.close()
        return written


class _DummyProgress:
    def __enter__(self): return self
    def __exit__(self, *_): pass
    def update(self, n): pass
