"""
Structured exception hierarchy for tcga_downloader.

All exceptions carry both a user-facing message and a suggested fix,
so they can be displayed cleanly in both the CLI and the GUI.
"""

from __future__ import annotations


class GDCError(Exception):
    """
    Base exception for all GDC-related errors.

    Parameters
    ----------
    message : str
        Plain-English description of what went wrong.
    fix : str
        Step-by-step instructions on how to resolve the issue.
    step : str
        Which pipeline step raised the error (for diagnostics).
    http_status : int | None
        HTTP status code if the error came from an API call.
    """

    def __init__(
        self,
        message: str,
        fix: str = "",
        step: str = "",
        http_status: int | None = None,
    ):
        self.message     = message
        self.fix         = fix
        self.step        = step
        self.http_status = http_status
        super().__init__(message)

    def __str__(self) -> str:
        parts = [self.message]
        if self.step:
            parts.append(f"(step: {self.step})")
        return " ".join(parts)

    def formatted(self) -> str:
        """Return a multi-line string suitable for printing to the user."""
        lines = [
            "─" * 60,
            f"ERROR: {self.message}",
        ]
        if self.http_status:
            lines.append(f"HTTP status: {self.http_status}")
        if self.fix:
            lines.append("")
            lines.append("How to fix it:")
            for line in self.fix.strip().splitlines():
                lines.append(f"  {line}")
        lines.append("─" * 60)
        return "\n".join(lines)


class AuthenticationError(GDCError):
    """Token is missing, expired, or rejected by the GDC."""

    def __init__(self, extra: str = ""):
        super().__init__(
            message="GDC authentication failed. Your token may be missing or expired.",
            fix=(
                "1. Go to https://portal.gdc.cancer.gov and log in.\n"
                "2. Click your username → 'Download Token'.\n"
                "3. Provide the new token file path when prompted.\n"
                "   (Tokens expire after 30 days.)"
                + (f"\n\nExtra detail: {extra}" if extra else "")
            ),
            step="authentication",
            http_status=401,
        )


class AccessDeniedError(GDCError):
    """File or project requires dbGaP approval."""

    def __init__(self, extra: str = ""):
        super().__init__(
            message="Access denied. This data requires dbGaP authorization.",
            fix=(
                "1. Check whether the files are open-access or controlled-access:\n"
                "   https://portal.gdc.cancer.gov\n"
                "2. For controlled-access, apply for dbGaP access:\n"
                "   https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?page=login\n"
                "3. Most TCGA STAR count files ARE open-access — check your filters."
                + (f"\n\nExtra detail: {extra}" if extra else "")
            ),
            step="download",
            http_status=403,
        )


class RateLimitError(GDCError):
    """GDC API rate limit hit."""

    def __init__(self):
        super().__init__(
            message="Too many requests — GDC API is rate-limiting this connection.",
            fix=(
                "1. Wait 60 seconds, then re-run the script.\n"
                "   (It will resume automatically from where it stopped.)\n"
                "2. Avoid running multiple downloads at the same time."
            ),
            step="api request",
            http_status=429,
        )


class ConnectionError(GDCError):  # noqa: A001
    """Network unreachable or GDC is down."""

    def __init__(self, context: str = ""):
        super().__init__(
            message=f"Cannot reach the GDC API{' during ' + context if context else ''}.",
            fix=(
                "1. Check your internet connection.\n"
                "2. Try opening https://api.gdc.cancer.gov/status in a browser.\n"
                "3. Check for GDC outages: https://status.gdc.cancer.gov\n"
                "4. Re-run the script when connectivity is restored."
            ),
            step=context,
        )


class TimeoutError(GDCError):  # noqa: A001
    """Request or download timed out."""

    def __init__(self, context: str = ""):
        super().__init__(
            message=f"Request timed out{' during ' + context if context else ''}.",
            fix=(
                "1. Re-run the script — it will resume from this step.\n"
                "2. For very large projects, consider the GDC Data Transfer Tool:\n"
                "   https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/\n"
                "3. Check that your internet connection is stable."
            ),
            step=context,
        )


class NoFilesFoundError(GDCError):
    """The GDC query returned zero files."""

    def __init__(self, project_id: str):
        super().__init__(
            message=f"No open-access STAR count files found for project '{project_id}'.",
            fix=(
                f"1. Verify '{project_id}' has RNA-Seq data:\n"
                "   https://portal.gdc.cancer.gov/repository\n"
                "   Filters: Data Type = Gene Expression Quantification,\n"
                "            Workflow = STAR - Counts, Access = open\n"
                "2. Some projects only have controlled-access RNA-seq data\n"
                "   (requires dbGaP approval).\n"
                "3. Double-check the project ID spelling (e.g. TCGA-BRCA)."
            ),
            step="file discovery",
        )


class DownloadCorruptedError(GDCError):
    """Downloaded file is not a valid archive."""

    def __init__(self, path: str):
        super().__init__(
            message="The downloaded archive is corrupted or incomplete.",
            fix=(
                f"1. Delete the corrupted file:\n   {path}\n"
                "2. Re-run the script — the download will restart.\n"
                "3. If this keeps happening, your connection may be dropping mid-download.\n"
                "   Try the GDC Data Transfer Tool for more reliable large downloads."
            ),
            step="zip extraction",
        )


class ColumnDetectionError(GDCError):
    """Cannot identify the unstranded count column in a STAR file."""

    def __init__(self, filepath: str, found_headers: list[str]):
        super().__init__(
            message=(
                f"Could not identify the 'unstranded' count column in:\n  {filepath}\n"
                f"  Found columns: {found_headers}"
            ),
            fix=(
                "1. Open the file in a text editor to inspect its columns.\n"
                "2. Expected: gene_id  unstranded  stranded_first  stranded_second\n"
                "3. If the column is named differently, report this at:\n"
                "   https://github.com/yourusername/tcga-gdc-downloader/issues\n"
                "4. Include the first 5 lines of the file in your report."
            ),
            step="column detection",
        )


class DiskSpaceError(GDCError):
    """Not enough disk space to write output."""

    def __init__(self, path: str, detail: str = ""):
        super().__init__(
            message=f"Could not write to disk at: {path}",
            fix=(
                "1. Check available disk space on your drive.\n"
                "2. Free up space or choose a different output directory.\n"
                "3. TCGA full projects can be 5–50 GB of raw files."
                + (f"\n\nSystem detail: {detail}" if detail else "")
            ),
            step="file write",
        )
