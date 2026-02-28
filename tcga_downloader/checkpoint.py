"""
checkpoint.py
=============
Saves and restores pipeline progress so any step can be safely restarted
after a network failure, crash, or user interruption.

State is stored as a simple JSON file inside the project output directory.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from tcga_downloader.constants import CHECKPOINT_FILE


class Checkpoint:
    """
    Persistent pipeline state tracker.

    Each pipeline step calls ``save(step_name)`` when it completes.
    On restart, ``is_done(step_name)`` is consulted before doing work.

    Parameters
    ----------
    output_dir : Path
        Directory where checkpoint JSON will be written.
    """

    # Canonical step order — used to reset from a given step
    STEPS = [
        "auth",
        "project",
        "files_found",
        "downloaded",
        "clinical",
        "cdr",
        "matrix",
        "merged",
    ]

    def __init__(self, output_dir: Path):
        self.path = output_dir / CHECKPOINT_FILE
        self._data: dict = self._load()

    # ── Persistence ───────────────────────────────────────────────────────────

    def _load(self) -> dict:
        if self.path.exists():
            try:
                return json.loads(self.path.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                return {}
        return {}

    def _write(self) -> None:
        self.path.write_text(
            json.dumps(self._data, indent=2, default=str),
            encoding="utf-8",
        )

    # ── Public interface ──────────────────────────────────────────────────────

    def save(self, step: str, payload: dict | None = None) -> None:
        """
        Mark a step as completed and store optional associated data.

        Parameters
        ----------
        step : str
            One of the values in STEPS.
        payload : dict, optional
            Any serialisable data to store alongside the completion record.
            Useful for caching file paths, counts, etc.
        """
        self._data[step] = {
            "done":    True,
            "ts":      datetime.now().isoformat(),
            "payload": payload or {},
        }
        self._write()

    def is_done(self, step: str) -> bool:
        """Return True if the step has been marked as complete."""
        return bool(self._data.get(step, {}).get("done"))

    def get(self, step: str) -> dict:
        """Return the payload dict stored for a completed step."""
        return self._data.get(step, {}).get("payload", {})

    def reset_from(self, step: str) -> None:
        """
        Remove completion records for step and all steps that follow it.
        Used when the user chooses to redo part of the pipeline.

        Parameters
        ----------
        step : str
            The step to reset from (inclusive).
        """
        if step not in self.STEPS:
            return
        idx = self.STEPS.index(step)
        for s in self.STEPS[idx:]:
            self._data.pop(s, None)
        self._write()

    def any_done(self) -> bool:
        """Return True if at least one step has been completed."""
        return any(self.is_done(s) for s in self.STEPS)

    def summary(self) -> str:
        """
        Return a multi-line string summarising which steps are done.

        Example output:
            auth            ✅  done  (completed 2025-06-01 10:32:14)
            project         ✅  done  (completed 2025-06-01 10:32:20)
            files_found     ✅  done  (completed 2025-06-01 10:32:35)
            downloaded      ⬜  pending
        """
        lines = ["  Pipeline step status:"]
        for step in self.STEPS:
            rec = self._data.get(step, {})
            done = rec.get("done", False)
            ts   = rec.get("ts", "")
            ts_str = f"  (completed {ts[:19].replace('T', ' ')})" if ts else ""
            icon  = "✅" if done else "⬜"
            status = "done" if done else "pending"
            lines.append(f"    {step:<16} {icon}  {status}{ts_str}")
        return "\n".join(lines)
