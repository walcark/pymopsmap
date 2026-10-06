"""One directory per MOPSMAP run.

The launch file, the refractive index files and the ascii outputs all carry
fixed names, so two runs sharing a directory overwrite each other. A single
process-wide directory was enough while runs were sequential; a sweep that
computes several points at once is not.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path
from types import TracebackType

KEEP_TEMP = "PYMOPSMAP_KEEP_TEMP"


class Workspace:
    """A temporary directory owned by one run.

    Examples
    --------
    >>> with Workspace() as workspace:
    ...     launch = workspace.file("mopsmap.txt")
    """

    def __init__(self) -> None:
        self.path = Path(tempfile.mkdtemp(prefix="pymopsmap-"))

    def file(self, name: str) -> Path:
        """The path a named artefact of this run takes."""
        return self.path / name

    def __enter__(self) -> Workspace:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        if _keep_requested():
            return
        shutil.rmtree(self.path, ignore_errors=True)


def _keep_requested() -> bool:
    return os.environ.get(KEEP_TEMP, "").lower() in ("1", "true", "yes")
