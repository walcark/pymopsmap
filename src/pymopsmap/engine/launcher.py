"""MOPSMAP binary execution and output capture."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from pymopsmap.exceptions import MopsmapError
from pymopsmap.utils import MOPSMAP_PATH, get_logger

logger = get_logger(__name__)


def launch_mopsmap(
    input_filename: Path,
    expected_output: Path | None = None,
) -> dict[str, Any]:
    """
    Launch the MOPSMAP binary and return captured artefacts.

    Parameters
    ----------
    input_filename : Path
        Path to the MOPSMAP launch file.
    expected_output : Path, optional
        Expected path of the NetCDF output file. When not supplied, the
        function searches for ``output.nc`` in the same directory as the
        launch file.
    """
    cmd = [str(MOPSMAP_PATH), str(input_filename)]
    logger.debug("Running MOPSMAP: %s", " ".join(cmd))

    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)

    stdout = proc.stdout or ""
    stderr = proc.stderr or ""

    if stdout.strip():
        logger.debug("[MOPSMAP STDOUT]\n%s", stdout)
    if stderr.strip():
        logger.warning("[MOPSMAP STDERR]\n%s", stderr)

    if proc.returncode != 0:
        raise MopsmapError(
            f"MOPSMAP exited with code {proc.returncode}",
            returncode=proc.returncode,
            stderr=stderr,
        )

    output_path: Path | None
    if expected_output is not None:
        output_path = expected_output
    else:
        candidates = list(input_filename.parent.glob("output.nc"))
        output_path = candidates[0] if candidates else None

    if output_path is None or not output_path.exists():
        raise MopsmapError(
            "MOPSMAP exited with code 0 but produced no output file.",
            returncode=0,
            stderr=stderr,
        )

    logger.debug("MOPSMAP finished. Output: %s", output_path)
    return {"output_path": output_path, "stdout": stdout}
