"""MOPSMAP binary execution and output capture."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Any

from pymopsmap.exceptions import MopsmapError
from pymopsmap.utils import MOPSMAP_PATH, get_logger

logger = get_logger(__name__)

# MOPSMAP reports several failures on stdout and still exits zero: a dataset
# file it cannot open, and every Fortran `stop` that prints its reason first.
# The exit code alone therefore says nothing, and a partially failed run
# otherwise returns a mixture of correct values and silent gaps.
_FAILURE_MARKERS = ("Error opening", "Error:")


def _runtime_env() -> dict[str, str]:
    """
    The environment the MOPSMAP binary needs.

    It is linked against the NetCDF Fortran bindings, which ship in the
    project environment rather than on the system, so the loader has to be
    told where to look.
    """
    env = dict(os.environ)
    library_dir = Path(sys.prefix) / "lib"
    existing = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = (
        f"{library_dir}{os.pathsep}{existing}"
        if existing
        else str(library_dir)
    )
    return env


def launch_mopsmap(input_filename: Path) -> dict[str, Any]:
    """
    Launch the MOPSMAP binary and return its captured output.

    Parameters
    ----------
    input_filename : Path
        Path to the MOPSMAP launch file.

    Returns
    -------
    dict
        The captured stdout, which carries the integrated results.
    """
    cmd = [str(MOPSMAP_PATH), str(input_filename)]
    logger.debug("Running MOPSMAP: %s", " ".join(cmd))

    proc = subprocess.run(
        cmd, text=True, capture_output=True, check=False, env=_runtime_env()
    )

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

    reported = _reported_failures(stdout)
    if reported:
        raise MopsmapError(
            "MOPSMAP exited with code 0 but reported: " + " | ".join(reported),
            returncode=0,
            stderr=stderr,
        )

    logger.debug("MOPSMAP finished.")
    return {"stdout": stdout}


def _reported_failures(stdout: str) -> list[str]:
    """Lines on which MOPSMAP announced a failure without saying so in code."""
    return [
        line.strip()
        for line in stdout.splitlines()
        if any(marker in line for marker in _FAILURE_MARKERS)
    ]
