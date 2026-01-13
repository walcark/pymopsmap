"""
mopsmap.py

Author  : Kévin Walcarius
Date    : 2025-11-19
Version : 1.0
License : MIT
Summary : Module used to define the Mopsmap launch command.
"""

import numpy as np
from typing import Any
import os
import subprocess
from pathlib import Path

from pymopsmap.utils import get_logger, MOPSMAP_PATH

logger = get_logger(__name__)


def launch_mopsmap(input_filename: Path) -> dict[str, Any]:
    """
    Launches the Mopsmap executable given a Mopsmap input file
    with all the execution directives. The filename is formatted
    by the `write_launching_file` method.
    """

    cmd = [str(MOPSMAP_PATH), str(input_filename)]
    logger.debug(f"Running Mopsmap command: {' '.join(cmd)}")

    # Run process and capture streams
    proc = subprocess.run(
        cmd,
        text=True,
        capture_output=True,  # keep capture
        check=False,  # we handle error ourselves
    )

    stdout = proc.stdout or ""
    stderr = proc.stderr or ""

    # Log everything in debug mode
    if stdout.strip():
        logger.debug("[MOPSMAP STDOUT]\n" + stdout)

    if stderr.strip():
        logger.warning("[MOPSMAP STDERR]\n" + stderr)

    # Detect process failure
    if proc.returncode != 0:
        raise RuntimeError(
            f"MOPSMAP exited with code {proc.returncode}\n"
            f"STDERR:\n{stderr}\n"
            f"STDOUT:\n{stdout}"
        )

    # Detect silent failure when output.nc is not in the /tmp/... directory.
    output_dir = input_filename.parent
    output = list(output_dir.glob("output.nc"))

    if not output:
        raise RuntimeError(
            "MOPSMAP finished with returncode 0 but did NOT produce "
            "output.nc. This indicates a silent failure. "
            "Check STDERR above."
        )

    logger.debug("Mopsmap finished successfully.")
    logger.debug(f"Output file: {output[0]}")

    return {"output_path": output[0], "stdout": stdout}
