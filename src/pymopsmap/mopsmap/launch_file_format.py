"""
mopsmap_launch_file.py

Author  : Kévin Walcarius
Date    : 2025-11-19
Version : 1.0
License : MIT
Summary : Tools to format the input .txt file used by Mopsmap.
"""

import os

from pymopsmap.utils import get_tempfile, get_logger
from pymopsmap.classes import MicroParameters
from .commands import microparams_command, wl_command
from pathlib import Path

logger = get_logger(__name__)


def write_launching_file(
    mp: MicroParameters | list[MicroParameters],
    n_angles: int = 2000,
    rh: float | None = None,
    mopsmap_data_path: Path | None = None,
) -> dict[str, Path]:
    """
    Generate a temporary launching file for MOPSMAP v1.0.

    This function assembles the complete configuration required by MOPSMAP:
        - the path of the Optical properties dataset,
        - the list of aerosol modes,
        - additionnal parameters (wavelength, rh, number of scattering angles, etc.)

    The resulting configuration file is written to a unique file in the system
    temporary directory.

    Parameters
    ----------
    mp: MicroParameters
        A list of aerosol micro-parameters.
    n_angles : int, optional
        Number of scattering angles (default: 2000).
    rh : float, optional
        Relative humidity in percent (0–100). Not included if ``None``.

    Returns
    -------
    dict[str, Path]
        Path to the Mopsmap Launching file and the output netcdf file.
    """

    logger.debug("Writing Mopsmap input file.")

    paths = _generate_paths()

    # File content
    file_prefix = _file_prefix(mopsmap_data_path)
    file_content = microparams_command(mp)
    file_suffix = _file_suffix(paths["netcdf"], n_angles, rh)

    # Final file content
    file_content: str = "\n".join([file_prefix, file_content, file_suffix])
    logger.debug("File content:\n", file_content)

    # Create file in /tmp
    with open(paths["mopsmap"], "w") as f:
        f.write(file_content)

    logger.debug(f"Mopsmap input file writen: {paths['mopsmap']}")

    return paths


def _generate_paths() -> dict[str, Path]:
    """
    Generates the path of the Mopsmap launch file and of the output
    netcdf file.
    """
    netcdf_path: Path = get_tempfile("output.nc")
    mopsmap_path: Path = get_tempfile("mopsmap.txt")

    return {"mopsmap": mopsmap_path, "netcdf": netcdf_path}


def _file_prefix(path_to_dataset: Path | None = None) -> str:
    """
    Write the prefix of the input .txt file used by Mopsmap.
    """

    if path_to_dataset is None:
        default_path = os.getenv("MOPSMAP_DATASET_PATH", "")
        assert default_path != "", (
            "Mosmap dataset path should be provided as input or "
            "with the MOPSMAP_DATASET_PATH environment variable."
        )
        path_to_dataset = Path(default_path)

    return f"scatlib '{path_to_dataset}'"


def _file_suffix(
    nc_path: Path,
    n_angles: int,
    rh: float | None = None,
) -> str:
    """
    Write the suffix of the input .txt file used by Mopsmap.
    """
    file_suffix = [f"output num_theta {n_angles}"]
    file_suffix.append(wl_command())
    if rh is not None:
        file_suffix.append("rH {rh}")
    file_suffix.append("output integrated")
    file_suffix.append(f"output netcdf '{str(nc_path)}'")

    return "\n".join(file_suffix)
