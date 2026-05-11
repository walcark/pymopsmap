"""MOPSMAP launch file generation."""

from __future__ import annotations

from pathlib import Path

from pymopsmap.models import MicroParameters
from pymopsmap.models.output_request import (
    DEFAULT_OUTPUT,
    OutputRequest,
    OutputType,
)
from pymopsmap.utils import (
    DATASET_CACHE_DIR,
    MOPSMAP_PATH,
    get_logger,
    get_tempfile,
)

from .commands import microparams_command, wl_command

logger = get_logger(__name__)

_ASCII_TYPES = {
    OutputType.PHASE_FUNCTION,
    OutputType.SCATTERING_MATRIX,
    OutputType.VOLUME_SCATTERING_FUNCTION,
    OutputType.LIDAR,
    OutputType.COEFF,
}


def write_launching_file(
    mp: MicroParameters | list[MicroParameters],
    output_types: OutputRequest = DEFAULT_OUTPUT,
    n_angles: int = 2000,
    rh: float | None = None,
    mopsmap_data_path: Path | None = None,
) -> dict[str, Path]:
    """
    Generate a MOPSMAP launch file and return paths to the generated artefacts.

    Returns a dict with:
      - mopsmap   : path to the launch .txt file
      - netcdf    : path to the expected NetCDF output
      - ascii_base: base path for ASCII output files
    """
    logger.debug("Writing MOPSMAP input file.")

    paths = _generate_paths()

    dataset_path = mopsmap_data_path or DATASET_CACHE_DIR

    mp_list = [mp] if isinstance(mp, MicroParameters) else mp

    water_refr = MOPSMAP_PATH.parent / "data" / "refr_water_segelstein"
    file_prefix = f"scatlib '{dataset_path}'\nwater_refrac_file '{water_refr}'"
    file_content = microparams_command(mp)
    file_suffix = _file_suffix(
        nc_path=paths["netcdf"],
        ascii_base=paths.get("ascii_base"),
        output_types=output_types,
        n_angles=n_angles,
        rh=rh,
        wavelengths=mp_list[0].wavelength,
    )

    content = "\n".join([file_prefix, file_content, file_suffix])

    with open(paths["mopsmap"], "w") as f:
        f.write(content)

    logger.debug("MOPSMAP input file written: %s", paths["mopsmap"])
    return paths


def _generate_paths() -> dict[str, Path]:
    paths: dict[str, Path] = {
        "netcdf": get_tempfile("output.nc"),
        "mopsmap": get_tempfile("mopsmap.txt"),
    }
    paths["ascii_base"] = get_tempfile("mopsmap_out")
    return paths


def _file_suffix(
    nc_path: Path,
    ascii_base: Path | None,
    output_types: OutputRequest,
    n_angles: int,
    rh: float | None,
    wavelengths: list | None = None,
) -> str:
    lines = [f"output num_theta {n_angles}", wl_command(wavelengths)]
    if rh is not None:
        lines.append(f"rH {rh}")

    lines.append("output integrated")
    lines.append(f"output netcdf '{nc_path}'")

    ascii_needed = output_types & _ASCII_TYPES
    if ascii_needed and ascii_base is not None:
        lines.append(f"output ascii_file '{ascii_base}'")
        for otype in sorted(ascii_needed, key=lambda x: x.value):
            lines.append(f"output {otype.value}")

    return "\n".join(lines)
