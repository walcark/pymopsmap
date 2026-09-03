"""MOPSMAP launch file generation."""

from __future__ import annotations

from pathlib import Path

from pymopsmap.engine.outputs import DEFAULT_OUTPUT, OutputRequest, OutputType
from pymopsmap.microparams import MicroParameters
from pymopsmap.utils import DATASET_CACHE_DIR, MOPSMAP_PATH, get_logger

from .commands import microparams_command, wl_command
from .workspace import Workspace

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
    workspace: Workspace | None = None,
    output_types: OutputRequest = DEFAULT_OUTPUT,
    n_angles: int = 2000,
    rh: float | None = None,
    mopsmap_data_path: Path | None = None,
) -> dict[str, Path]:
    """
    Generate a MOPSMAP launch file and return paths to the generated artefacts.

    Returns a dict with:
      - mopsmap   : path to the launch .txt file
      - ascii_base: base path for ASCII output files
    """
    logger.debug("Writing MOPSMAP input file.")

    workspace = workspace or Workspace()
    paths = _generate_paths(workspace)

    dataset_path = mopsmap_data_path or DATASET_CACHE_DIR

    mp_list = [mp] if isinstance(mp, MicroParameters) else mp

    water_refr = MOPSMAP_PATH.parent / "data" / "refr_water_segelstein"
    file_prefix = f"scatlib '{dataset_path}'\nwater_refrac_file '{water_refr}'"
    file_content = microparams_command(mp, workspace)
    file_suffix = _file_suffix(
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


def _generate_paths(workspace: Workspace) -> dict[str, Path]:
    return {
        "mopsmap": workspace.file("mopsmap.txt"),
        "ascii_base": workspace.file("mopsmap_out"),
    }


def _file_suffix(
    ascii_base: Path | None,
    output_types: OutputRequest,
    n_angles: int,
    rh: float | None,
    wavelengths: list | None = None,
) -> str:
    lines = [f"output num_theta {n_angles}", wl_command(wavelengths)]
    if rh is not None:
        lines.append(f"rH {rh}")

    # No netcdf output is requested: nothing parses it, the results are read
    # from stdout and the ascii files. Asking for it also segfaults a binary
    # built against a different netcdf-fortran, after the computation has
    # already succeeded.
    lines.append("output integrated")

    ascii_needed = output_types & _ASCII_TYPES
    if ascii_needed and ascii_base is not None:
        lines.append(f"output ascii_file '{ascii_base}'")
        for otype in sorted(ascii_needed, key=lambda x: x.value):
            lines.append(f"output {otype.value}")

    return "\n".join(lines)
