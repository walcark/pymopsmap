from pathlib import Path
from typing import List
import subprocess
import tempfile
import logging
import uuid
import os

from pymopsmap.classes import Wavelength, Mode
from .temp import get_tempdir

logger = logging.getLogger()


def netcdf_from_mopsmap(
    wl: Wavelength,
    modes: List[Mode],
    n_angles: int = 2000,
    rh: float | None = None,
) -> str:
    file = write_launching_file(wl=wl, modes=modes, n_angles=n_angles, rh=rh)
    launch_mopsmap(file)
    return get_tempdir() / "output.nc"


def launch_mopsmap(input_filename: str) -> None:
    """
    Launch MOPSMAP v1.0 using a prepared launching file.

    Parameters
    ----------
    input_filename : str or pathlib.Path
        Path to the launching file previously created with
        :func:`write_launching_file`.

    Raises
    ------
    RuntimeError
        If the MOPSMAP executable path is not defined.
    subprocess.CalledProcessError
        If MOPSMAP returns a non-zero exit code.
    """
    path_to_mopsmap: str = os.getenv("MOPSMAP_PATH")
    if path_to_mopsmap is None:
        raise RuntimeError(
            "Environment variable MOPSMAP_PATH is not set. "
            "It must point to a directory containing the `mopsmap` "
            "executable."
        )

    input_filename: Path = Path(tempfile.gettempdir()) / input_filename
    cmd: List[str] = [f"{path_to_mopsmap}/mopsmap", f"{input_filename}"]

    try:
        subprocess.run(cmd, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        logger.error("[ERROR] MOPSMAP v1.0 failed.")
        logger.error(f"[COMMAND] {' '.join(e.cmd)}")
        logger.error(f"[RETURN CODE] {e.returncode}")
        logger.error(f"[STDOUT]\n{e.stdout}")
        logger.error(f"[STDERR]\n{e.stderr}")
        raise


def write_launching_file(
    wl: Wavelength,
    modes: List[Mode],
    n_angles: int = 2000,
    rh: float | None = None,
) -> str:
    """
    Generate a temporary launching file for MOPSMAP v1.0.

    This function assembles the complete configuration required by MOPSMAP:
    the scattering library path, the list of aerosol modes, the number of
    scattering angles, wavelength configuration, and optional relative
    humidity specification. The resulting configuration file is written to a
    unique file in the system temporary directory.

    Parameters
    ----------
    wl : Wavelength
        A Wavelength object containing one or multiple wavelengths.
    modes : list of Mode
        A list of aerosol modes. Each mode provides its own MOPSMAP command
        block via ``Mode.command(i)``.
    n_angles : int, optional
        Number of scattering angles (default: 2000).
    rh : float, optional
        Relative humidity in percent (0–100). Not included if ``None``.

    Returns
    -------
    pathlib.Path
        Path to the generated launching file.

    Notes
    -----
    - The scattering dataset path and output paths are currently hardcoded
      and may need to be customized by the user.
    - The file is written into the system temporary directory and is not
      automatically deleted.
    """
    # Harcoded dataset path (#TODO should be override in the future)
    path_to_dataset = "/home/kwalcarius/bin/mopsmap/optical_dataset"
    file_prefix = f"scatlib '{path_to_dataset}'"

    # Modes
    file_content = "\n".join([mode.command(i + 1) for i, mode in enumerate(modes)])

    # Suffix after modes
    file_suffix = [f"output num_theta {n_angles}", wl.command]

    if rh:
        file_suffix += "rH {rh}"

    file_suffix.extend(
        [
            f"output netcdf '{str(get_tempdir())}/output.nc'",
            "output integrated",
        ]
    )

    # Final file content
    file_content: str = "\n".join([file_prefix, file_content, "\n".join(file_suffix)])
    print(file_content)

    # Create file in /tmp
    filepath = get_tempdir() / f"launching_file_{uuid.uuid4().hex}.txt"
    with open(filepath, "w") as f:
        f.write(file_content)

    return filepath


if __name__ == "__main__":
    from pymopsmap.classes import (
        SpheroidLognormal,
        wavelength,
        LognormalPSD,
        RefractiveIndex,
        Mode,
    )

    wl = wavelength([0.44, 0.55, 0.67])

    mode = Mode(
        shape=SpheroidLognormal(zeta1=0.5, zeta2=0.5, aspect_ratio=1.5, sigma_ar=2.0),
        psd=LognormalPSD(rm=0.2, sigma=1.5, n=0.1, rmin=0.01, rmax=2.0),
        refr_index=RefractiveIndex(
            wl=[0.44, 0.55, 0.67], n_real=[1.5, 1.6, 1.9], n_imag=[0.01, 0.001, 0.0001]
        ),
    )

    netcdf_from_mopsmap(wl=wl, modes=[mode])

    import xarray as xr

    ds = xr.load_dataset(f"{tempfile.gettempdir()}/output.nc")
    print(ds)
