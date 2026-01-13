"""
optiprops_to_smartg.py

Author  : Kévin Walcarius
Date    : 2025-01-08
Version : 1.0
License : MIT
Summary : Module to transform optical parameters derived with
          Mopsmap into an xarray dataset readable by the Smart-G
          AerOPAC class.
"""

from einops import rearrange
from datetime import datetime
from pathlib import Path

from pymopsmap.classes import OptiProps
from pymopsmap.utils import get_logger
import xarray as xr

logger = get_logger(__name__)


def create_lut_for_smartg(
    op: OptiProps, specie: str, output_directory: Path
) -> xr.Dataset:
    """
    Creates a LUT that can be used as input in the AerOPAC class of SMART-G.
    This LUT is computed for a given set of wavelengths and relative humidity.

    Parameters
    ----------
    op: OptiProps
        Output optical properties from Mopsmap.
    specie: str
        Name of the file to create. Filename will be <filename>.nc,
        and will be read as `AerOPAC(filename)` in Smart-G.

    Returns
    -------
        xr.Dataset
            The LUT for the specie of interest
    """
    logger.debug(
        "Processing Mopsmap output LUT to create filename: %s", specie
    )

    # Extract wavelength, humidity, theta and mueller index from op
    wls = op.coord("wl")
    rhs = op.coord("rh")
    thetas = op.coord("theta")

    # Ensure ndim for each DataArray
    assert len(op.sel("phase").dims) == 4
    assert len(op.sel("kext").dims) == 2
    assert len(op.sel("ssa").dims) == 2

    # Shaping output Smart-G DataArrays
    def transform(var: str, output_shape: str) -> str:
        input_shape = " ".join(map(str, op.sel(var).dims))
        return f"{input_shape} -> {output_shape}"

    transform_phase = transform("phase", "rh wl mueller_idx theta")
    phase_darr = xr.DataArray(
        rearrange(op.sel("phase").values, transform_phase),
        coords={"theta": thetas, "wav": wls, "hum": rhs},
        dims=["hum", "wav", "stk", "theta"],
    )
    transform_kext = transform("kext", "rh wl")
    ext_darr = xr.DataArray(
        rearrange(op.sel("kext").values, transform_kext),
        coords={"hum": rhs, "wav": wls},
        dims=["hum", "wav"],
    )
    transform_ssa = transform("ssa", "rh wl")
    ssa_darr = xr.DataArray(
        rearrange(op.sel("ssa").values, transform_ssa),
        coords={"hum": rhs, "wav": wls},
        dims=["hum", "wav"],
    )

    dataset: xr.Dataset = xr.Dataset(
        {"phase": phase_darr, "ext": ext_darr, "ssa": ssa_darr},
        attrs={
            "name": specie,
            "H_mix_min": 0,
            "H_mix_max": 99,
            "H_stra_min": 0,
            "H_stra_max": 0,
            "H_free_min": 0,
            "H_free_max": 0,
            "Z_mix": 2.0,
            "Z_free": 0.0,
            "Z_stra": 0.0,
            "date": datetime.today().strftime("%Y-%m-%d"),
            "source": "Created using MOPSMAP v1.0.",
        },
    )
    logger.info("Mopsmap output processed. Output: %s", dataset)
    output_filename = output_directory / f"{specie}_sol.nc"
    dataset.to_netcdf(output_filename)
    logger.info("Mopsmap output store in file %s.", output_filename)
    return dataset
