from .cams_to_microparams import (
    read_aerosol_microphysical_parameters,
    CamsAerosol,
    CamsVersion,
)
from .optiprops_to_smartg import create_lut_for_smartg

from pymopsmap.mopsmap import compute_optical_properties
from pymopsmap.classes import extend_optiprops
from pymopsmap.utils import SortedPosFloat64List

from numpy.typing import ArrayLike
from pathlib import Path
import xarray as xr


_all_ = [
    "read_aerosol_microphysical_parameters",
    "create_lut_for_smartg",
    "cams_to_smartg",
    "CamsVersion",
    "CamsAerosol",
]


def cams_to_smartg(
    aerosol: CamsAerosol,
    version: CamsVersion,
    wl_microns: SortedPosFloat64List,
    rh: ArrayLike,
    output_directory: Path,
) -> xr.Dataset:
    """
    Format an opticam properties dataset for Smart-G from CAMS
    microphysical parameters.
    """

    index, mps = read_aerosol_microphysical_parameters(
        aerosol=aerosol, version=version, wl_microns=wl_microns, rh=rh
    )

    ops = [compute_optical_properties(mp=mp) for mp in mps]
    op_tot = extend_optiprops(index, ops)

    ds = create_lut_for_smartg(
        op_tot, specie=aerosol.value, output_directory=output_directory
    )

    return ds
