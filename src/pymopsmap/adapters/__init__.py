"""Adapters — input converters (CAMS) and output converters (SMART-G)."""

from __future__ import annotations

from pathlib import Path

import xarray as xr
from numpy.typing import ArrayLike

from pymopsmap.models import OptiProps
from pymopsmap.models.output_request import DEFAULT_OUTPUT, OutputRequest
from pymopsmap.utils import SortedPosFloat64List

from .input.cams import (
    CamsAerosol,
    CamsVersion,
    read_aerosol_microphysical_parameters,
)
from .output.smartg import create_lut_for_smartg


def cams_to_smartg(
    aerosol: CamsAerosol,
    version: CamsVersion,
    wl_microns: SortedPosFloat64List,
    rh: ArrayLike,
    output_directory: Path,
) -> xr.Dataset:
    return create_lut_for_smartg(
        op=cams_to_optiprops(aerosol, version, wl_microns, rh),
        specie=aerosol.value,
        output_directory=output_directory,
    )


def cams_to_kext(
    aerosol: CamsAerosol,
    version: CamsVersion,
    wl_microns: SortedPosFloat64List,
    rh: ArrayLike,
    output_types: OutputRequest = DEFAULT_OUTPUT,
    quiet: bool = False,
) -> xr.DataArray:
    op: OptiProps = cams_to_optiprops(
        aerosol, version, wl_microns, rh, output_types=output_types, quiet=quiet
    )
    return op.sel("kext")


def cams_to_optiprops(
    aerosol: CamsAerosol,
    version: CamsVersion,
    wl_microns: SortedPosFloat64List,
    rh: ArrayLike,
    output_types: OutputRequest = DEFAULT_OUTPUT,
    quiet: bool = False,
) -> OptiProps:
    from pymopsmap.engine import compute_optical_properties

    dispatch = read_aerosol_microphysical_parameters(
        aerosol=aerosol, version=version, wl_microns=wl_microns, rh=rh
    )
    return compute_optical_properties(dispatch, output_types=output_types, quiet=quiet)


__all__ = [
    "CamsAerosol",
    "CamsVersion",
    "read_aerosol_microphysical_parameters",
    "create_lut_for_smartg",
    "cams_to_smartg",
    "cams_to_kext",
    "cams_to_optiprops",
]
