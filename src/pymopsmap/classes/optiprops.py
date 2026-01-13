"""
optiprops.py

Author  : Kévin Walcarius
Date    : 2026-01-08
Version : 1.0
License : MIT
Summary : Definition of the OptiProps class, encapsulating an
          xr.Dataset with all the optical properties computed
          by Mopsmap.
"""

from dataclasses import dataclass
from typing import Union
import xarray as xr
import numpy as np


@dataclass(frozen=True)
class OptiProps:
    ds: xr.Dataset

    def save(self, path: str) -> None:
        self.ds.to_netcdf(path)

    def sel(self, prop: str, **kwargs) -> xr.DataArray:
        return self.ds[prop].sel(**kwargs)

    def coord(self, coord: str) -> np.ndarray:
        return self.ds[coord].values


def extend_optiprops(
    index: list[dict[str, float]], optiprops_li: list[OptiProps]
) -> OptiProps:
    """
    Concatenate multiple OptiProps given an `index` of parameters
    used to generate each partial OptiProps.

    """
    if len(index) != len(optiprops_li):
        raise ValueError(
            "index and optiprops_li must have same length, "
            f"got {len(index)} vs {len(optiprops_li)}",
        )
    if not optiprops_li:
        raise ValueError("optiprops_li is empty")

    ds = xr.concat([op.ds for op in optiprops_li], dim="run")

    # attach index as coords on run (supports float/int/str)
    keys = list(index[0].keys())
    for k in keys:
        vals = [idx[k] for idx in index]
        ds = ds.assign_coords({k: ("run", np.asarray(vals))})

    if len(keys) == 1:
        ds = ds.swap_dims({"run": k})
    else:
        ds = ds.set_index(run=keys).unstack("run")

    return OptiProps(ds=ds)
