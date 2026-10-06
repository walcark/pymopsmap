"""
Size-parameter limits of the optical dataset, read from its index.

MOPSMAP refuses a particle larger than the file covering its refractive index,
and that limit is not one number: index.nc carries max_sizepara over
(mimag, mreal, eps, np). For spheres it barely moves, but for the merged
spheroids it spans three orders of magnitude, so guessing one value for all of
them lets a request through that the engine then rejects.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import xarray as xr

if TYPE_CHECKING:
    from pymopsmap.shapes import Shape

# The shape codes MOPSMAP uses internally (get_nc_filename.f90).
SPHERE = -1
SPHEROID_MERGED = -7
IRREGULAR = {
    "A": -109,
    "B": -110,
    "C": -111,
    "D": -118,
    "E": -126,
    "F": -127,
}

# Published limits of Gasteiger and Wiegner (2018), Tables 1 and 2, used when
# no index is at hand. They are conservative, and blind to the refractive
# index.
PUBLISHED_MAXIMUM: dict[str, float] = {
    "sphere": 1005.0,
    "spheroid": 1005.0,
    "spheroid-lognormal": 1005.0,
    "spheroid-distr-file": 1005.0,
    "irregular": 30.2,
    "irregular-distr-file": 30.2,
    "irregular-overlay": 30.2,
}


class SizeParameterLimits:
    """The largest size parameter the dataset covers, per shape and index."""

    def __init__(self, index_path: str | Path | None):
        self._table: xr.DataArray | None = None
        if index_path is not None:
            with xr.open_dataset(index_path) as ds:
                if "max_sizepara" in ds:
                    self._table = ds["max_sizepara"].load()

    def maximum(self, shape: Shape, n_real: float, n_imag: float) -> float:
        """
        The limit for one shape at one refractive index.

        Parameters
        ----------
        shape : Shape
            The particle shape.
        n_real, n_imag : float
            The refractive index MOPSMAP will work on.

        Returns
        -------
        float
            The largest size parameter covered, falling back to the published
            table when the index is unavailable.
        """
        published = PUBLISHED_MAXIMUM.get(shape.type, float("inf"))
        if self._table is None:
            return published

        code = _shape_code(shape)
        if code is None or code not in self._table["np"].values:
            return published

        entry = self._table.sel(np=code).sel(
            mreal=n_real, mimag=n_imag, method="nearest"
        )
        if "eps" in entry.dims:
            ratio = 1.0 if code == SPHERE else _aspect_ratio(shape)
            if ratio is None:
                # The aspect ratios come from a distribution: the run is only
                # as covered as the least covered ratio it may reach.
                entry = entry.where(entry > 0.0).min()
            else:
                entry = entry.sel(eps=ratio, method="nearest")

        limit = float(entry)
        # eps = 1 sits on the grid with a limit of zero: a perfect sphere has
        # no spheroid file. Reading it as a limit rejects every size.
        return limit if limit > 0.0 else published


def _shape_code(shape: Shape) -> int | None:
    """Map a shape onto the code MOPSMAP files it under."""
    if shape.type == "sphere":
        return SPHERE
    if shape.type.startswith("spheroid"):
        # An aspect ratio of one is a sphere, and is read from the sphere
        # files (make_contributions.f90, i_np = 1).
        ratio = _aspect_ratio(shape)
        if ratio is not None and abs(ratio - 1.0) < 1e-3:
            return SPHERE
        return SPHEROID_MERGED
    if shape.type == "irregular":
        return IRREGULAR.get(shape.shape_id)  # type: ignore[union-attr]
    return None


def _aspect_ratio(shape: Shape) -> float | None:
    """
    The one aspect ratio of a shape, or None when a distribution sets it.

    A sphere is a spheroid of aspect ratio one. A shape reading its ratios
    from a file spans many, so no single row of the index applies.
    """
    if shape.type == "sphere":
        return 1.0
    ratio = getattr(shape, "aspect_ratio", None)
    return None if ratio is None else float(ratio)
