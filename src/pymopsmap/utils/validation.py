"""Validation of interpolation targets against source grids."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from pymopsmap.exceptions import DomainError


def check_within_grid(
    axis: str,
    values: ArrayLike,
    grid: ArrayLike,
    context: str | None = None,
) -> None:
    """
    Raise if any requested value falls outside the source grid.

    xarray's ``interp`` does not extrapolate: it returns NaN silently, and a
    NaN refractive index then resolves to an arbitrary optical dataset file.
    Checking the domain up front turns that into a named error.

    Parameters
    ----------
    axis : str
        Name of the interpolated coordinate, used in the error message.
    values : ArrayLike
        Requested values.
    grid : ArrayLike
        Source grid the values are interpolated on.
    context : str, optional
        What is being read, for instance the specie name.

    Raises
    ------
    DomainError
        If any value lies outside ``[grid.min(), grid.max()]``.
    """
    requested = np.atleast_1d(np.asarray(values, dtype=float))
    bounds = (float(np.min(grid)), float(np.max(grid)))

    outside = requested[(requested < bounds[0]) | (requested > bounds[1])]
    if outside.size:
        raise DomainError(axis, float(outside[0]), bounds, context)
