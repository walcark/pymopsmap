"""
Normalisation of a parameter space into the points to compute.

Kept apart from the iteration and the assembly on purpose: those two are what a
sweep engine replaces, while this normalisation stays as it is and only grows
the input forms it accepts.
"""

from __future__ import annotations

from itertools import product
from typing import Any

import numpy as np
import xarray as xr


def as_space(**axes: Any) -> tuple[list[dict[str, Any]], list[str]]:
    """
    Turn keyword arguments into the points to compute.

    A scalar, or ``None``, fixes an axis without sweeping it. A list or an
    array sweeps it, and axes sweep independently, so distinct axes multiply.
    A single-element list still sweeps: passing a list expresses the intent to
    keep that dimension in the result.

    Parameters
    ----------
    **axes : Any
        One entry per parameter, scalar or sequence.

    Returns
    -------
    points : list of dict
        The cartesian product of the swept axes, each point carrying every
        axis by name.
    dims : list of str
        The axes that are swept, in the order they were given. They become the
        dimensions of the result.
    """
    values: dict[str, list[Any]] = {}
    dims: list[str] = []

    for name, value in axes.items():
        if _is_sequence(value):
            values[name] = [_scalar(v) for v in value]
            dims.append(name)
        else:
            values[name] = [value]

    points = [
        dict(zip(values, combination))
        for combination in product(*values.values())
    ]
    return points, dims


def _is_sequence(value: Any) -> bool:
    return isinstance(value, (list, tuple, np.ndarray))


def _scalar(value: Any) -> Any:
    return value.item() if isinstance(value, np.generic) else value


def assemble(
    index: list[dict[str, Any]], results: list[xr.Dataset]
) -> xr.Dataset:
    """
    Stack per-point results onto the dimensions that were swept.

    Parameters
    ----------
    index : list of dict
        One entry per result, holding the swept axes and their values.
    results : list of xr.Dataset
        The per-point results, in the same order.

    Returns
    -------
    xr.Dataset
        A single dataset carrying one dimension per swept axis.
    """
    if not results:
        raise ValueError("Cannot assemble an empty sweep.")
    if len(index) != len(results):
        raise ValueError(
            f"index and results must have the same length, got "
            f"{len(index)} and {len(results)}."
        )

    axes = list(index[0])
    ds = xr.concat(results, dim="run")
    for axis in axes:
        ds = ds.assign_coords(
            {axis: ("run", np.asarray([point[axis] for point in index]))}
        )

    if len(axes) == 1:
        return ds.swap_dims({"run": axes[0]})
    return ds.set_index(run=axes).unstack("run")
