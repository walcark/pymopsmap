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
