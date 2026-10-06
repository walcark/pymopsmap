"""Annotated numeric list types with pydantic validation."""

import math
from typing import Annotated, TypeAlias

import numpy as np
from pydantic import AfterValidator, BeforeValidator


# --------------------------------------------------------------------------
# Type validators
# --------------------------------------------------------------------------
def coerce_as_float_list(
    value: float | list[float] | np.ndarray,
) -> list[float]:
    """
    Bring a scalar or array onto a plain list of floats, unrounded.

    Values are kept at full double precision. Rounding here would be absolute,
    while the only tolerance MOPSMAP applies is relative (1e-6, in
    interpolate_linear when a single grid point is loaded), so an absolute
    round is coarser than that tolerance for small imaginary indices rather
    than protective of it.
    """
    return np.atleast_1d(np.asarray(value)).tolist()


def assert_finite(value: list[float]) -> list[float]:
    if any(not math.isfinite(v) for v in value):
        raise ValueError("Input value should be finite (no NaN, no infinity).")
    return value


def assert_strictly_positive(value: list[float]) -> list[float]:
    if any(v <= 0 for v in value):
        raise ValueError("Input value should be strictly positive.")
    return value


def assert_sorted(value: list[float]) -> list[float]:
    value_sorted = sorted(value)
    if any(v != vs for v, vs in zip(value, value_sorted)):
        raise ValueError("Input list should be sorted.")
    return value


# --------------------------------------------------------------------------
# Types definition
# --------------------------------------------------------------------------
# assert_finite guards every downstream type: a NaN passes the positivity
# check (nan <= 0 is False) and would otherwise reach the dataset resolver,
# where it silently selects an arbitrary refractive index grid point.
Float64List: TypeAlias = Annotated[
    list[float],
    BeforeValidator(coerce_as_float_list),
    AfterValidator(assert_finite),
]

PosFloat64List: TypeAlias = Annotated[
    Float64List, AfterValidator(assert_strictly_positive)
]

SortedFloat64List: TypeAlias = Annotated[
    Float64List, AfterValidator(assert_sorted)
]

SortedPosFloat64List: TypeAlias = Annotated[
    SortedFloat64List, AfterValidator(assert_strictly_positive)
]
