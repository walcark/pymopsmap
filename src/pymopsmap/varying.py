"""
Numbers that may instead be a field of numbers.

A size distribution plays two roles. It describes a species, where a
parameter is allowed to vary, and it describes one point of a computation,
where every parameter is a number. The first role is what these types carry:
a scalar passes the constraint it declares, an array passes through, and the
constraint is applied again to every value the array holds once a point is
materialised.

``MicroParameters`` never sees one of these. It is built point by point, from
numbers, and keeps its own strict fields.
"""

from __future__ import annotations

from typing import Annotated, Any

import xarray as xr
from pydantic import AfterValidator, BeforeValidator


def _bounds(
    low: float | None = None,
    high: float | None = None,
    *,
    strict_low: bool = False,
    strict_high: bool = False,
) -> Any:
    """Build the validator for one constraint, skipped for an array."""

    def check(value: Any) -> Any:
        if isinstance(value, xr.DataArray):
            return value
        number = float(value)
        if low is not None and (number <= low if strict_low else number < low):
            raise ValueError(
                f"should be {'greater than' if strict_low else 'at least'}"
                f" {low}, got {number}"
            )
        if high is not None and (
            number >= high if strict_high else number > high
        ):
            raise ValueError(
                f"should be {'less than' if strict_high else 'at most'}"
                f" {high}, got {number}"
            )
        return number

    return AfterValidator(check)


def _keep_array(value: Any) -> Any:
    """Let an array through before pydantic tries to make a float of it."""
    return value


Varying = Annotated[
    float | xr.DataArray, BeforeValidator(_keep_array), _bounds()
]
"""A number with no constraint, or a field of them."""

Positive = Annotated[
    float | xr.DataArray,
    BeforeValidator(_keep_array),
    _bounds(low=0.0, strict_low=True),
]
"""Strictly positive, or a field of such numbers."""

NonNegative = Annotated[
    float | xr.DataArray, BeforeValidator(_keep_array), _bounds(low=0.0)
]
"""Zero or more, or a field of such numbers."""

AboveOne = Annotated[
    float | xr.DataArray,
    BeforeValidator(_keep_array),
    _bounds(low=1.0, strict_low=True),
]
"""Strictly above one, as a geometric standard deviation must be."""

AtLeastOne = Annotated[
    float | xr.DataArray, BeforeValidator(_keep_array), _bounds(low=1.0)
]
"""One or more, as an aspect ratio must be."""

UnitInterval = Annotated[
    float | xr.DataArray,
    BeforeValidator(_keep_array),
    _bounds(low=0.0, high=1.0),
]
"""Between zero and one inclusive."""


def aspect_ratio_range(low: float, high: float) -> Any:
    """The bounded aspect ratio a tabulated distribution is defined over."""
    return Annotated[
        float | xr.DataArray,
        BeforeValidator(_keep_array),
        _bounds(low=low, high=high),
    ]
