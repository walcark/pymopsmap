"""The MOPSMAP command each shape writes."""

from __future__ import annotations

import pytest

from pymopsmap.shapes import (
    Irregular,
    IrregularDistrFile,
    IrregularOverlay,
    Sphere,
    Spheroid,
    SpheroidDistrFile,
    SpheroidLognormal,
)

# A path with a slash in it, which is what every real one has.
PATH = "/home/someone/mopsmap/data/ar_kandler"


@pytest.mark.parametrize(
    "shape",
    [
        SpheroidDistrFile(distr_filename=PATH),
        IrregularDistrFile(distr_filename=PATH),
        IrregularOverlay(distr_filename=PATH, xmin=1.0, xmax=30.0),
    ],
)
def test_a_file_path_is_quoted(shape) -> None:
    """
    A list-directed Fortran read ends the record on an unquoted slash.

    Without the quotes MOPSMAP leaves its filename variable untouched and
    reports a missing file whose name is uninitialised memory.
    """
    assert f"'{PATH}'" in shape.command


def test_the_overlay_writes_its_bounds() -> None:
    """The overlay command interpolates its size parameter range."""
    shape = IrregularOverlay(distr_filename=PATH, xmin=1.5, xmax=28.0)

    assert shape.command.endswith("1.5 28.0")


@pytest.mark.parametrize(
    ("shape", "expected"),
    [
        (Sphere(), "shape sphere"),
        (
            Spheroid(mode="prolate", aspect_ratio=2.0),
            "shape spheroid prolate 2.0",
        ),
        (Irregular(shape_id="D"), "shape irregular D"),
        (
            SpheroidLognormal(
                zeta1=0.5, zeta2=0.5, aspect_ratio=2.0, sigma_ar=0.6
            ),
            "shape spheroid log_normal 0.5 0.5 2.0 0.6",
        ),
    ],
)
def test_the_other_shapes_keep_their_command(shape, expected: str) -> None:
    assert shape.command == expected
