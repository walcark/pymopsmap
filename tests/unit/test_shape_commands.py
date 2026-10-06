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


class TestNonAbsorbingFraction:
    """
    A fraction of the mode that does not absorb at all.

    Section 3.1 of Gasteiger and Wiegner (2018) introduces it, and MOPSMAP
    raises the imaginary index of the rest so the average is unchanged
    (init_wavelength_refr.f90, line 223).
    """

    def _mode(self, fraction: float):
        from pymopsmap import MicroParameters
        from pymopsmap.psd import LognormalPSD
        from pymopsmap.shapes import Sphere

        return MicroParameters(
            wavelength=[0.45],
            n_real=1.53,
            n_imag=0.0083,
            shape=Sphere(),
            psd=LognormalPSD(rm=0.1, sigma=1.6, n=1e6, rmin=0.001, rmax=5.0),
            nonabs_fraction=fraction,
        )

    def test_it_is_written_only_when_asked_for(self) -> None:
        from pymopsmap.engine.commands import microparams_command
        from pymopsmap.engine.workspace import Workspace

        with Workspace() as workspace:
            without = microparams_command(self._mode(0.0), workspace)
            with_it = microparams_command(self._mode(0.5), workspace)

        assert "nonabs_fraction" not in without
        assert "mode 1 refrac nonabs_fraction 0.5" in with_it

    def test_one_is_refused(self) -> None:
        import pytest
        from pydantic import ValidationError

        # Every particle non-absorbing would divide by zero in MOPSMAP.
        with pytest.raises(ValidationError):
            self._mode(1.0)
