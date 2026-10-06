"""How a nonspherical size is read, carried by the mode that is read."""

from __future__ import annotations

import pytest

from pymopsmap import MicroParameters, Specie
from pymopsmap.psd import LognormalPSD
from pymopsmap.shapes import Irregular, Sphere
from pymopsmap.species import Mode


def _mode(size_equ: str = "cs") -> MicroParameters:
    return MicroParameters(
        wavelength=[0.532],
        n_real=1.54,
        n_imag=0.005,
        shape=Irregular(shape_id="D"),
        psd=LognormalPSD(rm=0.5, sigma=2.0, n=1e8, rmin=0.001, rmax=1.75),
        size_equ=size_equ,  # type: ignore[arg-type]
    )


class TestOnAMode:
    def test_it_defaults_to_the_cross_section(self):
        """MOPSMAP's own default (read_input.f90, size_equ = 0)."""
        assert _mode().size_equ == "cs"

    @pytest.mark.parametrize("value", ["cs", "vol", "vol_cs_ratio"])
    def test_the_three_mopsmap_takes_are_accepted(self, value):
        assert _mode(value).size_equ == value

    def test_anything_else_is_refused(self):
        from pydantic import ValidationError

        with pytest.raises(ValidationError):
            _mode("radius")


class TestInTheLaunchFile:
    def test_it_is_written_from_the_modes(self, tmp_path):
        from pymopsmap.engine.launch_file import write_launching_file
        from pymopsmap.engine.workspace import Workspace

        with Workspace() as workspace:
            paths = write_launching_file([_mode("vol")], workspace)
            content = paths["mopsmap"].read_text()

        assert "size_equ vol" in content

    def test_modes_that_disagree_are_refused(self):
        """MOPSMAP reads one size_equ for the whole run, not one per mode."""
        from pymopsmap.engine.launch_file import write_launching_file
        from pymopsmap.engine.workspace import Workspace

        with Workspace() as workspace:
            with pytest.raises(ValueError, match="one size equivalence"):
                write_launching_file([_mode("cs"), _mode("vol")], workspace)


class TestOnASpecies:
    def test_it_survives_a_round_trip_through_the_catalogue(self, tmp_path):
        """
        Table 5 of the article cannot be expressed without it.

        Its four columns are one size distribution read three ways, and the
        mass-to-backscatter factor differs between them by a factor of two.
        """
        specie = Specie.custom(
            Mode(
                shape=Irregular(shape_id="D"),
                psd=LognormalPSD(
                    rm=0.5, sigma=2.0, n=1e8, rmin=0.001, rmax=1.75
                ),
                n_real=1.54,
                n_imag=0.005,
                size_equ="vol_cs_ratio",
            ),
            name="aggregate",
        )
        path = tmp_path / "aggregate.nc"
        specie.to_netcdf(path)

        from pymopsmap import load

        reloaded = load(path)

        assert reloaded.at(wl=[0.532]).modes[0].size_equ == "vol_cs_ratio"

    def test_a_sphere_keeps_the_default(self):
        specie = Specie.custom(
            Mode(
                shape=Sphere(),
                psd=LognormalPSD(
                    rm=0.5, sigma=2.0, n=1e8, rmin=0.001, rmax=1.75
                ),
                n_real=1.54,
                n_imag=0.005,
            ),
            name="drop",
        )

        assert specie.at(wl=[0.532]).modes[0].size_equ == "cs"
