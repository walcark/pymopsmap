"""The resolved concentrations travel as a variable of the mixture."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap as pm
from tests.conftest import integrated_result

WL = [0.44, 0.55]


@pytest.fixture
def engine(monkeypatch):
    def fake_run_point(modes, output_types, rh, quiet):
        area = sum(mode.psd.n * mode.psd.rm**2 for mode in modes)
        return integrated_result(modes[0].wavelength, kext=area * 1e-9)

    monkeypatch.setattr("pymopsmap.engine.run_point", fake_run_point)


class TestScalarWeights:
    def test_the_concentrations_are_a_variable(self, engine):
        op = pm.Mix({pm.CAMS.SULPHATE: 3.2e9}).compute(wl=WL, rh=50.0)

        assert op["concentration"].dims == ("specie",)

    def test_each_species_is_named_on_its_axis(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 3.2e9, pm.CAMS.SEA_SALT: 1.1e8})

        op = mix.compute(wl=WL, rh=50.0)

        assert list(op["specie"].values) == ["sulphate", "sea_salt"]

    def test_the_requested_value_is_what_is_reported(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 3.2e9, pm.CAMS.SEA_SALT: 1.1e8})

        op = mix.compute(wl=WL, rh=50.0)

        assert float(
            op["concentration"].sel(specie="sulphate")
        ) == pytest.approx(3.2e9)

    def test_it_carries_its_unit(self, engine):
        op = pm.Mix({pm.CAMS.SULPHATE: 3.2e9}).compute(wl=WL, rh=50.0)

        assert op["concentration"].attrs["units"] == "m-3"


class TestFieldWeights:
    """A composition that varies from pixel to pixel."""

    @pytest.fixture
    def fields(self):
        ramp = xr.DataArray(
            np.linspace(0.0, 1.0, 4)[None, :].repeat(3, 0), dims=["y", "x"]
        )
        return (1.0 - ramp) * 1e9, ramp * 1e9

    def test_the_concentrations_keep_the_image_dimensions(
        self, engine, fields
    ):
        sulphate, sea_salt = fields
        mix = pm.Mix({pm.CAMS.SULPHATE: sulphate, pm.CAMS.SEA_SALT: sea_salt})

        op = mix.compute(wl=WL, rh=50.0)

        assert set(op["concentration"].dims) == {"specie", "y", "x"}

    def test_a_pixel_reports_its_own_composition(self, engine, fields):
        sulphate, sea_salt = fields
        mix = pm.Mix({pm.CAMS.SULPHATE: sulphate, pm.CAMS.SEA_SALT: sea_salt})

        op = mix.compute(wl=WL, rh=50.0)

        west = op["concentration"].isel(y=0, x=0)
        assert float(west.sel(specie="sulphate")) == pytest.approx(1e9)
        assert float(west.sel(specie="sea_salt")) == pytest.approx(0.0)

    def test_the_result_can_be_written(self, engine, fields, tmp_path):
        """An attribute holding an array could not be serialised."""
        sulphate, sea_salt = fields
        mix = pm.Mix({pm.CAMS.SULPHATE: sulphate, pm.CAMS.SEA_SALT: sea_salt})

        op = mix.compute(wl=WL, rh=50.0)
        op.to_netcdf(tmp_path / "scene.nc")

        assert (tmp_path / "scene.nc").exists()


class TestSpeciesAlone:
    def test_a_species_result_carries_no_concentration(self, engine):
        """It is a property of the mixture, not of the species."""
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=50.0)

        assert "concentration" not in op
        assert "specie" not in op.dims
