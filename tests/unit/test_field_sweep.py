"""A swept parameter can be a field, one value per pixel of an image."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap as pm
from tests.conftest import integrated_result

WL = [0.44, 0.55]

# Humidity over a small image, with values that repeat between pixels. They
# straddle the deliquescence step of CAMS sulphate, whose modal radius does not
# move below 30 percent, so each distinct value gives a distinct result.
RH_FIELD = xr.DataArray(
    [[10.0, 50.0, 10.0], [90.0, 50.0, 10.0]], dims=["y", "x"]
)


@pytest.fixture
def engine(monkeypatch):
    calls: list[float] = []

    def fake_run_point(modes, output_types, rh, quiet):
        calls.append(modes[0].psd.rm)
        return integrated_result(
            modes[0].wavelength, kext=1e-6 * modes[0].psd.rm
        )

    monkeypatch.setattr("pymopsmap.engine.run_point", fake_run_point)
    return calls


class TestFieldShape:
    def test_the_result_keeps_the_image_dimensions(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=RH_FIELD)

        assert op["kext"].dims == ("y", "x", "wl")
        assert op.sizes["y"] == 2 and op.sizes["x"] == 3

    def test_each_pixel_carries_its_own_value(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=RH_FIELD)

        driest = float(op["kext"].isel(y=0, x=0, wl=0))
        wettest = float(op["kext"].isel(y=1, x=0, wl=0))
        assert wettest > driest

    def test_pixels_sharing_a_value_share_a_result(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=RH_FIELD)

        assert float(op["kext"].isel(y=0, x=0, wl=0)) == pytest.approx(
            float(op["kext"].isel(y=1, x=2, wl=0))
        )


class TestDeduplication:
    def test_a_repeated_value_costs_one_run(self, engine):
        """Six pixels, three distinct humidities."""
        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=RH_FIELD)

        assert len(engine) == 3

    def test_a_uniform_field_costs_a_single_run(self, engine):
        uniform = xr.DataArray(np.full((4, 4), 50.0), dims=["y", "x"])

        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=uniform)

        assert len(engine) == 1


class TestNaming:
    def test_the_dimension_names_are_the_caller_s(self, engine):
        named = xr.DataArray([50.0, 70.0, 90.0], dims="rh_nominal")

        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=named)

        assert op["kext"].dims == ("rh_nominal", "wl")

    def test_a_plain_list_still_names_the_dimension_after_the_parameter(
        self, engine
    ):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[0.0, 50.0])

        assert op["kext"].dims == ("rh", "wl")


class TestMixtures:
    def test_a_mixture_accepts_a_field_too(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 1.0, pm.CAMS.SEA_SALT: 1.0})

        op = mix.compute(wl=WL, rh=RH_FIELD)

        assert op["kext"].dims == ("y", "x", "wl")
