"""A parameter declared on the species is swept by compute."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap as pm
from tests.conftest import integrated_result

WL = [0.44, 0.55]


@pytest.fixture
def engine(monkeypatch):
    """Record what each run was given, so the points can be told apart."""
    calls: list[dict] = []

    def fake_run_point(modes, output_types, rh, quiet, n_angles=2000):
        psd = modes[0].psd
        calls.append({"radius": getattr(psd, "radius", None), "rh": rh})
        return integrated_result(modes[0].wavelength, kext=1e-6 * psd.radius)

    monkeypatch.setattr("pymopsmap.engine.run_point", fake_run_point)
    return calls


def _single(radius) -> pm.Specie:
    return pm.Specie.custom(
        pm.Mode(
            shape=pm.shapes.Sphere(),
            psd=pm.psd.FixedPSD(radius=radius, n=1.0),
            n_real=1.56,
            n_imag=0.00215,
        )
    )


class TestDeclaredAxis:
    def test_it_becomes_a_dimension_of_the_result(self, engine):
        radii = xr.DataArray([0.1, 0.5, 2.0], dims="radius")

        op = _single(radii).compute(wl=WL)

        assert op["kext"].dims == ("radius", "wl")
        assert op.sizes["radius"] == 3

    def test_one_run_per_value(self, engine):
        _single(xr.DataArray([0.1, 0.5, 2.0], dims="radius")).compute(wl=WL)

        assert [call["radius"] for call in engine] == [0.1, 0.5, 2.0]

    def test_each_point_keeps_its_own_result(self, engine):
        op = _single(xr.DataArray([0.1, 2.0], dims="radius")).compute(wl=WL)

        small = float(op["kext"].isel(radius=0, wl=0))
        large = float(op["kext"].isel(radius=1, wl=0))
        assert large > small

    def test_a_scalar_parameter_adds_no_dimension(self, engine):
        op = _single(0.5).compute(wl=WL)

        assert op["kext"].dims == ("wl",)

    def test_the_caller_names_the_dimension(self, engine):
        radii = xr.DataArray([0.1, 0.5], dims="size_parameter")

        op = _single(radii).compute(wl=WL)

        assert op["kext"].dims == ("size_parameter", "wl")


class TestSeveralAxes:
    def test_distinct_dimensions_multiply(self, engine):
        aer = pm.Specie.custom(
            pm.Mode(
                shape=pm.shapes.Sphere(),
                psd=pm.psd.FixedPSD(
                    radius=xr.DataArray([0.1, 0.5, 2.0], dims="radius"),
                    n=1.0,
                ),
                n_real=1.56,
                n_imag=xr.DataArray([1e-4, 1e-2], dims="absorption"),
            )
        )

        op = aer.compute(wl=WL)

        assert set(op["kext"].dims) == {"radius", "absorption", "wl"}
        assert len(engine) == 6

    def test_a_shared_dimension_walks_a_trajectory(self, engine):
        """Two parameters on one dimension vary together, not as a grid."""
        aging = np.linspace(0.05, 0.5, 4)
        aer = pm.Specie.custom(
            pm.Mode(
                shape=pm.shapes.Sphere(),
                psd=pm.psd.LognormalPSD(
                    rm=xr.DataArray(aging, dims="aging"),
                    sigma=xr.DataArray(np.linspace(1.4, 2.1, 4), dims="aging"),
                    n=1e9,
                    rmin=0.001,
                    rmax=40.0,
                ),
                n_real=1.45,
                n_imag=1e-3,
            )
        )

        assert aer.swept == {"aging": 4}


class TestWithHumidity:
    def test_a_declared_axis_and_a_humidity_multiply(self, engine):
        aer = pm.Specie.custom(
            pm.Mode(
                shape=pm.shapes.Sphere(),
                psd=pm.psd.FixedPSD(
                    radius=xr.DataArray([0.1, 0.5], dims="radius"), n=1.0
                ),
                n_real=1.45,
                n_imag=1e-3,
                kappa=0.3,
            )
        )

        op = aer.compute(wl=WL, rh=[50.0, 90.0])

        assert set(op["kext"].dims) == {"radius", "rh", "wl"}
        assert len(engine) == 4
