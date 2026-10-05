"""Building a species by hand, and saving it in the catalogue format."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap as pm

WL = [0.44, 0.55, 0.67]


def _mode(**overrides) -> pm.Mode:
    defaults = {
        "shape": pm.shapes.Sphere(),
        "psd": pm.psd.LognormalPSD(
            rm=0.05, sigma=1.5, n=1e9, rmin=0.001, rmax=40.0
        ),
        "n_real": 1.45,
        "n_imag": 1e-3,
        "density_dry": 1.8,
    }
    return pm.Mode(**(defaults | overrides))


class TestScalarSpecie:
    def test_it_computes_like_a_catalogue_one(self, monkeypatch):
        aer = pm.Specie.custom(_mode())

        point = aer.at(wl=WL)

        assert len(point.modes) == 1
        assert point.modes[0].psd.rm == pytest.approx(0.05)

    def test_it_is_dry_by_default(self):
        aer = pm.Specie.custom(_mode())

        assert aer.growth is pm.species.schema.Growth.NONE
        with pytest.raises(ValueError, match="does not model humidity"):
            aer.at(wl=WL, rh=50.0)

    def test_a_kappa_makes_it_hygroscopic(self):
        aer = pm.Specie.custom(_mode(kappa=0.3))

        assert aer.growth is pm.species.schema.Growth.KAPPA
        assert aer.at(wl=WL, rh=80.0).engine_rh == 80.0

    def test_several_modes_keep_their_order(self):
        fine = _mode(
            psd=pm.psd.LognormalPSD(
                rm=0.05, sigma=1.5, n=7e10, rmin=0.001, rmax=40.0
            )
        )
        coarse = _mode(
            psd=pm.psd.LognormalPSD(
                rm=1.0, sigma=2.0, n=3e9, rmin=0.001, rmax=40.0
            )
        )

        aer = pm.Specie.custom([fine, coarse])

        assert [m.psd.rm for m in aer.at(wl=WL).modes] == [0.05, 1.0]

    def test_the_amplitude_sums_the_modes(self):
        fine = _mode(
            psd=pm.psd.LognormalPSD(
                rm=0.05, sigma=1.5, n=7e10, rmin=0.001, rmax=40.0
            )
        )
        coarse = _mode(
            psd=pm.psd.LognormalPSD(
                rm=1.0, sigma=2.0, n=3e9, rmin=0.001, rmax=40.0
            )
        )

        assert pm.Specie.custom([fine, coarse]).amplitude == pytest.approx(
            7.3e10
        )


class TestSweptParameters:
    def test_a_data_array_becomes_a_swept_axis(self):
        aer = pm.Specie.custom(
            _mode(n_imag=xr.DataArray([1e-4, 1e-3, 1e-2], dims="absorption"))
        )

        assert aer.swept == {"absorption": 3}

    def test_a_swept_parameter_materialises_point_by_point(self):
        aer = pm.Specie.custom(
            _mode(n_imag=xr.DataArray([1e-4, 1e-2], dims="absorption"))
        )

        low = aer.at(wl=[0.55], absorption=0)
        high = aer.at(wl=[0.55], absorption=1)

        assert low.modes[0].n_imag[0] < high.modes[0].n_imag[0]

    def test_distinct_dimensions_multiply(self):
        """A swept parameter goes where its scalar would, in the model."""
        aer = pm.Specie.custom(
            _mode(
                psd=pm.psd.LognormalPSD(
                    rm=xr.DataArray(np.linspace(0.05, 0.5, 4), dims="rm"),
                    sigma=1.5,
                    n=1e9,
                    rmin=0.001,
                    rmax=40.0,
                ),
                n_imag=xr.DataArray([1e-4, 1e-2], dims="absorption"),
            )
        )

        assert aer.swept == {"rm": 4, "absorption": 2}

    def test_a_swept_psd_field_reaches_the_materialised_point(self):
        aer = pm.Specie.custom(
            _mode(
                psd=pm.psd.LognormalPSD(
                    rm=xr.DataArray([0.05, 0.5], dims="rm"),
                    sigma=1.5,
                    n=1e9,
                    rmin=0.001,
                    rmax=40.0,
                )
            )
        )

        assert aer.at(wl=WL, rm=1).modes[0].psd.rm == pytest.approx(0.5)

    def test_a_swept_shape_field_reaches_it_too(self):
        """Figure 2 of the article compares aspect ratios 1.4 and 3.0."""
        aer = pm.Specie.custom(
            _mode(
                shape=pm.shapes.Spheroid(
                    mode="prolate",
                    aspect_ratio=xr.DataArray([1.4, 3.0], dims="eps"),
                )
            )
        )

        assert aer.swept == {"eps": 2}
        assert aer.at(wl=WL, eps=1).modes[0].shape.aspect_ratio == 3.0


class TestRoundTrip:
    def test_it_saves_and_reloads(self, tmp_path):
        aer = pm.Specie.custom(_mode(kappa=0.3), name="my_aerosol")
        path = tmp_path / "my_aerosol.nc"

        aer.to_netcdf(path)
        back = pm.load(path)

        assert back.name == "my_aerosol"
        assert back.growth is pm.species.schema.Growth.KAPPA

    def test_the_reloaded_species_gives_the_same_point(self, tmp_path):
        aer = pm.Specie.custom(_mode(kappa=0.3))
        path = tmp_path / "aer.nc"
        aer.to_netcdf(path)

        original = aer.at(wl=WL, rh=50.0).modes[0]
        reloaded = pm.load(path).at(wl=WL, rh=50.0).modes[0]

        assert reloaded.psd.rm == pytest.approx(original.psd.rm)
        np.testing.assert_allclose(reloaded.n_real, original.n_real)
        assert reloaded.kappa == pytest.approx(original.kappa)

    def test_a_multi_mode_species_round_trips(self, tmp_path):
        aer = pm.Specie.custom([_mode(), _mode()])
        path = tmp_path / "aer.nc"

        aer.to_netcdf(path)

        assert len(pm.load(path).modes) == 2

    def test_a_catalogue_species_round_trips_too(self, tmp_path):
        path = tmp_path / "sulphate.nc"
        pm.load(pm.CAMS.SULPHATE).to_netcdf(path)

        back = pm.load(path)

        assert back.modes == ["fine", "coarse"]
        assert back.rh_range == (0.0, 95.0)
