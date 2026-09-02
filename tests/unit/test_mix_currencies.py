"""Expressing mixture weights as masses or as optical depth fractions."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap as pm

WL = [0.44, 0.55, 0.67]


@pytest.fixture
def engine(monkeypatch):
    """
    A result that depends on the modes, so species respond to humidity
    differently and their optical depth fractions actually move.
    """

    def fake_run_point(modes, output_types, rh, quiet):
        wl = np.asarray(modes[0].wavelength, dtype=float)
        area = sum(mode.psd.n * mode.psd.rm**2 for mode in modes)
        volume = sum(mode.psd.n * mode.psd.rm**3 for mode in modes)
        kext = np.full(len(wl), area * 1e-6)
        ds = xr.Dataset(
            {
                "kext": (("wl",), kext),
                "ssa": (("wl",), np.full(len(wl), 0.9)),
                "ksca": (("wl",), kext * 0.9),
                "mass_conc": volume * 1e-3,
            },
            coords={"wl": wl},
        )
        ds.attrs["shape_types"] = sorted({m.shape.type for m in modes})
        return ds

    monkeypatch.setattr("pymopsmap.engine.run_point", fake_run_point)


def _fraction(op: xr.Dataset, part: xr.Dataset, wl: float) -> float:
    return float(
        part["kext"].sel(wl=wl) / op["kext"].sel(wl=wl)  # noqa: PD011
    )


class TestOpticalDepthFractions:
    def test_the_fractions_hold_at_the_reference_point(self, engine):
        mix = pm.Mix.from_optical_depth(
            {pm.CAMS.SULPHATE: 0.30, pm.CAMS.DUST: 0.70},
            wl_ref=0.55,
            rh_ref=50.0,
        )

        op = mix.compute(wl=WL, rh=50.0)
        shares = mix.contributions(op)

        assert shares["sulphate"] == pytest.approx(0.30, rel=1e-9)
        assert shares["dust"] == pytest.approx(0.70, rel=1e-9)

    def test_the_extinction_is_normalised_at_the_reference(self, engine):
        mix = pm.Mix.from_optical_depth(
            {pm.CAMS.SULPHATE: 0.30, pm.CAMS.DUST: 0.70},
            wl_ref=0.55,
            rh_ref=50.0,
        )

        op = mix.compute(wl=WL, rh=50.0)

        assert float(op["kext"].sel(wl=0.55)) == pytest.approx(1.0)

    def test_the_fractions_move_at_another_humidity(self, engine):
        """Species do not grow alike, so a fraction is not an invariant."""
        mix = pm.Mix.from_optical_depth(
            {pm.CAMS.SULPHATE: 0.30, pm.CAMS.SEA_SALT: 0.70},
            wl_ref=0.55,
            rh_ref=0.0,
        )

        op = mix.compute(wl=WL, rh=90.0)

        assert mix.contributions(op)["sulphate"] != pytest.approx(0.30)

    def test_fractions_are_normalised(self, engine):
        mix = pm.Mix.from_optical_depth(
            {pm.CAMS.SULPHATE: 3.0, pm.CAMS.DUST: 7.0},
            wl_ref=0.55,
            rh_ref=50.0,
        )

        op = mix.compute(wl=WL, rh=50.0)

        assert mix.contributions(op)["sulphate"] == pytest.approx(0.30)

    def test_the_reference_humidity_is_required(self):
        with pytest.raises(TypeError):
            pm.Mix.from_optical_depth({pm.CAMS.SULPHATE: 1.0}, wl_ref=0.55)

    def test_the_reference_wavelength_must_be_covered(self, engine):
        mix = pm.Mix.from_optical_depth(
            {pm.CAMS.SULPHATE: 1.0}, wl_ref=5.0, rh_ref=50.0
        )

        with pytest.raises(pm.DomainError, match="wl"):
            mix.compute(wl=WL, rh=50.0)


class TestMassConcentrations:
    def test_the_requested_mass_is_reached(self, engine):
        mix = pm.Mix.from_mass(
            {pm.CAMS.SULPHATE: 4.1e-9, pm.CAMS.DUST: 2.2e-8}, rh_ref=0.0
        )

        op = mix.compute(wl=WL, rh=0.0)

        assert float(op["mass_conc"]) == pytest.approx(4.1e-9 + 2.2e-8)

    def test_the_reference_humidity_is_required(self):
        with pytest.raises(TypeError):
            pm.Mix.from_mass({pm.CAMS.SULPHATE: 1e-9})


class TestResolvedConcentrations:
    def test_they_travel_with_the_result(self, engine):
        mix = pm.Mix.from_optical_depth(
            {pm.CAMS.SULPHATE: 0.30, pm.CAMS.DUST: 0.70},
            wl_ref=0.55,
            rh_ref=50.0,
        )

        op = mix.compute(wl=WL, rh=50.0)

        assert op.attrs["species"] == ["sulphate", "dust"]
        assert len(op.attrs["concentrations"]) == 2
        assert all(c > 0 for c in op.attrs["concentrations"])

    def test_a_number_mixture_reports_what_was_asked(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 3.2e9})

        op = mix.compute(wl=WL, rh=50.0)

        assert op.attrs["concentrations"] == pytest.approx([3.2e9])

    def test_weights_report_the_request_in_its_own_currency(self, engine):
        mix = pm.Mix.from_optical_depth(
            {pm.CAMS.SULPHATE: 0.30, pm.CAMS.DUST: 0.70},
            wl_ref=0.55,
            rh_ref=50.0,
        )

        assert mix.weights == {"sulphate": 0.30, "dust": 0.70}
        assert mix.currency.value == "optical_depth"
