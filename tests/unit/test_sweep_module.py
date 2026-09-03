"""The sweep engine: one MOPSMAP run per point, memoised in a store."""

from __future__ import annotations

import pytest

import pymopsmap as pm
from tests.conftest import integrated_result

WL = [0.44, 0.55, 0.67]


@pytest.fixture
def engine(monkeypatch):
    """Count the MOPSMAP runs a sweep actually performs."""
    calls: list[dict] = []

    def fake_run_point(modes, output_types, rh, quiet):
        calls.append({"rh": rh, "rm": modes[0].psd.rm})
        # A tabulated species is already wet when it reaches the engine, so
        # what varies between points is its size distribution, not `rh`.
        return integrated_result(
            modes[0].wavelength, kext=1e-6 * modes[0].psd.rm
        )

    monkeypatch.setattr("pymopsmap.engine.run_point", fake_run_point)
    return calls


class TestShape:
    def test_a_scalar_humidity_gives_no_extra_dimension(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=50.0)

        assert op["kext"].dims == ("wl",)

    def test_a_list_adds_its_dimension(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[0.0, 50.0, 90.0])

        assert op["kext"].dims == ("rh", "wl")
        assert list(op["rh"].values) == [0.0, 50.0, 90.0]

    def test_one_run_per_point(self, engine):
        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[0.0, 50.0, 90.0])

        assert len(engine) == 3

    def test_the_wavelengths_are_computed_in_one_run(self, engine):
        """MOPSMAP takes the whole spectral grid, so wl never loops."""
        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=50.0)

        assert len(engine) == 1


class TestMemoisation:
    def test_recomputing_the_same_grid_runs_nothing(self, engine):
        specie = pm.load(pm.CAMS.SULPHATE)
        specie.compute(wl=WL, rh=[0.0, 50.0])
        engine.clear()

        specie.compute(wl=WL, rh=[0.0, 50.0])

        assert engine == []

    def test_a_wider_grid_is_a_different_sweep(self, engine):
        """
        A store holds one grid: its axes are fixed when it is created. Asking
        for a wider one is a new sweep, not an extension of the first.
        """
        specie = pm.load(pm.CAMS.SULPHATE)
        specie.compute(wl=WL, rh=[0.0, 50.0])
        engine.clear()

        specie.compute(wl=WL, rh=[0.0, 50.0, 90.0])

        assert len(engine) == 3

    def test_two_species_do_not_share_their_entries(self, engine):
        """The store key carries the species, not only the swept point."""
        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=50.0)
        engine.clear()

        pm.load(pm.CAMS.SEA_SALT).compute(wl=WL, rh=50.0)

        assert len(engine) == 1

    def test_a_different_wavelength_grid_is_a_different_sweep(self, engine):
        specie = pm.load(pm.CAMS.SULPHATE)
        specie.compute(wl=WL, rh=50.0)
        engine.clear()

        specie.compute(wl=[0.35, 0.55], rh=50.0)

        assert len(engine) == 1


class TestValues:
    def test_each_point_keeps_its_own_result(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[0.0, 90.0])

        dry = float(op["kext"].sel(rh=0.0).isel(wl=0))
        wet = float(op["kext"].sel(rh=90.0).isel(wl=0))
        assert wet > dry

    def test_a_mixture_still_combines(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 1.0, pm.CAMS.SEA_SALT: 1.0})

        op = mix.compute(wl=WL, rh=[0.0, 50.0])

        assert op["kext"].dims == ("rh", "wl")
