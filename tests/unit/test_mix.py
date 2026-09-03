"""External mixtures of species."""

from __future__ import annotations

import numpy as np
import pytest

import pymopsmap as pm
from tests.conftest import integrated_result

WL = [0.44, 0.55, 0.67]


@pytest.fixture
def engine(monkeypatch):
    """A recorder returning a result that depends on the modes it is given."""
    calls: list[dict] = []

    def fake_run_point(modes, output_types, rh, quiet):
        calls.append({"modes": modes, "rh": rh})
        wl = np.asarray(modes[0].wavelength, dtype=float)
        # kext proportional to the total number, so scaling is observable.
        total = sum(mode.psd.n for mode in modes)
        ds = integrated_result(wl, kext=1e-6 * max(total, 1.0))
        ds.attrs["shape_types"] = sorted({m.shape.type for m in modes})
        return ds

    monkeypatch.setattr("pymopsmap.engine.run_point", fake_run_point)
    return calls


class TestConstruction:
    def test_accepts_catalogue_entries(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 3.2e9, pm.CAMS.SEA_SALT: 1.1e8})

        assert [s.name for s in mix.species] == ["sulphate", "sea_salt"]

    def test_accepts_an_already_loaded_specie(self, engine):
        specie = pm.load(pm.CAMS.DUST)

        assert pm.Mix({specie: 1.0}).species == [specie]

    def test_weights_are_inspectable(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 3.2e9})

        assert mix.weights == {"sulphate": 3.2e9}

    def test_an_empty_mixture_is_refused(self):
        with pytest.raises(ValueError, match="empty"):
            pm.Mix({})


class TestScaling:
    def test_a_single_specie_at_its_stored_amplitude_is_itself(self, engine):
        specie = pm.load(pm.CAMS.SULPHATE)
        alone = specie.compute(wl=WL, rh=50.0)

        mixed = pm.Mix({specie: specie.amplitude}).compute(wl=WL, rh=50.0)

        np.testing.assert_allclose(
            mixed["kext"].values, alone["kext"].values, rtol=1e-12
        )

    def test_doubling_the_weight_doubles_the_extinction(self, engine):
        specie = pm.load(pm.CAMS.SULPHATE)
        single = pm.Mix({specie: specie.amplitude}).compute(wl=WL, rh=50.0)

        double = pm.Mix({specie: 2 * specie.amplitude}).compute(wl=WL, rh=50.0)

        np.testing.assert_allclose(
            double["kext"].values, 2 * single["kext"].values, rtol=1e-12
        )

    def test_the_intensive_properties_are_left_alone(self, engine):
        specie = pm.load(pm.CAMS.SULPHATE)

        double = pm.Mix({specie: 2 * specie.amplitude}).compute(wl=WL, rh=50.0)

        np.testing.assert_allclose(double["ssa"].values, 0.9)

    def test_the_ratio_between_modes_is_preserved(self, engine):
        """Sea salt is bimodal; scaling must not flatten fine on coarse."""
        specie = pm.load(pm.CAMS.SEA_SALT)
        stored = [
            float(specie.tree[mode].to_dataset()["n"]) for mode in specie.modes
        ]

        pm.Mix({specie: 10 * specie.amplitude}).compute(wl=WL, rh=50.0)

        used = [mode.psd.n for mode in engine[-1]["modes"]]
        assert used[0] / used[1] == pytest.approx(stored[0] / stored[1])


class TestCombination:
    def test_extinction_is_the_sum_of_the_species(self, engine):
        a, b = pm.load(pm.CAMS.SULPHATE), pm.load(pm.CAMS.SEA_SALT)
        alone_a = pm.Mix({a: a.amplitude}).compute(wl=WL, rh=50.0)
        alone_b = pm.Mix({b: b.amplitude}).compute(wl=WL, rh=50.0)

        mixed = pm.Mix({a: a.amplitude, b: b.amplitude}).compute(
            wl=WL, rh=50.0
        )

        np.testing.assert_allclose(
            mixed["kext"].values,
            alone_a["kext"].values + alone_b["kext"].values,
            rtol=1e-12,
        )

    def test_each_specie_runs_on_its_own(self, engine):
        pm.Mix({pm.CAMS.SULPHATE: 1.0, pm.CAMS.SEA_SALT: 1.0}).compute(
            wl=WL, rh=50.0
        )

        assert len(engine) == 2

    def test_a_swept_humidity_is_preserved(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 1.0, pm.CAMS.SEA_SALT: 1.0})

        out = mix.compute(wl=WL, rh=[0.0, 50.0, 90.0])

        assert out["kext"].dims == ("rh", "wl")
        assert len(engine) == 6


class TestAmplitude:
    def test_it_sums_the_modes(self, engine):
        specie = pm.load(pm.CAMS.SEA_SALT)

        assert specie.amplitude == pytest.approx(73.0)  # 70 fine + 3 coarse

    def test_a_monomodal_specie_reports_its_single_mode(self, engine):
        assert pm.load(pm.CAMS.SULPHATE).amplitude == pytest.approx(1.0)
