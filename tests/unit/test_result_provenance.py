"""What produced a result travels with it, through the sweep store."""

from __future__ import annotations

import numpy as np
import pytest

import pymopsmap as pm
from tests.conftest import integrated_result

WL = [0.44, 0.55]


@pytest.fixture
def engine(monkeypatch):
    def fake_run_point(modes, output_types, rh, quiet, n_angles=2000):
        return integrated_result(modes[0].wavelength)

    monkeypatch.setattr("pymopsmap.engine.run_point", fake_run_point)


class TestShapeProvenance:
    def test_it_survives_the_sweep(self, engine):
        """The store keeps the numbers, not the attributes around them."""
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=50.0)

        assert op.attrs["shape_types"] == ["sphere"]

    def test_it_survives_a_swept_axis(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[0.0, 50.0])

        assert op.attrs["shape_types"] == ["sphere"]

    def test_the_effective_radius_is_not_lost_in_a_mixture(self, engine):
        """Without the provenance the combination cannot rebuild it."""
        mix = pm.Mix({pm.CAMS.SULPHATE: 1.0, pm.CAMS.SEA_SALT: 1.0})

        op = mix.compute(wl=WL, rh=50.0)

        assert np.isfinite(op["reff"]).all()

    def test_a_mixture_records_every_shape_it_came_from(self, engine):
        mix = pm.Mix({pm.CAMS.SULPHATE: 1.0, pm.CAMS.SEA_SALT: 1.0})

        assert mix.compute(wl=WL, rh=50.0).attrs["shape_types"] == ["sphere"]
