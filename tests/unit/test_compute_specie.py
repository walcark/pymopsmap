"""Computing a species over a scalar or a list of humidities."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap as pm
from pymopsmap.sweep import as_space

WL = [0.44, 0.55, 0.67]


@pytest.fixture
def engine(monkeypatch):
    """Replace the MOPSMAP run by a recorder returning a synthetic result."""
    calls: list[dict] = []

    def fake_run_point(modes, output_types, rh, quiet):
        calls.append({"modes": modes, "rh": rh, "outputs": output_types})
        wl = np.asarray(modes[0].wavelength, dtype=np.float32)
        return xr.Dataset(
            {"kext": (("wl",), np.arange(len(wl), dtype=np.float32))},
            coords={"wl": wl},
        )

    monkeypatch.setattr("pymopsmap.engine.run_point", fake_run_point)
    return calls


class TestAsSpace:
    def test_a_scalar_axis_is_not_swept(self):
        points, dims = as_space(rh=50.0)

        assert points == [{"rh": 50.0}]
        assert dims == []

    def test_an_absent_axis_is_not_swept(self):
        points, dims = as_space(rh=None)

        assert points == [{"rh": None}]
        assert dims == []

    def test_a_list_axis_becomes_a_dimension(self):
        points, dims = as_space(rh=[0.0, 50.0, 90.0])

        assert points == [{"rh": 0.0}, {"rh": 50.0}, {"rh": 90.0}]
        assert dims == ["rh"]

    def test_an_array_axis_becomes_a_dimension(self):
        points, dims = as_space(rh=np.array([0.0, 50.0]))

        assert len(points) == 2
        assert dims == ["rh"]

    def test_axes_multiply(self):
        points, dims = as_space(rh=[0.0, 50.0], kappa=[0.2, 0.4, 0.6])

        assert len(points) == 6
        assert dims == ["rh", "kappa"]

    def test_a_single_element_list_still_sweeps(self):
        """A list expresses intent: the axis is kept even with one value."""
        points, dims = as_space(rh=[50.0])

        assert points == [{"rh": 50.0}]
        assert dims == ["rh"]


class TestComputeScalar:
    def test_returns_a_result_indexed_by_wavelength(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=50.0)

        assert op["kext"].dims == ("wl",)
        assert len(engine) == 1

    def test_both_modes_reach_the_engine(self, engine):
        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=50.0)

        assert len(engine[0]["modes"]) == 2

    def test_a_tabulated_specie_is_not_grown_by_the_engine(self, engine):
        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=50.0)

        assert engine[0]["rh"] is None


class TestComputeList:
    def test_a_list_adds_a_humidity_dimension(self, engine):
        op = pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[0.0, 50.0, 90.0])

        assert op["kext"].dims == ("rh", "wl")
        assert list(op["rh"].values) == [0.0, 50.0, 90.0]

    def test_one_engine_run_per_humidity(self, engine):
        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[0.0, 50.0, 90.0])

        assert len(engine) == 3

    def test_each_run_carries_its_own_micro_parameters(self, engine):
        pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[0.0, 90.0])

        dry, wet = (call["modes"][0].psd.rm for call in engine)
        assert wet > dry


class TestValidationStillApplies:
    def test_an_out_of_range_humidity_raises_before_any_run(self, engine):
        from pymopsmap.exceptions import DomainError

        with pytest.raises(DomainError):
            pm.load(pm.CAMS.SULPHATE).compute(wl=WL, rh=[50.0, 99.0])

        assert engine == []

    def test_a_tabulated_specie_still_requires_a_humidity(self, engine):
        with pytest.raises(ValueError, match="rh"):
            pm.load(pm.CAMS.SULPHATE).compute(wl=WL)
