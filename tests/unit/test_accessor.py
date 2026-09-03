"""The .mopsmap accessor on a result dataset."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap  # noqa: F401  registers the accessor

RH = [0.0, 50.0, 90.0]
WL = [0.44, 0.55]
THETA = np.linspace(0.0, 180.0, 7)


@pytest.fixture
def result() -> xr.Dataset:
    shape = (len(RH), len(WL))
    return xr.Dataset(
        {
            "kext": (("rh", "wl"), np.full(shape, 1e-4)),
            "ssa": (("rh", "wl"), np.full(shape, 0.95)),
            "phase": (
                ("rh", "wl", "theta"),
                np.ones((*shape, len(THETA))),
            ),
        },
        coords={"rh": RH, "wl": WL, "theta": THETA},
    )


class TestSmartgExport:
    def test_it_writes_a_lut(self, result, tmp_path):
        path = tmp_path / "sulphate_sol.nc"

        result.mopsmap.to_smartg(path, name="sulphate")

        assert path.exists()

    def test_the_lut_carries_the_smartg_dimension_names(
        self, result, tmp_path
    ):
        path = tmp_path / "lut.nc"
        result.mopsmap.to_smartg(path, name="sulphate")

        lut = xr.open_dataset(path)
        assert lut["ext"].dims == ("hum", "wav")
        assert lut["phase"].dims == ("hum", "wav", "stk", "theta")

    def test_the_humidity_dimension_can_be_named(self, result, tmp_path):
        """v1 hardcoded 'rh', which broke as soon as a user named it."""
        renamed = result.rename(rh="rh_nominal")
        path = tmp_path / "lut.nc"

        renamed.mopsmap.to_smartg(
            path, name="sulphate", humidity_dim="rh_nominal"
        )

        assert xr.open_dataset(path).sizes["hum"] == len(RH)

    def test_a_missing_humidity_dimension_is_reported(self, result, tmp_path):
        with pytest.raises(KeyError, match="humidity"):
            result.mopsmap.to_smartg(
                tmp_path / "lut.nc", name="sulphate", humidity_dim="absent"
            )

    def test_a_result_without_a_phase_function_is_reported(
        self, result, tmp_path
    ):
        with pytest.raises(KeyError, match="phase"):
            result.drop_vars("phase").mopsmap.to_smartg(
                tmp_path / "lut.nc", name="sulphate"
            )

    def test_the_values_survive_the_transposition(self, result, tmp_path):
        path = tmp_path / "lut.nc"
        marked = result.copy()
        marked["kext"][2, 1] = 7.0

        marked.mopsmap.to_smartg(path, name="sulphate")

        lut = xr.open_dataset(path)
        assert float(lut["ext"].isel(hum=2, wav=1)) == pytest.approx(7.0)
