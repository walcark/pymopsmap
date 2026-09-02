"""Assembling per-point results onto the swept dimensions."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from pymopsmap.sweep import assemble

WL = np.array([0.44, 0.55], dtype=np.float32)


def _result(value: float) -> xr.Dataset:
    return xr.Dataset(
        {"kext": (("wl",), np.full(2, value, dtype=np.float32))},
        coords={"wl": WL},
    )


class TestAssemble:
    def test_one_axis_becomes_a_dimension(self):
        out = assemble(
            [{"rh": 0.0}, {"rh": 50.0}], [_result(1.0), _result(2.0)]
        )

        assert out["kext"].dims == ("rh", "wl")
        assert list(out["rh"].values) == [0.0, 50.0]
        assert out["kext"].sel(rh=50.0).values.tolist() == [2.0, 2.0]

    def test_two_axes_unstack_into_a_grid(self):
        index = [
            {"rh": rh, "kappa": k}
            for rh in (0.0, 50.0)
            for k in (0.2, 0.4, 0.6)
        ]
        out = assemble(index, [_result(float(i)) for i in range(6)])

        assert set(out["kext"].dims) == {"rh", "kappa", "wl"}
        assert out.sizes["rh"] == 2
        assert out.sizes["kappa"] == 3

    def test_extra_dimensions_survive(self):
        phase = xr.Dataset(
            {"phase": (("wl", "theta"), np.ones((2, 5), dtype=np.float32))},
            coords={"wl": WL},
        )
        out = assemble([{"rh": 0.0}, {"rh": 50.0}], [phase, phase])

        assert out["phase"].dims == ("rh", "wl", "theta")

    def test_a_mismatched_index_is_rejected(self):
        with pytest.raises(ValueError, match="same length"):
            assemble([{"rh": 0.0}], [_result(1.0), _result(2.0)])

    def test_an_empty_sweep_is_rejected(self):
        with pytest.raises(ValueError, match="empty"):
            assemble([], [])
