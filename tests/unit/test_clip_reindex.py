"""Clipped wavelengths are restored on the coordinate, not by position."""

from __future__ import annotations

import numpy as np
import xarray as xr

from pymopsmap.engine.coverage import reindex_to_full_grid

FULL_WL = [0.4, 0.5, 0.6, 0.7]
MASK = np.array([False, True, True, False])


def _kept(**extra_vars) -> xr.Dataset:
    """A result over the two kept wavelengths only."""
    return xr.Dataset(
        data_vars={"kext": (("wl",), np.array([1.0, 2.0], dtype=np.float32))}
        | extra_vars,
        coords={"wl": np.array([0.5, 0.6], dtype=np.float32)},
    )


class TestReindexToFullGrid:
    def test_scalar_variable_is_padded_with_nan(self):
        out = reindex_to_full_grid(_kept(), FULL_WL, MASK)

        assert list(out["wl"].values) == [np.float32(w) for w in FULL_WL]
        np.testing.assert_array_equal(
            np.isnan(out["kext"].values), [True, False, False, True]
        )
        assert out["kext"].sel(wl=0.5, method="nearest") == 1.0

    def test_extra_dimensions_are_preserved(self):
        phase = (("wl", "theta"), np.ones((2, 3), dtype=np.float32))
        out = reindex_to_full_grid(_kept(phase=phase), FULL_WL, MASK)

        assert out["phase"].dims == ("wl", "theta")
        assert out["phase"].shape == (4, 3)
        assert np.isnan(out["phase"].values[0]).all()
        assert (out["phase"].values[1] == 1.0).all()

    def test_three_dimensional_variable_is_preserved(self):
        coeff = (("wl", "l", "element"), np.ones((2, 5, 6), dtype=np.float32))
        out = reindex_to_full_grid(_kept(coeff=coeff), FULL_WL, MASK)

        assert out["coeff"].shape == (4, 5, 6)

    def test_nothing_clipped_returns_the_full_grid(self):
        ds = xr.Dataset(
            data_vars={"kext": (("wl",), np.ones(4, dtype=np.float32))},
            coords={"wl": np.asarray(FULL_WL, dtype=np.float32)},
        )

        out = reindex_to_full_grid(ds, FULL_WL, np.ones(4, dtype=bool))

        assert not np.isnan(out["kext"].values).any()
