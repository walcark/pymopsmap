"""Shared helpers: a faithful stand-in for one MOPSMAP run."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from pymopsmap.engine.outputs import OUTPUT_VARIABLES, OutputType


def integrated_result(
    wl, kext: float = 1e-6, ssa: float = 0.9, **overrides
) -> xr.Dataset:
    """
    A result carrying every variable the integrated block produces.

    The contract declares them all, so a stand-in that returns a subset is
    rejected exactly as a real partial run would be.
    """
    wl = np.asarray(wl, dtype=float)
    kext_a = np.full(len(wl), kext)
    ksca_a = kext_a * ssa
    defaults = {
        "kext": kext_a,
        "ssa": np.full(len(wl), ssa),
        "ksca": ksca_a,
        "g": np.full(len(wl), 0.7),
        # The parser gives every integrated column the wavelength axis, even
        # the ones whose value does not vary with it.
        "reff": np.full(len(wl), 0.15),
        "n": np.full(len(wl), 1.0e9),
        "cross_dens": np.full(len(wl), 2.0e-6),
        "vol_dens": np.full(len(wl), 3.0e-12),
        "mass_conc": np.full(len(wl), 4.0e-9),
        "angstrom_ext": np.full(len(wl), 1.2),
        "angstrom_sca": np.full(len(wl), 1.1),
        "angstrom_abs": np.full(len(wl), 0.5),
    }
    defaults.update(overrides)

    variables = {
        name: ((("wl",), value) if np.ndim(value) else value)
        for name, value in defaults.items()
    }
    ds = xr.Dataset(variables, coords={"wl": wl})
    ds.attrs["shape_types"] = ["sphere"]
    return ds


@pytest.fixture(autouse=True)
def isolated_sweep_store(tmp_path, monkeypatch):
    """Each test sweeps into its own store, so none inherits another's."""
    monkeypatch.setenv("PYMOPSMAP_SWEEP_STORE", str(tmp_path / "sweeps"))


@pytest.fixture
def integrated():
    """The helper, as a fixture, for tests that stub the engine."""
    return integrated_result


def variables_of(*families: OutputType) -> list[str]:
    names: list[str] = []
    for family in families:
        names.extend(OUTPUT_VARIABLES[family])
    return names
