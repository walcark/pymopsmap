"""Dataset files are resolved for the index MOPSMAP will actually use."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pymopsmap.models.microparams import LognormalPSD, MicroParameters, Sphere
from pymopsmap.scatlib.resolver import NCFileResolver

GRID_REAL = [1.36, 1.40, 1.44, 1.48, 1.52, 1.56, 1.60]
GRID_IMAG = [
    0.001075,
    0.00152,
    0.00215,
    0.003041,
    0.004300,
    0.006081,
    0.0086,
    0.012162,
]


@pytest.fixture
def resolver(tmp_path) -> NCFileResolver:
    path: Path = tmp_path / "index.nc"
    xr.Dataset(
        {
            "mreal": ("mreal", np.array(GRID_REAL)),
            "mimag": ("mimag", np.array(GRID_IMAG)),
        }
    ).to_netcdf(path)
    return NCFileResolver(path)


def _waso(kappa: float | None = 0.249) -> MicroParameters:
    """Water soluble at 532 nm, dry."""
    return MicroParameters(
        wavelength=[0.532],
        n_real=[1.530],
        n_imag=[0.0086],
        shape=Sphere(),
        psd=LognormalPSD(rm=0.0212, sigma=2.24, n=1e9, rmin=0.005, rmax=20.0),
        kappa=kappa,
    )


def _reals(files: list[str]) -> set[str]:
    return {f.split("_")[1] for f in files if f.startswith("spheres/")}


class TestDryModes:
    def test_a_dry_mode_resolves_around_its_own_index(self, resolver):
        files = resolver.resolve([_waso(kappa=None)])

        assert _reals(files) == {"1.5200", "1.5600"}

    def test_a_kappa_mode_without_humidity_stays_dry(self, resolver):
        files = resolver.resolve([_waso()], rh=None)

        assert _reals(files) == {"1.5200", "1.5600"}


class TestGrownModes:
    def test_growth_moves_the_index_the_engine_reads(self, resolver):
        """
        MOPSMAP mixes the dry index with water before computing, so the files
        it opens are not the ones the dry index points at. waso at 50 percent
        lands on 1.4907, not 1.530.
        """
        files = resolver.resolve([_waso()], rh=50.0)

        assert _reals(files) == {"1.4800", "1.5200"}

    def test_both_states_are_covered_when_they_differ(self, resolver):
        files_dry = resolver.resolve([_waso()], rh=None)
        files_wet = resolver.resolve([_waso()], rh=50.0)

        assert files_dry != files_wet

    def test_a_zero_kappa_mode_does_not_move(self, resolver):
        """Soot takes up no water, so humidity leaves it where it is."""
        files = resolver.resolve([_waso(kappa=0.0)], rh=90.0)

        assert _reals(files) == {"1.5200", "1.5600"}

    def test_stronger_growth_moves_further(self, resolver):
        low = resolver.resolve([_waso()], rh=50.0)
        high = resolver.resolve([_waso()], rh=90.0)

        assert min(_reals(high)) < min(_reals(low))
