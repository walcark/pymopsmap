"""Unit tests for NCFileResolver."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pymopsmap.models.microparams import (
    Irregular,
    LognormalPSD,
    MicroParameters,
    Sphere,
)
from pymopsmap.scatlib.resolver import NCFileResolver


def _make_index(mreal_vals, mimag_vals, eps_vals=None, tmp_path=None) -> Path:
    """Create a minimal synthetic index.nc fixture."""
    ds = xr.Dataset(
        {
            "mreal": ("mreal", np.array(mreal_vals, dtype=np.float64)),
            "mimag": ("mimag", np.array(mimag_vals, dtype=np.float64)),
        }
    )
    if eps_vals is not None:
        ds["eps"] = ("eps", np.array(eps_vals, dtype=np.float64))
    if tmp_path is None:
        tmp_path = Path(tempfile.mkdtemp())
    p = tmp_path / "index.nc"
    ds.to_netcdf(p)
    return p


@pytest.fixture
def index_path(tmp_path):
    return _make_index(
        mreal_vals=[1.3, 1.4, 1.5, 1.6],
        mimag_vals=[0.000001, 0.00001, 0.0001, 0.001],
        eps_vals=[1.0, 1.5, 2.0],
        tmp_path=tmp_path,
    )


@pytest.fixture
def resolver(index_path):
    return NCFileResolver(index_path)


class TestSphericalFilenames:
    def test_within_grid(self, resolver):
        mp = MicroParameters(
            wavelength=[0.55],
            n_real=1.45,
            n_imag=5e-5,
            shape=Sphere(),
            psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
        )
        files = resolver.resolve(mp)
        assert "index.nc" in files
        sphere_files = [f for f in files if f.startswith("spheres/")]
        assert len(sphere_files) == 4  # 2 mreal × 2 mimag brackets

    def test_exactly_on_grid_point(self, resolver):
        # mreal=1.4, mimag=0.001 are exact grid points in the fixture
        mp = MicroParameters(
            wavelength=[0.55],
            n_real=1.4,
            n_imag=0.001,
            shape=Sphere(),
            psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
        )
        files = resolver.resolve(mp)
        sphere_files = [f for f in files if f.startswith("spheres/")]
        # Exact grid points collapse to 1 value each → 1 × 1 = 1 file
        assert len(sphere_files) == 1

    def test_index_always_included(self, resolver):
        mp = MicroParameters(
            wavelength=[0.55],
            n_real=1.45,
            n_imag=5e-5,
            shape=Sphere(),
            psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
        )
        assert "index.nc" in resolver.resolve(mp)

    def test_filename_format(self, resolver):
        mp = MicroParameters(
            wavelength=[0.55],
            n_real=1.4,
            n_imag=0.0001,
            shape=Sphere(),
            psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
        )
        sphere_files = [
            f for f in resolver.resolve(mp) if f.startswith("spheres/")
        ]
        for f in sphere_files:
            assert f.startswith("spheres/sphere_")
            assert f.endswith(".nc")


class TestIrregularFilenames:
    def test_irregular_shape_id(self, resolver):
        mp = MicroParameters(
            wavelength=[0.55],
            n_real=1.45,
            n_imag=5e-5,
            shape=Irregular(shape_id="B"),
            psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
        )
        files = resolver.resolve(mp)
        irr_files = [f for f in files if f.startswith("irregular/")]
        assert all("shapeB" in f for f in irr_files)


class TestMultipleWavelengths:
    def test_deduplication_across_wavelengths(self, resolver):
        # Same n_real/n_imag → same files regardless of wavelength count
        mp_two = MicroParameters(
            wavelength=[0.55, 0.67],
            n_real=[1.45, 1.45],
            n_imag=[5e-5, 5e-5],
            shape=Sphere(),
            psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
        )
        mp_one = MicroParameters(
            wavelength=[0.55],
            n_real=1.45,
            n_imag=5e-5,
            shape=Sphere(),
            psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
        )
        # Files from two identical (n_real, n_imag) pairs should deduplicate
        files_two = resolver.resolve(mp_two)
        files_one = resolver.resolve(mp_one)
        assert sorted(files_two) == sorted(files_one)
