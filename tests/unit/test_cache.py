"""Unit tests for OpticalDatasetCache, ResultCache, and cache_key."""

import numpy as np
import xarray as xr

from pymopsmap.engine.outputs import DEFAULT_OUTPUT, OutputType
from pymopsmap.microparams import MicroParameters
from pymopsmap.psd import FixedPSD
from pymopsmap.scatlib.cache import OpticalDatasetCache
from pymopsmap.scatlib.results import ResultCache
from pymopsmap.shapes import Sphere
from pymopsmap.utils.caching import cache_key

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_result() -> xr.Dataset:
    wl = np.array([0.55, 0.67], dtype=np.float32)
    ds = xr.Dataset(
        data_vars={
            "kext": (("wl",), np.array([1e-3, 2e-3], dtype=np.float32))
        },
        coords={"wl": wl},
    )
    return ds


def _make_mp(**kwargs) -> MicroParameters:
    defaults = dict(
        wavelength=[0.55],
        n_real=1.4,
        n_imag=0.001,
        shape=Sphere(),
        psd=FixedPSD(radius=0.1, n=1e6),
    )
    defaults.update(kwargs)
    return MicroParameters(**defaults)


# ---------------------------------------------------------------------------
# OpticalDatasetCache
# ---------------------------------------------------------------------------


class TestOpticalDatasetCache:
    def test_is_cached_missing(self, tmp_path):
        c = OpticalDatasetCache(root_dir=tmp_path)
        assert not c.is_cached("spheres/test.nc")

    def test_full_path(self, tmp_path):
        c = OpticalDatasetCache(root_dir=tmp_path)
        assert (
            c.full_path("spheres/test.nc") == tmp_path / "spheres" / "test.nc"
        )

    def test_register_atomic_write(self, tmp_path):
        c = OpticalDatasetCache(root_dir=tmp_path)
        tmp_file = tmp_path / "transfer.tmp"
        tmp_file.write_bytes(b"nc content")
        c.register("spheres/test.nc", tmp_file)
        assert c.is_cached("spheres/test.nc")
        assert not tmp_file.exists()

    def test_is_cached_empty_file(self, tmp_path):
        c = OpticalDatasetCache(root_dir=tmp_path)
        (tmp_path / "empty.nc").touch()
        assert not c.is_cached("empty.nc")

    def test_cached_files_lists_nc(self, tmp_path):
        c = OpticalDatasetCache(root_dir=tmp_path)
        (tmp_path / "a.nc").write_bytes(b"x")
        files = c.cached_files()
        assert "a.nc" in files


# ---------------------------------------------------------------------------
# ResultCache
# ---------------------------------------------------------------------------


class TestResultCache:
    def test_miss_returns_none(self, tmp_path):
        rc = ResultCache(root_dir=tmp_path)
        mp = _make_mp()
        k = rc.key(mp, DEFAULT_OUTPUT)
        assert rc.get(k) is None

    def test_put_and_get(self, tmp_path):
        rc = ResultCache(root_dir=tmp_path)
        mp = _make_mp()
        op = _make_result()
        k = rc.key(mp, DEFAULT_OUTPUT)
        rc.put(k, op)
        loaded = rc.get(k)
        assert loaded is not None
        np.testing.assert_allclose(loaded["kext"].values, op["kext"].values)

    def test_key_determinism(self, tmp_path):
        rc = ResultCache(root_dir=tmp_path)
        mp = _make_mp()
        k1 = rc.key(mp, DEFAULT_OUTPUT)
        k2 = rc.key(mp, DEFAULT_OUTPUT)
        assert k1 == k2

    def test_different_output_types_give_different_keys(self, tmp_path):
        rc = ResultCache(root_dir=tmp_path)
        mp = _make_mp()
        k1 = rc.key(mp, DEFAULT_OUTPUT)
        k2 = rc.key(mp, frozenset({OutputType.INTEGRATED, OutputType.LIDAR}))
        assert k1 != k2

    def test_different_params_give_different_keys(self, tmp_path):
        rc = ResultCache(root_dir=tmp_path)
        mp1 = _make_mp(n_real=1.4)
        mp2 = _make_mp(n_real=1.5)
        k1 = rc.key(mp1, DEFAULT_OUTPUT)
        k2 = rc.key(mp2, DEFAULT_OUTPUT)
        assert k1 != k2

    def test_atomic_write(self, tmp_path):
        rc = ResultCache(root_dir=tmp_path)
        mp = _make_mp()
        op = _make_result()
        k = rc.key(mp, DEFAULT_OUTPUT)
        rc.put(k, op)
        # tmp file must be gone
        assert not (tmp_path / f"{k}.tmp").exists()
        assert (tmp_path / f"{k}.nc").exists()


# ---------------------------------------------------------------------------
# cache_key utility
# ---------------------------------------------------------------------------


class TestCacheKey:
    def test_deterministic(self):
        mp = _make_mp()
        assert cache_key(mp) == cache_key(mp)

    def test_different_models_differ(self):
        mp1 = _make_mp(n_real=1.4)
        mp2 = _make_mp(n_real=1.5)
        assert cache_key(mp1) != cache_key(mp2)

    def test_extra_changes_key(self):
        mp = _make_mp()
        k1 = cache_key(mp)
        k2 = cache_key(mp, extra="something")
        assert k1 != k2
