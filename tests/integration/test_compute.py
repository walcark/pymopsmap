"""Integration tests for end-to-end compute pipeline.

These tests require:
  - MOPSMAP binary at bin/mopsmap/mopsmap
  - PYMOPSMAP_DATASET_SOURCE set to a local dataset dir
  - The required .nc files present (or downloadable)

Run with:
  DATA_PATH=tests/data pixi run -e test pytest tests/integration/test_compute.py
"""

import os

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    os.getenv("PYMOPSMAP_DATASET_SOURCE") is None,
    reason="PYMOPSMAP_DATASET_SOURCE not set — skipping integration tests",
)


def _sphere_mp(wl_count=3):
    from pymopsmap import LognormalPSD, MicroParameters, Sphere

    return MicroParameters(
        wavelength=list(np.linspace(0.4, 0.8, wl_count)),
        n_real=1.45,
        n_imag=1e-4,
        shape=Sphere(),
        psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
    )


class TestSingleSphere:
    def test_basic_compute_returns_optiprops(self):
        from pymopsmap import compute

        mp = _sphere_mp()
        op = compute(mp, quiet=True)
        kext = op.sel("kext")
        assert kext.dims == ("wl",)
        assert len(kext) == 3
        assert (kext.values > 0).all()

    def test_result_cache_hit(self, tmp_path):
        from pymopsmap import compute
        from pymopsmap.cache.results import ResultCache
        from pymopsmap.models.output_request import DEFAULT_OUTPUT

        mp = _sphere_mp()
        compute(mp, quiet=True)

        rc = ResultCache()
        k = rc.key(mp, DEFAULT_OUTPUT)
        assert rc.get(k) is not None

    def test_invalid_params_raises_before_io(self):
        from pydantic import ValidationError

        from pymopsmap.models.microparams import LognormalPSD, MicroParameters, Sphere

        with pytest.raises(ValidationError):
            MicroParameters(
                wavelength=[0.55],
                n_real=1.4,
                n_imag=1e-4,
                shape=Sphere(),
                psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=10.0, rmax=1.0),
            )


class TestDispatchRh:
    def test_dispatch_returns_rh_dimension(self):
        from pymopsmap import MicroParameters, MicroParametersDispatch, Sphere, compute
        from pymopsmap.models.microparams import LognormalPSD

        dispatch = MicroParametersDispatch()
        for rh in [0.0, 50.0, 80.0]:
            mp = MicroParameters(
                wavelength=[0.55],
                n_real=1.45,
                n_imag=1e-4,
                shape=Sphere(),
                psd=LognormalPSD(rm=0.1, sigma=2.0, n=1e6, rmin=0.005, rmax=20.0),
                kappa=0.3,
            )
            dispatch.append(modes=mp, params={"rh": rh})

        op = compute(dispatch, quiet=True)
        assert "rh" in op.ds.dims or "rh" in op.ds.coords
