"""Integration tests for end-to-end compute pipeline.

These tests require:
  - MOPSMAP binary at bin/mopsmap/mopsmap
  - PYMOPSMAP_DATASET_SOURCE set to a local dataset dir
  - The required .nc files present (or downloadable)

Run with:
  pixi run -e dev pytest tests/integration/test_compute.py
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
    def test_basic_run_returns_optiprops(self):
        from pymopsmap.engine import run_point
        from pymopsmap.models.output_request import DEFAULT_OUTPUT

        op = run_point([_sphere_mp()], DEFAULT_OUTPUT, quiet=True)
        kext = op["kext"]
        assert kext.dims == ("wl",)
        assert len(kext) == 3
        assert (kext.values > 0).all()

    def test_result_cache_hit(self, tmp_path):
        from pymopsmap.engine import run_point
        from pymopsmap.models.output_request import DEFAULT_OUTPUT
        from pymopsmap.scatlib.results import ResultCache

        mp = _sphere_mp()
        run_point([mp], DEFAULT_OUTPUT, quiet=True)

        rc = ResultCache()
        assert rc.get(rc.key([mp], DEFAULT_OUTPUT)) is not None

    def test_invalid_params_raises_before_io(self):
        from pydantic import ValidationError

        from pymopsmap.models.microparams import (
            LognormalPSD,
            MicroParameters,
            Sphere,
        )

        with pytest.raises(ValidationError):
            MicroParameters(
                wavelength=[0.55],
                n_real=1.4,
                n_imag=1e-4,
                shape=Sphere(),
                psd=LognormalPSD(
                    rm=0.1, sigma=2.0, n=1e6, rmin=10.0, rmax=1.0
                ),
            )
