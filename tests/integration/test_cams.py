"""Integration tests for CAMS aerosol adapter.

Run with:
  DATA_PATH=tests/data pixi run -e test pytest tests/integration/test_cams.py
"""

import os

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    os.getenv("PYMOPSMAP_DATASET_SOURCE") is None,
    reason="PYMOPSMAP_DATASET_SOURCE not set — skipping integration tests",
)


class TestCamsContinenal:
    def test_cams_kext_shape(self):
        from pymopsmap.adapters import CamsAerosol, CamsVersion, cams_to_kext

        kext = cams_to_kext(
            aerosol=CamsAerosol.CONTINENTAL,
            version=CamsVersion.V49_R1,
            wl_microns=list(np.linspace(0.4, 0.8, 10)),
            rh=[0.0, 50.0, 80.0],
            quiet=True,
        )
        assert "wl" in kext.dims
        assert (kext.values > 0).all()

    def test_unsupported_combination_raises(self):
        from pymopsmap.adapters import cams_to_optiprops

        with pytest.raises(Exception):
            cams_to_optiprops(
                aerosol="NONEXISTENT_SPECIE",
                version="v99",
                wl_microns=[0.55],
                rh=[0.0],
                quiet=True,
            )
