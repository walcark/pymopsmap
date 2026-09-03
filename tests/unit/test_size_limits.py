"""Size-parameter limits read from the dataset index, not guessed."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from pymopsmap.scatlib.limits import SizeParameterLimits
from pymopsmap.shapes import Irregular, Sphere, Spheroid

MREAL = [1.40, 1.50]
MIMAG = [0.001, 0.010]
EPS = [1.0, 2.0]
NP = [-1, -7, -109]


@pytest.fixture
def limits(tmp_path) -> SizeParameterLimits:
    """A synthetic index where the spheroid limit collapses at high mimag."""
    values = np.empty((len(MIMAG), len(MREAL), len(EPS), len(NP)))
    values[...] = 1013.0
    values[..., NP.index(-109)] = 31.6
    values[MIMAG.index(0.010), :, :, NP.index(-7)] = 3.5

    path = tmp_path / "index.nc"
    xr.Dataset(
        {"max_sizepara": (("mimag", "mreal", "eps", "np"), values)},
        coords={"mimag": MIMAG, "mreal": MREAL, "eps": EPS, "np": NP},
    ).to_netcdf(path)
    return SizeParameterLimits(path)


class TestPerShape:
    def test_a_sphere_uses_the_spherical_entry(self, limits):
        assert limits.maximum(Sphere(), 1.40, 0.001) == pytest.approx(1013.0)

    def test_an_irregular_shape_has_its_own_limit(self, limits):
        assert limits.maximum(
            Irregular(shape_id="A"), 1.40, 0.001
        ) == pytest.approx(31.6)


class TestPerRefractiveIndex:
    def test_a_spheroid_limit_follows_the_imaginary_part(self, limits):
        """It collapses by a factor of nearly 300 across the grid."""
        shape = Spheroid(mode="oblate", aspect_ratio=2.0)

        assert limits.maximum(shape, 1.40, 0.001) == pytest.approx(1013.0)
        assert limits.maximum(shape, 1.40, 0.010) == pytest.approx(3.5)

    def test_the_nearest_grid_point_is_used(self, limits):
        shape = Spheroid(mode="oblate", aspect_ratio=2.0)

        assert limits.maximum(shape, 1.41, 0.0099) == pytest.approx(3.5)


class TestFallback:
    def test_without_an_index_the_published_table_is_used(self):
        limits = SizeParameterLimits(None)

        assert limits.maximum(Sphere(), 1.45, 0.001) == pytest.approx(1005.0)
        assert limits.maximum(
            Irregular(shape_id="A"), 1.45, 0.001
        ) == pytest.approx(30.2)


class TestClipping:
    def test_a_collapsed_spheroid_limit_clips(self, limits):
        from pymopsmap.engine.coverage import clip_modes_to_coverage
        from pymopsmap.microparams import MicroParameters
        from pymopsmap.psd import LognormalPSD

        mode = MicroParameters(
            wavelength=[0.55],
            n_real=[1.40],
            n_imag=[0.010],
            shape=Spheroid(mode="oblate", aspect_ratio=2.0),
            psd=LognormalPSD(rm=0.1, sigma=1.6, n=1e9, rmin=0.005, rmax=1.0),
        )

        with pytest.warns(UserWarning, match="clipped"):
            _, mask = clip_modes_to_coverage([mode], limits=limits)

        assert not mask.any()

    def test_the_same_mode_passes_where_the_limit_is_high(self, limits):
        from pymopsmap.engine.coverage import clip_modes_to_coverage
        from pymopsmap.microparams import MicroParameters
        from pymopsmap.psd import LognormalPSD

        mode = MicroParameters(
            wavelength=[0.55],
            n_real=[1.40],
            n_imag=[0.001],
            shape=Spheroid(mode="oblate", aspect_ratio=2.0),
            psd=LognormalPSD(rm=0.1, sigma=1.6, n=1e9, rmin=0.005, rmax=1.0),
        )

        _, mask = clip_modes_to_coverage([mode], limits=limits)

        assert mask.all()
