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
        from pymopsmap.exceptions import OutsideCoverageError
        from pymopsmap.microparams import MicroParameters
        from pymopsmap.psd import LognormalPSD

        mode = MicroParameters(
            wavelength=[0.55],
            n_real=[1.40],
            n_imag=[0.010],
            shape=Spheroid(mode="oblate", aspect_ratio=2.0),
            psd=LognormalPSD(rm=0.1, sigma=1.6, n=1e9, rmin=0.005, rmax=1.0),
        )

        # Nothing survives the clip, so there is no run to make: MOPSMAP would
        # be handed an empty wavelength file and misreport the reason.
        with pytest.raises(OutsideCoverageError, match="outside the dataset"):
            clip_modes_to_coverage([mode], limits=limits)

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


class TestAspectRatioOfOne:
    """
    A spheroid of aspect ratio one is a sphere, and MOPSMAP files it as one.

    ``make_contributions.f90`` sends every eps within 0.001 of one to i_np = 1,
    the sphere code. There is no spheroid file at eps = 1 to read instead: the
    index carries the row, with a maximum size parameter of zero.
    """

    def test_the_resolver_asks_for_the_sphere_file(self, tmp_path):
        import xarray as xr

        from pymopsmap.scatlib.resolver import NCFileResolver

        index = tmp_path / "index.nc"
        xr.Dataset(
            coords={
                "mreal": [1.52, 1.56],
                "mimag": [0.0, 0.0043],
                "eps": [0.833, 1.0, 1.2],
            }
        ).to_netcdf(index)
        resolver = NCFileResolver(index)

        files = resolver._files_for_params(
            Spheroid(mode="prolate", aspect_ratio=1.0), 1.52, 0.0043
        )

        assert files == ["spheres/sphere_1.5200_0.004300.nc"]

    def test_the_limit_comes_from_the_sphere_row(self, limits):
        sphere = limits.maximum(Sphere(), 1.40, 0.001)
        spheroid = limits.maximum(
            Spheroid(mode="prolate", aspect_ratio=1.0), 1.40, 0.001
        )

        assert spheroid == sphere


class TestRepeatedLookups:
    def test_the_same_index_is_looked_up_once(self, limits, monkeypatch):
        """
        An ensemble built from measured particles asks the same question
        thousands of times.

        Every mode of it carries the same refractive index on the same
        wavelength grid, so the limit is one selection, not one per mode and
        wavelength. Nine thousand modes over forty-nine wavelengths used to
        mean half a million xarray selections.
        """
        import xarray as xr

        shape = Spheroid(mode="oblate", aspect_ratio=2.0)
        selections = []
        original = xr.DataArray.sel

        def counting(self, *args, **kwargs):
            selections.append(kwargs)
            return original(self, *args, **kwargs)

        monkeypatch.setattr(xr.DataArray, "sel", counting)
        limits._known.clear()

        first = limits.maximum(shape, 1.40, 0.010)
        for _ in range(50):
            assert limits.maximum(shape, 1.40, 0.010) == first

        assert len(selections) <= 3
