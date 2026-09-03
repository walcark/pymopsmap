"""Refractive indices are checked against the dataset before running."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pymopsmap.exceptions import CoverageError, DomainError
from pymopsmap.microparams import MicroParameters
from pymopsmap.psd import LognormalPSD
from pymopsmap.scatlib.resolver import NCFileResolver
from pymopsmap.shapes import Sphere

GRID_REAL = [1.28, 1.32, 1.60, 1.64]
GRID_IMAG = [0.0, 0.000538, 0.00152, 0.147]


@pytest.fixture
def index_path(tmp_path) -> Path:
    ds = xr.Dataset(
        {
            "mreal": ("mreal", np.array(GRID_REAL)),
            "mimag": ("mimag", np.array(GRID_IMAG)),
        }
    )
    path = tmp_path / "index.nc"
    ds.to_netcdf(path)
    return path


@pytest.fixture
def resolver(index_path) -> NCFileResolver:
    return NCFileResolver(index_path)


def _mode(n_real: float, n_imag: float = 0.001) -> MicroParameters:
    return MicroParameters(
        wavelength=[0.532],
        n_real=[n_real],
        n_imag=[n_imag],
        shape=Sphere(),
        psd=LognormalPSD(rm=0.1, sigma=1.6, n=1e9, rmin=0.005, rmax=20.0),
    )


class TestGridBounds:
    def test_a_value_inside_the_grid_resolves(self, resolver):
        assert resolver.resolve([_mode(1.45)])

    def test_a_real_part_above_the_grid_raises(self, resolver):
        with pytest.raises(DomainError, match="n_real"):
            resolver.resolve([_mode(1.75)])

    def test_a_real_part_below_the_grid_raises(self, resolver):
        with pytest.raises(DomainError, match="n_real"):
            resolver.resolve([_mode(1.0)])

    def test_an_imaginary_part_above_the_grid_raises(self, resolver):
        with pytest.raises(DomainError, match="n_imag"):
            resolver.resolve([_mode(1.45, n_imag=0.5)])

    def test_the_error_names_the_available_range(self, resolver):
        with pytest.raises(DomainError) as exc:
            resolver.resolve([_mode(1.75)])

        assert exc.value.available == (1.28, 1.64)

    def test_the_check_runs_before_any_file_is_named(self, resolver):
        """MOPSMAP itself stops on an out-of-grid value, so never get there."""
        with pytest.raises(DomainError):
            resolver.resolve([_mode(1.75), _mode(1.45)])


class TestMissingFiles:
    """A grid value whose file this copy of the dataset does not ship."""

    def test_the_missing_file_is_named(self):
        from pymopsmap.scatlib.coverage import require_available

        with pytest.raises(CoverageError) as exc:
            require_available(
                missing=["spheres/sphere_1.7600_0.147000.nc"],
                modes=[_mode(1.75, n_imag=0.147)],
                source="/somewhere",
            )

        assert "sphere_1.7600_0.147000.nc" in str(exc.value)

    def test_a_value_outside_the_main_archive_points_at_the_extended_one(self):
        """Soot sits at 1.75, which only the extended archive covers."""
        from pymopsmap.scatlib.coverage import require_available

        with pytest.raises(CoverageError) as exc:
            require_available(
                missing=["spheres/sphere_1.7600_0.147000.nc"],
                modes=[_mode(1.75, n_imag=0.147)],
                source="/somewhere",
            )

        message = str(exc.value)
        assert "extended" in message
        assert "1.75" in message

    def test_a_value_inside_it_reports_an_incomplete_copy(self):
        """Nothing justifies the gap, so the dataset itself is truncated."""
        from pymopsmap.scatlib.coverage import require_available

        with pytest.raises(CoverageError, match="incomplete"):
            require_available(
                missing=["spheres/sphere_1.6400_0.147000.nc"],
                modes=[_mode(1.62, n_imag=0.147)],
                source="/somewhere",
            )

    def test_nothing_missing_passes(self):
        from pymopsmap.scatlib.coverage import require_available

        require_available(missing=[], modes=[_mode(1.45)], source="/somewhere")
