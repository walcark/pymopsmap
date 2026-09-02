"""Interpolation targets are validated against the source grids."""

from __future__ import annotations

import numpy as np
import pytest

from pymopsmap.exceptions import DomainError
from pymopsmap.utils.validation import check_within_grid

GRID = np.array([0, 5, 10, 50, 90, 95])


class TestCheckWithinGrid:
    def test_value_inside_the_grid_passes(self):
        check_within_grid("relative_humidity", [0, 42.5, 95], GRID)

    def test_value_above_the_grid_raises(self):
        with pytest.raises(DomainError) as exc:
            check_within_grid("relative_humidity", 99, GRID)

        assert "relative_humidity" in str(exc.value)
        assert "99" in str(exc.value)
        assert "95" in str(exc.value)

    def test_value_below_the_grid_raises(self):
        with pytest.raises(DomainError):
            check_within_grid("wavelength", -1, GRID)

    def test_error_carries_the_axis_and_the_bounds(self):
        with pytest.raises(DomainError) as exc:
            check_within_grid(
                "relative_humidity", 99, GRID, context="sulphate"
            )

        err = exc.value
        assert err.axis == "relative_humidity"
        assert err.available == (0.0, 95.0)
        assert err.context == "sulphate"


class TestCamsExtraction:
    def test_out_of_range_humidity_raises_a_domain_error(self):
        from pymopsmap.adapters.input.cams import (
            CamsAerosol,
            CamsVersion,
            read_aerosol_microphysical_parameters,
        )

        with pytest.raises(DomainError) as exc:
            read_aerosol_microphysical_parameters(
                aerosol=CamsAerosol.SULPHATE_CAMS,
                version=CamsVersion.V49_R1,
                wl_microns=[0.55],
                rh=[0.0, 99.0],
            )

        assert "relative_humidity" in str(exc.value)

    def test_out_of_range_wavelength_raises_a_domain_error(self):
        from pymopsmap.adapters.input.cams import (
            CamsAerosol,
            CamsVersion,
            read_aerosol_microphysical_parameters,
        )

        with pytest.raises(DomainError) as exc:
            read_aerosol_microphysical_parameters(
                aerosol=CamsAerosol.SULPHATE_CAMS,
                version=CamsVersion.V49_R1,
                wl_microns=[2.5],
                rh=[50.0],
            )

        assert "wavelength" in str(exc.value)

    def test_in_range_request_still_works(self):
        from pymopsmap.adapters.input.cams import (
            CamsAerosol,
            CamsVersion,
            read_aerosol_microphysical_parameters,
        )

        sweep = read_aerosol_microphysical_parameters(
            aerosol=CamsAerosol.SULPHATE_CAMS,
            version=CamsVersion.V49_R1,
            wl_microns=[0.44, 0.55],
            rh=[0.0, 50.0, 95.0],
        )

        assert len(sweep.mixtures) == 3
        assert all(np.isfinite(sweep.mixtures[0].modes[0].n_real))
