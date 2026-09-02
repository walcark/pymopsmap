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
