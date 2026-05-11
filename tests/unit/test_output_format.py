"""Unit tests for MOPSMAP output format parsers."""

import numpy as np
import pytest

from pymopsmap.engine.output_format import (
    format_stdout,
    parse_coeff,
    parse_lidar,
    parse_phase_function,
    parse_scattering_matrix,
    parse_volume_scattering_function,
)

# ---------------------------------------------------------------------------
# Fixtures: sample ASCII output strings matching MOPSMAP format
# ---------------------------------------------------------------------------

STDOUT_SAMPLE = """\
some header
integrated 0.550000 1.23e-3 0.85 0.70 0.15 1e6 2e-7 1e-5 1e-7 1.2 1.1 0.5
integrated 0.670000 2.34e-3 0.80 0.65 0.18 9e5 2.1e-7 1.1e-5 1.1e-7 1.3 1.2 0.6
"""

PHASE_FUNCTION_SAMPLE = """\
0.550000 0.0   4.500000
0.550000 90.0  0.700000
0.550000 180.0 0.300000
0.670000 0.0   4.200000
0.670000 90.0  0.650000
0.670000 180.0 0.280000
"""

SCATTERING_MATRIX_SAMPLE = """\
0.550000 0.0   4.5 1.1 1.0 1.0 0.1 0.05
0.550000 90.0  0.7 0.3 0.3 0.3 0.0 0.00
0.670000 0.0   4.2 1.0 1.0 1.0 0.1 0.04
0.670000 90.0  0.6 0.3 0.3 0.3 0.0 0.00
"""

VOL_SCA_SAMPLE = """\
0.550000 0.0   1.2e-5
0.550000 90.0  2.3e-6
0.670000 0.0   1.1e-5
0.670000 90.0  2.0e-6
"""

LIDAR_SAMPLE = """\
0.550000 1.23e-3 4.56e-5 27.0 0.03 1.2 1.1 1.5e3 6.5e-3
0.670000 2.34e-3 5.67e-5 28.0 0.04 1.3 1.2 1.4e3 6.0e-3
"""

COEFF_SAMPLE = """\
0.550000 0 1.5 0.5 0.5 0.5 0.1 0.05
0.550000 1 0.8 0.3 0.3 0.3 0.0 0.00
0.670000 0 1.4 0.5 0.5 0.5 0.1 0.04
0.670000 1 0.7 0.3 0.3 0.3 0.0 0.00
"""


# ---------------------------------------------------------------------------
# format_stdout
# ---------------------------------------------------------------------------


class TestFormatStdout:
    def test_basic_parse(self):
        ds = format_stdout(STDOUT_SAMPLE)
        assert "kext" in ds
        assert len(ds["wl"]) == 2

    def test_variable_names(self):
        ds = format_stdout(STDOUT_SAMPLE)
        expected = {
            "kext",
            "ssa",
            "g",
            "reff",
            "n",
            "cross_dens",
            "vol_dens",
            "mass_conc",
            "angstrom_ext",
            "angstrom_sca",
            "angstrom_abs",
        }
        assert expected.issubset(set(ds.data_vars))

    def test_wavelength_values(self):
        ds = format_stdout(STDOUT_SAMPLE)
        np.testing.assert_allclose(ds["wl"].values, [0.55, 0.67], rtol=1e-3)

    def test_empty_stdout_raises(self):
        with pytest.raises(ValueError):
            format_stdout("no integrated blocks here")


# ---------------------------------------------------------------------------
# parse_phase_function
# ---------------------------------------------------------------------------


class TestParsePhaseFunction:
    def test_shape(self):
        ds = parse_phase_function(PHASE_FUNCTION_SAMPLE)
        assert ds["phase"].dims == ("wl", "theta")
        assert ds["phase"].shape == (2, 3)

    def test_coords(self):
        ds = parse_phase_function(PHASE_FUNCTION_SAMPLE)
        assert "wl" in ds.coords
        assert "theta" in ds.coords

    def test_values(self):
        ds = parse_phase_function(PHASE_FUNCTION_SAMPLE)
        val = float(ds["phase"].sel(wl=0.55, theta=0.0))
        assert val == pytest.approx(4.5, rel=1e-4)


# ---------------------------------------------------------------------------
# parse_scattering_matrix
# ---------------------------------------------------------------------------


class TestParseScatteringMatrix:
    def test_shape(self):
        ds = parse_scattering_matrix(SCATTERING_MATRIX_SAMPLE)
        assert ds["scattering_matrix"].dims == ("wl", "theta", "element")
        assert ds["scattering_matrix"].shape == (2, 2, 6)

    def test_element_coord(self):
        ds = parse_scattering_matrix(SCATTERING_MATRIX_SAMPLE)
        assert list(ds["element"].values) == [
            "a1",
            "a2",
            "a3",
            "a4",
            "b1",
            "b2",
        ]


# ---------------------------------------------------------------------------
# parse_volume_scattering_function
# ---------------------------------------------------------------------------


class TestParseVolumeSca:
    def test_shape(self):
        ds = parse_volume_scattering_function(VOL_SCA_SAMPLE)
        assert "vol_sca_func" in ds
        assert ds["vol_sca_func"].shape == (2, 2)


# ---------------------------------------------------------------------------
# parse_lidar
# ---------------------------------------------------------------------------


class TestParseLidar:
    def test_variables(self):
        ds = parse_lidar(LIDAR_SAMPLE)
        expected = {
            "backscatter",
            "lidar_ratio",
            "depol_ratio",
            "angstrom_back",
            "angstrom_ext",
            "ext_to_mass",
            "back_to_mass",
        }
        assert expected.issubset(set(ds.data_vars))

    def test_shape(self):
        ds = parse_lidar(LIDAR_SAMPLE)
        assert len(ds["wl"]) == 2

    def test_lidar_ratio_value(self):
        ds = parse_lidar(LIDAR_SAMPLE)
        assert float(ds["lidar_ratio"].isel(wl=0)) == pytest.approx(
            27.0, rel=1e-4
        )


# ---------------------------------------------------------------------------
# parse_coeff
# ---------------------------------------------------------------------------


class TestParseCoeff:
    def test_shape(self):
        ds = parse_coeff(COEFF_SAMPLE)
        assert ds["coeff"].dims == ("wl", "l", "coeff_element")
        assert ds["coeff"].shape == (2, 2, 6)

    def test_coeff_element_names(self):
        ds = parse_coeff(COEFF_SAMPLE)
        assert list(ds["coeff_element"].values) == [
            "alpha1",
            "alpha2",
            "alpha3",
            "alpha4",
            "beta1",
            "beta2",
        ]
