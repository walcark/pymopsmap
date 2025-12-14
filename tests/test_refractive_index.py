from pathlib import Path

import pytest
from pydantic import ValidationError

from pymopsmap.classes import RefractiveIndex


# ==============================================================================
# Coercion tests
# ==============================================================================
def test_scalar_coercion_to_list():
    ri = RefractiveIndex(wl=0.55, n_real=1.5, n_imag=0.01)
    assert ri.wl == [0.55]
    assert ri.n_real == [1.5]
    assert ri.n_imag == [0.01]


def test_list_is_accepted():
    ri = RefractiveIndex(
        wl=[0.44, 0.55], n_real=[1.5, 1.48], n_imag=[0.01, 0.02]
    )
    assert ri.wl == [0.44, 0.55]
    assert ri.n_real == [1.5, 1.48]
    assert ri.n_imag == [0.01, 0.02]


def test_wl_not_specified():
    with pytest.raises(ValueError):
        RefractiveIndex(n_real=1.5, n_imag=0.01)


def test_invalid_type_raises_typeerror():
    with pytest.raises(ValidationError):
        RefractiveIndex(wl=(0.44, 0.55), n_real=1.5, n_imag=0.01)


# ==============================================================================
# Constant refractive index
# ==============================================================================
def test_constant_valid():
    ri = RefractiveIndex(wl=0.55, n_real=1.5, n_imag=0.01)
    assert ri.wl == [0.55]
    assert ri.n_real == [1.5]
    assert ri.n_imag == [0.01]


# ==============================================================================
# Test file content
# ==============================================================================
def test_file_based_valid():
    ri = RefractiveIndex(filename="my_ri.txt")
    assert ri.filename == "my_ri.txt"
    assert ri.command == "refrac file 'my_ri.txt'"


def test_file_based_allows_other_fields_but_ignores_them():
    ri = RefractiveIndex(
        wl=[0.44], n_real=[1.5], n_imag=[0.01], filename="file.txt"
    )
    assert ri.filename == "file.txt"
    assert ri.command == "refrac file 'file.txt'"


def test_arg_based_check_file_name_and_content():
    ri = RefractiveIndex(
        wl=[0.44, 0.55], n_real=[1.5, 1.48], n_imag=[0.01, 0.02]
    )
    filepath: Path = Path(ri.filename)

    assert filepath.exists()
    assert filepath.name == "ri.txt"

    command_filepath = Path(ri.command.split(sep=" ")[2][1:-1])
    assert command_filepath.exists()

    with open(filepath) as f:
        lines = f.readlines()

    assert len(lines) == 2
    assert "0.440000" in lines[0]
    assert "1.500000" in lines[0]
    assert "0.010000" in lines[0]

    filepath.unlink()


# ==============================================================================
# Wavelength-dependent refractive index
# ==============================================================================
def test_wldependent_valid():
    ri = RefractiveIndex(
        wl=[0.44, 0.55], n_real=[1.5, 1.48], n_imag=[0.01, 0.02]
    )
    assert ri.wl == [0.44, 0.55]


def test_wldependent_length_mismatch():
    with pytest.raises(ValueError):
        RefractiveIndex(wl=[0.44], n_real=[1.5, 1.48], n_imag=[0.01, 0.02])


def test_wldependent_wavelengths_sorted():
    with pytest.raises(ValueError):
        RefractiveIndex(
            wl=[0.55, 0.44], n_real=[1.5, 1.48], n_imag=[0.01, 0.02]
        )


def test_wldependent_negative_wavelength():
    with pytest.raises(ValueError):
        RefractiveIndex(
            wl=[-0.55, 0.44], n_real=[1.5, 1.48], n_imag=[0.01, 0.02]
        )


def test_wldependent_negative_imag_part():
    with pytest.raises(ValueError):
        RefractiveIndex(
            wl=[0.44, 0.55], n_real=[1.5, 1.48], n_imag=[-0.01, 0.02]
        )
