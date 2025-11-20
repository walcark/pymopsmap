from pymopsmap.classes import RefractiveIndex
from pathlib import Path
import pytest


# ==============================================================================
# Coercion tests
# ==============================================================================
def test_scalar_coercion_to_list():
    ri = RefractiveIndex(wl=0.55, n_real=1.5, n_imag=0.01)
    assert ri.wl == [0.55]
    assert ri.n_real == [1.5]
    assert ri.n_imag == [0.01]


def test_list_is_accepted_as_is():
    ri = RefractiveIndex(wl=[0.44, 0.55], n_real=[1.5, 1.48], n_imag=[0.01, 0.02])
    assert ri.wl == [0.44, 0.55]
    assert ri.n_real == [1.5, 1.48]
    assert ri.n_imag == [0.01, 0.02]


def test_invalid_type_raises_typeerror():
    with pytest.raises(TypeError):
        RefractiveIndex(wl=(0.44, 0.55), n_real=1.5, n_imag=0.01)


# ==============================================================================
# Constant refractive index
# ==============================================================================
def test_constant_valid():
    ri = RefractiveIndex(n_real=1.5, n_imag=0.01)
    assert ri.wl is None
    assert ri.n_real == [1.5]  # coerced
    assert ri.n_imag == [0.01]


def test_command_constant():
    ri = RefractiveIndex(n_real=1.5, n_imag=0.01)
    assert ri.command == "refrac 1.5 0.01"


# ==============================================================================
# File-based refractive index
# ==============================================================================
def test_file_based_valid():
    ri = RefractiveIndex(filename="my_ri.txt")
    assert ri.filename == "my_ri.txt"
    assert ri.command == "refrac file 'my_ri.txt'"


def test_file_based_allows_other_fields_but_ignores_them():
    ri = RefractiveIndex(wl=[0.44], n_real=[1.5], n_imag=[0.01], filename="file.txt")
    assert ri.filename == "file.txt"
    assert ri.command == "refrac file 'file.txt'"


# ==============================================================================
# Wavelength-dependent refractive index
# ==============================================================================
def test_wldependent_valid():
    ri = RefractiveIndex(wl=[0.44, 0.55], n_real=[1.5, 1.48], n_imag=[0.01, 0.02])
    assert ri.wl == [0.44, 0.55]


def test_wldependent_length_mismatch():
    with pytest.raises(ValueError):
        RefractiveIndex(wl=[0.44], n_real=[1.5, 1.48], n_imag=[0.01, 0.02])


def test_wldependent_wavelengths_sorted():
    with pytest.raises(ValueError):
        RefractiveIndex(wl=[0.55, 0.44], n_real=[1.5, 1.48], n_imag=[0.01, 0.02])


def test_wldependent_negative_wavelength():
    with pytest.raises(ValueError):
        RefractiveIndex(wl=[-0.55, 0.44], n_real=[1.5, 1.48], n_imag=[0.01, 0.02])


def test_wldependent_negative_imag_part():
    with pytest.raises(ValueError):
        RefractiveIndex(wl=[0.44, 0.55], n_real=[1.5, 1.48], n_imag=[-0.01, 0.02])


# ==============================================================================
# Temporary file generation
# ==============================================================================
def test_tempfile_is_created():
    ri = RefractiveIndex(wl=[0.44, 0.55], n_real=[1.5, 1.48], n_imag=[0.01, 0.02])

    tmpfile = ri._write_temp_file()
    assert Path(tmpfile).exists()

    with open(tmpfile) as f:
        lines = f.readlines()

    assert len(lines) == 2
    assert "0.440000" in lines[0]
    assert "1.500000" in lines[0]
    assert "0.010000" in lines[0]

    Path(tmpfile).unlink()


def test_command_wldependent_creates_file_and_returns_path():
    ri = RefractiveIndex(wl=[0.44], n_real=[1.5], n_imag=[0.01])

    out = ri.command
    assert out.startswith("refrac file '")

    # The file path extracted must exist
    filepath = out.split("'")[1]
    assert Path(filepath).exists()

    Path(filepath).unlink()
