from pymopsmap.classes.wavelength import Wavelength, wavelength
import pytest


# =====================================================================================
# 1) Coercion tests
# =====================================================================================
def test_scalar_float_is_coerced():
    wl = Wavelength(values=0.55)
    assert wl.values == [0.55]


def test_scalar_int_is_coerced():
    wl = Wavelength(values=1)
    assert wl.values == [1.0]


def test_list_is_kept():
    wl = Wavelength(values=[0.44, 0.55])
    assert wl.values == [0.44, 0.55]


def test_invalid_type_raises():
    with pytest.raises(TypeError):
        Wavelength(values={"bad": "type"})


# =====================================================================================
# 2) Validations tests
# =====================================================================================
def test_wavelength_must_be_positive():
    with pytest.raises(ValueError, match="Wavelengths must be > 0"):
        Wavelength(values=[0.55, -0.1])


def test_wavelength_must_be_sorted():
    with pytest.raises(ValueError, match="Wavelengths must be sorted ascending"):
        Wavelength(values=[0.67, 0.44])


def test_wavelength_single_element_ok():
    wl = Wavelength(values=[0.55])
    assert wl.values == [0.55]


# =====================================================================================
# 3) to_section() tests
# =====================================================================================
def test_to_section_single():
    wl = Wavelength(values=0.55)
    assert wl.to_section() == "wavelength 0.55"


def test_to_section_list():
    wl = Wavelength(values=[0.44, 0.55, 0.67])
    assert wl.to_section() == "wavelength list 0.44 0.55 0.67"


# =====================================================================================
# 4) Tests of the wrapper wavelength()
# =====================================================================================
def test_wrapper_passes_through():
    wl = Wavelength(values=[0.44])
    assert wavelength(wl) is wl


def test_wrapper_with_float():
    wl = wavelength(0.55)
    assert isinstance(wl, Wavelength)
    assert wl.values == [0.55]


def test_wrapper_with_list():
    wl = wavelength([0.44, 0.55])
    assert wl.values == [0.44, 0.55]
