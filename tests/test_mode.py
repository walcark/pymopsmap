from pymopsmap.classes import Mode, Sphere, RefractiveIndex, LognormalPSD
from textwrap import dedent
import pytest


# =================================================================================================
# Build base objects for test
# =================================================================================================
shape = Sphere()
psd = LognormalPSD(rm=0.2, sigma=1.6, n=1.0, rmin=0.01, rmax=200.0)
refr_index = RefractiveIndex(filename="test.txt")


# =================================================================================================
# Test positivity of (kappa, density) and command
# =================================================================================================
def test_kappa_and_density_are_positive():
    with pytest.raises(ValueError):
        Mode(shape=shape, psd=psd, refr_index=refr_index, kappa=-0.01)
    with pytest.raises(ValueError):
        Mode(shape=shape, psd=psd, refr_index=refr_index, kappa=-0.01, density=1.6)
    with pytest.raises(ValueError):
        Mode(shape=shape, psd=psd, refr_index=refr_index, kappa=0.01, density=-1.6)
    with pytest.raises(ValueError):
        Mode(shape=shape, psd=psd, refr_index=refr_index, density=-1.6)


def test_mode_command():
    mode = Mode(shape=shape, psd=psd, refr_index=refr_index)
    expected = dedent("""
        mode 3 shape sphere
        mode 3 size log_normal 0.2 1.6 1.0 0.01 200.0
        mode 3 refrac file 'test.txt'
    """).strip()
    assert mode.command(3) == expected

    mode = Mode(shape=shape, psd=psd, refr_index=refr_index, kappa=0.05, density=1.7)
    expected = dedent("""
        mode 2 shape sphere
        mode 2 size log_normal 0.2 1.6 1.0 0.01 200.0
        mode 2 refrac file 'test.txt'
        mode 2 kappa 0.05
        mode 2 density 1.7
    """).strip()
    assert mode.command(2) == expected
