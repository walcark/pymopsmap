import pytest

from pymopsmap.classes import (
    FixedPSD,
    LognormalPSD,
    ModifiedGammaPSD,
    FileDefinedPSD,
    Sphere,
    Spheroid,
    MicroParameters,
)
from pymopsmap.mopsmap.commands import psd_command, shape_command


def test_FixedPSD_command():
    psd = FixedPSD(radius=1.0, n=0.05)
    assert psd_command(psd) == "size 1.0 0.05"


def test_LognormalPSD_command():
    psd = LognormalPSD(rm=0.05, sigma=1.3, n=1.0, rmin=0.01, rmax=200.0)
    assert psd_command(psd) == "size log_normal 0.05 1.3 1.0 0.01 200.0"


def test_ModifiedGammaPSD_command():
    psd = ModifiedGammaPSD(
        A=0.05, B=1.0, alpha=0.03, gamma=2.0, rmin=0.01, rmax=0.5
    )
    assert psd_command(psd) == "size mod_gamma 0.05 0.01 0.5 0.03 1.0 2.0"


def test_FileDefinedPSD_command():
    psd = FileDefinedPSD(filename="test.txt")
    assert psd_command(psd) == "size bin_file 'test.txt'"


def test_sphere_command():
    sphere = Sphere()
    assert shape_command(sphere) == "shape sphere"


def test_spheroid_command():
    for mode in ["prolate", "oblate"]:
        spheroid = Spheroid(mode=mode, aspect_ratio=1.3)
        assert shape_command(spheroid) == f"shape spheroid {mode} 1.3"


@pytest.mark.parametrize(
    "k, d",
    [(0.0, 1.6), (-0.01, 1.6), (0.01, 0.0), (0.01, -0.01)],
)
def test_kappa_and_density_are_positive(k, d):
    with pytest.raises(ValueError):
        MicroParameters(
            shape=Sphere(),
            psd=LognormalPSD(rm=0.05, sigma=1.3, n=1.0, rmin=0.01, rmax=200.0),
            wavelength=[0.330, 0.550, 0.870],
            n_real=[1.0, 1.5, 2.0],
            n_imag=[1e-4, 1e-5, 1e-5],
            kappa=k,
            density=d,
        )
