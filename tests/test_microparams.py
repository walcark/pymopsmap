import pytest

from pydantic import ValidationError

from pymopsmap.classes import (
    FixedPSD,
    LognormalPSD,
    ModifiedGammaPSD,
    FileDefinedPSD,
    Sphere,
    Spheroid,
    MicroParameters,
)


# =====================================================================================
# Particle-size distribution tests
# =====================================================================================
@pytest.mark.parametrize(
    "radius, n",
    [
        (-1.0, 0.5),
        (0.0, 0.5),
        (1.0, -0.1),
        (-2.0, -1.0),
    ],
)
def test_FixedPSD_radius_and_n_must_be_positive(radius, n):
    with pytest.raises(ValidationError):
        FixedPSD(radius=radius, n=n)


def test_FixedPSD_radius_is_float():
    with pytest.raises(ValidationError):
        FixedPSD(radius=[1.0, 5.0], n=1.0)


def test_FixedPSD_command():
    psd = FixedPSD(radius=1.0, n=0.05)
    assert psd.command == "size 1.0 0.05"


@pytest.mark.parametrize(
    "rm, sigma, n, rmin, rmax",
    [
        (1.0, 2.0, 0, -0.01, 200.0),
        (-0.01, 2.0, 0.5, 0.01, 200.0),
        (1.0, 2.0, 0.5, 0.01, -0.01),
        (1.0, 2.0, 0.5, -0.01, 200.0),
    ],
)
def test_LognormalPSD_arguments_must_be_positive(rm, sigma, n, rmin, rmax):
    with pytest.raises(ValidationError):
        LognormalPSD(rm=rm, sigma=sigma, n=n, rmin=rmin, rmax=rmax)


def test_LognormalPSD_radius():
    with pytest.raises(ValueError, match="rmax > rmin required, got"):
        LognormalPSD(rm=0.05, sigma=1.3, n=0.05, rmin=1.0, rmax=0.5)


def test_LognormalPSD_command():
    psd = LognormalPSD(rm=0.05, sigma=1.3, n=1.0, rmin=0.01, rmax=200.0)
    assert psd.command == "size log_normal 0.05 1.3 1.0 0.01 200.0"


@pytest.mark.parametrize(
    "A, B, rmin, rmax",
    [
        (1.0, 0.0, 0.01, 200.0),
        (0.0, 1.0, 0.01, 200.0),
        (0.01, 0.5, -0.01, 200.0),
        (0.5, 0.01, 0.01, -0.2),
    ],
)
def test_ModifiedGammaPSD_arguments_must_be_positive(A, B, rmin, rmax):
    with pytest.raises(ValidationError):
        ModifiedGammaPSD(A=A, B=B, alpha=0.03, gamma=2.0, rmin=rmin, rmax=rmax)


def test_ModifiedGammaPSD_radius():
    with pytest.raises(ValueError, match="rmax > rmin required, got"):
        ModifiedGammaPSD(
            A=0.05, B=1.0, alpha=0.03, gamma=2.0, rmin=1.0, rmax=0.5
        )


def test_ModifiedGammaPSD_command():
    psd = ModifiedGammaPSD(
        A=0.05, B=1.0, alpha=0.03, gamma=2.0, rmin=0.01, rmax=0.5
    )
    assert psd.command == "size mod_gamma 0.05 0.01 0.5 0.03 1.0 2.0"


def test_FileDefinedPSD_filename_is_string():
    with pytest.raises(ValidationError):
        FileDefinedPSD(filename=1.0)
    with pytest.raises(ValidationError):
        FileDefinedPSD(filename={"a": "test.txt"})
    with pytest.raises(ValidationError):
        FileDefinedPSD(filename=("test.txt",))


def test_FileDefinedPSD_command():
    psd = FileDefinedPSD(filename="test.txt")
    assert psd.command == "size bin_file 'test.txt'"


# ==============================================================================
# Shape tests
# ==============================================================================
def test_spheroid_mode_is_oblate_or_prolate():
    with pytest.raises(ValueError):
        Spheroid(mode="unknown-mode", aspect_ratio=1.0)


def test_spheroid_aspect_ratio_is_greater_than_one():
    with pytest.raises(ValueError):
        Spheroid(mode="spheroid", aspect_ratio=0.99)


# ==============================================================================
# Kappa and density tests
# ==============================================================================
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
