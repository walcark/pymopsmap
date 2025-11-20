from pymopsmap.classes import Sphere, Spheroid
import pytest


# ==============================================================================
# Sphere tests
# ==============================================================================
def test_sphere_command():
    sphere = Sphere()
    assert sphere.command == "shape sphere"


# ==============================================================================
# Spheroïd tests
# ==============================================================================
def test_spheroid_mode_is_oblate_or_prolate():
    with pytest.raises(ValueError):
        Spheroid(mode="unknown-mode", aspect_ratio=1.0)


def test_spheroid_aspect_ratio_is_greater_than_one():
    with pytest.raises(ValueError):
        Spheroid(mode="spheroid", aspect_ratio=0.99)


def test_spheroid_command():
    for mode in ["prolate", "oblate"]:
        spheroid = Spheroid(mode=mode, aspect_ratio=1.3)
        assert spheroid.command == f"shape spheroid {mode} 1.3"
