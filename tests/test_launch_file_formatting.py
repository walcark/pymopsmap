import re

import pytest

from pymopsmap.engine.launch_file import write_launching_file
from pymopsmap.microparams import MicroParameters
from pymopsmap.psd import FixedPSD
from pymopsmap.shapes import Sphere
from pymopsmap.utils import DATASET_CACHE_DIR


def _normalize(content: str) -> str:
    tmp_pattern = r"/tmp/pymopsmap-[^/]+"
    content = re.sub(tmp_pattern, "/tmp/pymopsmap-TMP", content)
    return "\n".join(line.rstrip() for line in content.splitlines()).strip()


@pytest.fixture
def sphere_mp():
    return MicroParameters(
        wavelength=[0.500, 0.700, 0.900, 1.5],
        n_real=[1.00027] * 4,
        n_imag=[0.0] * 4,
        shape=Sphere(),
        psd=FixedPSD(radius=0.1, n=1.0),
    )


def test_launch_file_contains_scatlib(sphere_mp):
    paths = write_launching_file(sphere_mp)
    content = paths["mopsmap"].read_text()
    assert f"scatlib '{DATASET_CACHE_DIR}'" in content


def test_launch_file_contains_integrated(sphere_mp):
    paths = write_launching_file(sphere_mp)
    content = paths["mopsmap"].read_text()
    assert "output integrated" in content


def test_launch_file_with_rh(sphere_mp):
    paths = write_launching_file(sphere_mp, rh=50.0)
    content = paths["mopsmap"].read_text()
    assert "rH 50.0" in content


def test_launch_file_without_rh(sphere_mp):
    paths = write_launching_file(sphere_mp)
    content = paths["mopsmap"].read_text()
    assert "rH" not in content


def test_launch_file_no_ascii_for_integrated_only(sphere_mp):
    from pymopsmap.engine.outputs import DEFAULT_OUTPUT

    paths = write_launching_file(sphere_mp, output_types=DEFAULT_OUTPUT)
    content = paths["mopsmap"].read_text()
    assert "ascii_file" not in content


def test_launch_file_ascii_for_lidar(sphere_mp):
    from pymopsmap.engine.outputs import OutputType

    paths = write_launching_file(
        sphere_mp,
        output_types=frozenset({OutputType.INTEGRATED, OutputType.LIDAR}),
    )
    content = paths["mopsmap"].read_text()
    assert "ascii_file" in content
    assert "output lidar" in content


@pytest.fixture
def single_wl_mp():
    return MicroParameters(
        wavelength=[0.550],
        n_real=1.45,
        n_imag=1e-4,
        shape=Sphere(),
        psd=FixedPSD(radius=0.1, n=1.0),
    )


def test_single_wavelength_uses_wavelength_value_command(single_wl_mp):
    """Single wl: use 'wavelength <val>' (MOPSMAP interpolate_linear bug)."""
    paths = write_launching_file(single_wl_mp)
    content = paths["mopsmap"].read_text()
    match = re.search(r"^wavelength (\S+)$", content, re.MULTILINE)
    assert match is not None
    assert float(match.group(1)) == pytest.approx(0.55, rel=1e-9)
    assert "from_refrac_file" not in content


def test_single_wavelength_uses_constant_refrac_command(single_wl_mp):
    """Single wl must use 'refrac nr ni' (MOPSMAP interpolate_linear bug)."""
    paths = write_launching_file(single_wl_mp)
    content = paths["mopsmap"].read_text()
    match = re.search(r"refrac ([\d.eE+-]+) ([\d.eE+-]+)", content)
    assert match is not None
    assert float(match.group(1)) == pytest.approx(1.45, rel=1e-9)
    assert float(match.group(2)) == pytest.approx(1e-4, rel=1e-9)
    assert "refrac file" not in content


def test_multi_wavelength_uses_from_refrac_file(sphere_mp):
    """Multi-wl must use 'wavelength from_refrac_file' and 'refrac file'."""
    paths = write_launching_file(sphere_mp)
    content = paths["mopsmap"].read_text()
    assert "wavelength from_refrac_file" in content
    assert "refrac file" in content
