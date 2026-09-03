"""The launch file asks only for outputs the pipeline reads back."""

from __future__ import annotations

import pytest

from pymopsmap.engine.launch_file import write_launching_file
from pymopsmap.engine.outputs import OutputType
from pymopsmap.microparams import MicroParameters
from pymopsmap.psd import FixedPSD
from pymopsmap.shapes import Sphere


def _mp() -> MicroParameters:
    return MicroParameters(
        wavelength=[0.55],
        n_real=1.45,
        n_imag=1e-4,
        shape=Sphere(),
        psd=FixedPSD(radius=0.5, n=1e6),
    )


class TestNoUnreadOutput:
    def test_the_netcdf_output_is_not_requested(self):
        """
        Nothing parses it: the results come from stdout and the ascii files.
        Asking for it segfaults a binary built against another netcdf-fortran,
        after the computation has already succeeded.
        """
        paths = write_launching_file(mp=_mp())

        assert "output netcdf" not in paths["mopsmap"].read_text()

    def test_no_netcdf_path_is_advertised(self):
        assert "netcdf" not in write_launching_file(mp=_mp())

    def test_the_integrated_output_is_still_requested(self):
        content = write_launching_file(mp=_mp())["mopsmap"].read_text()

        assert "output integrated" in content

    def test_ascii_outputs_are_still_requested(self):
        paths = write_launching_file(
            mp=_mp(),
            output_types=frozenset(
                {OutputType.INTEGRATED, OutputType.PHASE_FUNCTION}
            ),
        )
        content = paths["mopsmap"].read_text()

        assert "output ascii_file" in content
        assert "output phase_function" in content


class TestLibraryPath:
    def test_the_environment_carries_the_fortran_netcdf_library(self):
        """
        The binary is linked against libnetcdff, which lives in the
        environment rather than on the system.
        """
        from pymopsmap.engine.launcher import _runtime_env

        env = _runtime_env()

        assert "LD_LIBRARY_PATH" in env

    def test_an_existing_library_path_is_preserved(self, monkeypatch):
        from pymopsmap.engine.launcher import _runtime_env

        monkeypatch.setenv("LD_LIBRARY_PATH", "/somewhere/else")

        assert "/somewhere/else" in _runtime_env()["LD_LIBRARY_PATH"]


@pytest.mark.skipif(
    not __import__("pathlib").Path("bin/mopsmap/mopsmap").exists(),
    reason="MOPSMAP binary not present",
)
class TestTheBinaryRuns:
    def test_it_loads_its_shared_libraries(self):
        """A missing libnetcdff shows up as exit code 127."""
        import subprocess

        from pymopsmap.engine.launcher import _runtime_env
        from pymopsmap.utils import MOPSMAP_PATH

        proc = subprocess.run(
            [str(MOPSMAP_PATH)],
            capture_output=True,
            text=True,
            env=_runtime_env(),
        )

        assert "cannot open shared object" not in proc.stderr
