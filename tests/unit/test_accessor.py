"""The .mopsmap accessor on a result dataset."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap  # noqa: F401  registers the accessor

RH = [0.0, 50.0, 90.0]
WL = [0.44, 0.55]
THETA = np.linspace(0.0, 180.0, 7)


@pytest.fixture
def result() -> xr.Dataset:
    shape = (len(RH), len(WL))
    return xr.Dataset(
        {
            "kext": (("rh", "wl"), np.full(shape, 1e-4)),
            "ssa": (("rh", "wl"), np.full(shape, 0.95)),
            "phase": (
                ("rh", "wl", "theta"),
                np.ones((*shape, len(THETA))),
            ),
        },
        coords={"rh": RH, "wl": WL, "theta": THETA},
    )


class TestSmartgExport:
    def test_it_writes_a_lut(self, result, tmp_path):
        path = tmp_path / "sulphate_sol.nc"

        result.mopsmap.to_smartg(path, name="sulphate")

        assert path.exists()

    def test_the_lut_carries_the_smartg_dimension_names(
        self, result, tmp_path
    ):
        path = tmp_path / "lut.nc"
        result.mopsmap.to_smartg(path, name="sulphate")

        lut = xr.open_dataset(path)
        assert lut["ext"].dims == ("hum", "wav")
        assert lut["phase"].dims == ("hum", "wav", "stk", "theta")

    def test_the_angles_come_out_increasing(self, result, tmp_path):
        """SMART-G 2.0 raises on a decreasing theta; MOPSMAP emits one."""
        flipped = result.isel(theta=slice(None, None, -1))
        path = tmp_path / "lut.nc"

        flipped.mopsmap.to_smartg(path, name="sulphate")

        theta = xr.open_dataset(path)["theta"].values
        assert (np.diff(theta) > 0).all()

    def test_sorting_the_angles_carries_the_phase_with_them(
        self, result, tmp_path
    ):
        """Reordering the coordinate without its data mirrors the phase."""
        peaked = result.copy()
        peaked["phase"] = peaked["phase"] * xr.DataArray(
            np.linspace(10.0, 1.0, len(THETA)), dims="theta"
        )
        path = tmp_path / "lut.nc"

        peaked.isel(theta=slice(None, None, -1)).mopsmap.to_smartg(
            path, name="sulphate"
        )

        lut = xr.open_dataset(path)
        forward = lut["phase"].isel(hum=0, wav=0, stk=0, theta=0)
        backward = lut["phase"].isel(hum=0, wav=0, stk=0, theta=-1)
        assert float(forward) > float(backward)

    def test_the_humidity_dimension_can_be_named(self, result, tmp_path):
        """v1 hardcoded 'rh', which broke as soon as a user named it."""
        renamed = result.rename(rh="rh_nominal")
        path = tmp_path / "lut.nc"

        renamed.mopsmap.to_smartg(
            path, name="sulphate", humidity_dim="rh_nominal"
        )

        assert xr.open_dataset(path).sizes["hum"] == len(RH)

    def test_a_missing_humidity_dimension_is_reported(self, result, tmp_path):
        with pytest.raises(KeyError, match="humidity"):
            result.mopsmap.to_smartg(
                tmp_path / "lut.nc", name="sulphate", humidity_dim="absent"
            )

    def test_a_result_without_a_phase_function_is_reported(
        self, result, tmp_path
    ):
        with pytest.raises(KeyError, match="phase"):
            result.drop_vars("phase").mopsmap.to_smartg(
                tmp_path / "lut.nc", name="sulphate"
            )

    def test_the_values_survive_the_transposition(self, result, tmp_path):
        path = tmp_path / "lut.nc"
        marked = result.copy()
        marked["kext"][2, 1] = 7.0

        marked.mopsmap.to_smartg(path, name="sulphate")

        lut = xr.open_dataset(path)
        assert float(lut["ext"].isel(hum=2, wav=1)) == pytest.approx(7.0)


@pytest.fixture
def matrix_result(result) -> xr.Dataset:
    """A result holding the full scattering matrix, one value per element.

    MOPSMAP emits the elements in the order a1, a2, a3, a4, b1, b2, so
    marking them 1 to 6 makes the reordering readable in the assertion.
    """
    marks = np.arange(1.0, 7.0)
    matrix = np.broadcast_to(marks, (len(RH), len(WL), len(THETA), 6)).copy()
    return result.assign(
        scattering_matrix=(("rh", "wl", "theta", "element"), matrix)
    )


class TestSmartgPhaseMatrix:
    def test_the_wavelengths_come_out_in_nanometres(self, result, tmp_path):
        """MOPSMAP works in micrometres, SMART-G reads these in nm."""
        path = tmp_path / "lut.nc"
        result.mopsmap.to_smartg(path, name="sulphate")

        wav = xr.open_dataset(path)["wav"].values
        assert wav == pytest.approx([w * 1e3 for w in WL])

    def test_the_matrix_elements_land_in_smartg_order(
        self, matrix_result, tmp_path
    ):
        """F11, F21, F33, F34, F22, F44 from a1, b1, a3, -b2, a2, a4."""
        path = tmp_path / "lut.nc"
        matrix_result.mopsmap.to_smartg(path, name="sulphate")

        stk = xr.open_dataset(path)["phase"].isel(hum=0, wav=0, theta=0)
        assert stk.values == pytest.approx([1.0, 5.0, 3.0, -6.0, 2.0, 4.0])

    def test_a_scalar_phase_function_fills_only_f11_and_f22(
        self, result, tmp_path
    ):
        """SMART-G zero-fills a short table; say which terms are lost."""
        path = tmp_path / "lut.nc"
        result.mopsmap.to_smartg(path, name="sulphate")

        stk = xr.open_dataset(path)["phase"].isel(hum=0, wav=0, theta=0)
        assert stk.values == pytest.approx([1.0, 0.0, 0.0, 0.0, 1.0, 0.0])

    def test_sorting_the_angles_carries_the_matrix_with_them(
        self, matrix_result, tmp_path
    ):
        """The theta sort has to reach scattering_matrix, not just phase."""
        peaked = matrix_result.copy()
        peaked["scattering_matrix"] = peaked[
            "scattering_matrix"
        ] * xr.DataArray(np.linspace(10.0, 1.0, len(THETA)), dims="theta")
        path = tmp_path / "lut.nc"

        peaked.isel(theta=slice(None, None, -1)).mopsmap.to_smartg(
            path, name="sulphate"
        )

        phase = xr.open_dataset(path)["phase"].isel(hum=0, wav=0, stk=0)
        assert float(phase.isel(theta=0)) > float(phase.isel(theta=-1))
