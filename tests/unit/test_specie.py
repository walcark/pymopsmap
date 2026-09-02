"""Loading a species and materialising one point of it."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import pymopsmap as pm
from pymopsmap.exceptions import DomainError, SchemaError
from pymopsmap.species.schema import Growth

WL = [0.44, 0.55, 0.67]


@pytest.fixture(scope="module")
def sulphate():
    return pm.load(pm.CAMS.SULPHATE)


def _dry_specie(growth: str, tmp_path, **extra) -> object:
    """A hand-built single-mode species, written and reloaded."""
    ds = xr.Dataset(
        data_vars={
            "n_real": (("wl",), np.full(3, 1.45)),
            "n_imag": (("wl",), np.full(3, 1e-3)),
            "rm": 0.05,
            "sigma": 1.5,
            "rmin": 0.001,
            "rmax": 40.0,
            "n": 1e9,
            "density_dry": 1.8,
            **extra,
        },
        coords={"wl": WL},
        attrs={"psd_type": "lognormal", "shape_type": "sphere"},
    )
    tree = xr.DataTree()
    tree.attrs.update(source="TEST", version="v1")
    tree["aer"] = xr.DataTree()
    tree["aer"].attrs["growth"] = growth
    tree["aer/only"] = xr.DataTree(ds)

    path = tmp_path / f"{growth}.nc"
    pm.species.schema.write(tree, path)
    return pm.load(path)


class TestLoad:
    def test_loads_from_the_builtin_catalogue(self, sulphate):
        assert sulphate.name == "sulphate"
        assert sulphate.source == "CAMS"
        assert sulphate.growth is Growth.TABULATED
        assert sulphate.modes == ["fine", "coarse"]

    def test_defaults_to_the_latest_version(self, sulphate):
        assert sulphate.version == "49r1"

    def test_an_explicit_version_is_honoured(self):
        assert pm.load(pm.CAMS.SULPHATE, version="47r1").version == "47r1"

    def test_a_specie_absent_from_a_version_is_reported(self):
        with pytest.raises(KeyError, match="47r1"):
            pm.load(pm.CAMS.SECONDARY_ORGANIC, version="47r1")

    def test_an_unknown_version_is_reported(self):
        with pytest.raises(ValueError, match="99r9"):
            pm.load(pm.CAMS.SULPHATE, version="99r9")

    def test_loads_from_a_path(self, tmp_path):
        assert _dry_specie("none", tmp_path).name == "aer"

    def test_a_file_holding_several_species_needs_a_name(self):
        path = pm.species.catalog.path_for(pm.CAMS.SULPHATE, "49r1")

        with pytest.raises(ValueError, match="several"):
            pm.load(path)


class TestRanges:
    def test_humidity_range_comes_from_the_table(self, sulphate):
        assert sulphate.rh_range == (0.0, 95.0)

    def test_wavelength_range_comes_from_the_table(self, sulphate):
        assert sulphate.wl_range == pytest.approx((0.2, 2.35))

    def test_a_dry_specie_has_no_humidity_range(self, tmp_path):
        assert _dry_specie("none", tmp_path).rh_range is None


class TestMaterialisation:
    def test_returns_one_micro_parameters_per_mode(self, sulphate):
        point = sulphate.at(wl=WL, rh=50.0)

        assert len(point.modes) == 2
        assert [m.psd.type for m in point.modes] == ["lognormal"] * 2
        assert [m.shape.type for m in point.modes] == ["sphere"] * 2

    def test_values_match_the_interpolated_table(self, sulphate):
        point = sulphate.at(wl=WL, rh=50.0)

        expected = sulphate.tree["fine"].to_dataset().interp(rh=50.0, wl=WL)
        fine = point.modes[0]
        assert fine.psd.rm == pytest.approx(float(expected["rm"]))
        assert fine.psd.n == pytest.approx(float(expected["n"]))
        np.testing.assert_allclose(fine.n_real, expected["n_real"].values)
        np.testing.assert_allclose(fine.n_imag, expected["n_imag"].values)

    def test_truncation_radii_come_from_the_file(self, sulphate):
        fine = sulphate.at(wl=WL, rh=50.0).modes[0]

        assert fine.psd.rmin == pytest.approx(0.001)
        assert fine.psd.rmax == pytest.approx(40.0)

    def test_a_tabulated_specie_is_already_wet(self, sulphate):
        """The engine must not grow it a second time."""
        assert sulphate.at(wl=WL, rh=50.0).engine_rh is None

    def test_wavelength_outside_the_table_raises(self, sulphate):
        with pytest.raises(DomainError, match="wl"):
            sulphate.at(wl=[2.5], rh=50.0)

    def test_humidity_outside_the_table_raises(self, sulphate):
        with pytest.raises(DomainError, match="rh"):
            sulphate.at(wl=WL, rh=99.0)

    def test_a_tabulated_specie_requires_a_humidity(self, sulphate):
        with pytest.raises(ValueError, match="rh"):
            sulphate.at(wl=WL)


class TestGrowthModes:
    def test_a_dry_specie_refuses_a_humidity(self, tmp_path):
        specie = _dry_specie("none", tmp_path)

        with pytest.raises(ValueError, match="does not model humidity"):
            specie.at(wl=WL, rh=50.0)

    def test_a_dry_specie_computes_without_one(self, tmp_path):
        point = _dry_specie("none", tmp_path).at(wl=WL)

        assert point.modes[0].psd.rm == pytest.approx(0.05)
        assert point.engine_rh is None

    def test_a_kappa_specie_hands_growth_to_the_engine(self, tmp_path):
        specie = _dry_specie("kappa", tmp_path, kappa=0.25)

        point = specie.at(wl=WL, rh=80.0)

        assert point.engine_rh == 80.0
        assert point.modes[0].kappa == pytest.approx(0.25)
        assert point.modes[0].density == pytest.approx(1.8)

    def test_a_kappa_specie_accepts_an_override(self, tmp_path):
        specie = _dry_specie("kappa", tmp_path, kappa=0.25)

        point = specie.at(wl=WL, rh=80.0, kappa=0.6)

        assert point.modes[0].kappa == pytest.approx(0.6)

    def test_the_engine_humidity_ceiling_is_enforced(self, tmp_path):
        specie = _dry_specie("kappa", tmp_path, kappa=0.25)

        with pytest.raises(DomainError, match="rh"):
            specie.at(wl=WL, rh=99.5)


class TestSchemaStillGuards:
    def test_a_file_failing_validation_is_refused(self, tmp_path):
        path = tmp_path / "bad.nc"
        tree = xr.DataTree()
        tree["aer"] = xr.DataTree()
        tree["aer"].attrs["growth"] = "tabulated"
        tree["aer/only"] = xr.DataTree(
            xr.Dataset(attrs={"psd_type": "lognormal", "shape_type": "sphere"})
        )
        tree.to_netcdf(path)

        with pytest.raises(SchemaError):
            pm.load(path)
