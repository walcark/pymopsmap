"""The canonical species NetCDF contract."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from pymopsmap.exceptions import SchemaError
from pymopsmap.species import schema

RH = [0.0, 50.0, 90.0]
WL = [0.44, 0.55, 0.67]


def _tabulated_mode() -> xr.Dataset:
    """A lognormal sphere whose wet state is tabulated against humidity."""
    return xr.Dataset(
        data_vars={
            "n_real": (("rh", "wl"), np.full((3, 3), 1.45)),
            "n_imag": (("rh", "wl"), np.full((3, 3), 1e-3)),
            "rm": (("rh",), np.array([0.05, 0.07, 0.12])),
            "sigma": (("rh",), np.full(3, 1.5)),
            "rmin": 0.001,
            "rmax": 40.0,
            "n": 3.2e9,
            "density_dry": 1.8,
        },
        coords={"rh": RH, "wl": WL},
        attrs={"psd_type": "lognormal", "shape_type": "sphere"},
    )


def _dry_mode() -> xr.Dataset:
    """The same mode, dry, with a hygroscopicity parameter."""
    return xr.Dataset(
        data_vars={
            "n_real": (("wl",), np.full(3, 1.45)),
            "n_imag": (("wl",), np.full(3, 1e-3)),
            "rm": 0.05,
            "sigma": 1.5,
            "rmin": 0.001,
            "rmax": 40.0,
            "n": 3.2e9,
            "density_dry": 1.8,
            "kappa": 0.25,
        },
        coords={"wl": WL},
        attrs={"psd_type": "lognormal", "shape_type": "sphere"},
    )


def _tree(growth: str, *modes: tuple[str, xr.Dataset]) -> xr.DataTree:
    tree = xr.DataTree()
    tree.attrs.update(source="TEST", version="v1")
    tree["specie"] = xr.DataTree()
    tree["specie"].attrs["growth"] = growth
    for name, ds in modes:
        tree[f"specie/{name}"] = xr.DataTree(ds)
    return tree


# ---------------------------------------------------------------------------
# Round trip
# ---------------------------------------------------------------------------


class TestRoundTrip:
    def test_written_tree_reads_back_identical(self, tmp_path):
        tree = _tree("tabulated", ("fine", _tabulated_mode()))
        path = tmp_path / "specie.nc"

        schema.write(tree, path)
        back = schema.read(path)

        assert list(back["specie"].children) == ["fine"]
        xr.testing.assert_allclose(
            back["specie/fine"].to_dataset(),
            tree["specie/fine"].to_dataset(),
        )

    def test_write_stamps_the_schema_revision(self, tmp_path):
        path = tmp_path / "specie.nc"
        schema.write(_tree("kappa", ("fine", _dry_mode())), path)

        assert schema.read(path).attrs["schema_rev"] == schema.SCHEMA_REV

    def test_multi_mode_species_round_trip(self, tmp_path):
        tree = _tree(
            "tabulated",
            ("fine", _tabulated_mode()),
            ("coarse", _tabulated_mode()),
        )
        path = tmp_path / "specie.nc"

        schema.write(tree, path)

        assert list(schema.read(path)["specie"].children) == [
            "fine",
            "coarse",
        ]


# ---------------------------------------------------------------------------
# Growth consistency
# ---------------------------------------------------------------------------


class TestGrowthValidation:
    def test_tabulated_without_humidity_axis_is_rejected(self):
        with pytest.raises(SchemaError, match="rh"):
            schema.validate(_tree("tabulated", ("fine", _dry_mode())))

    def test_kappa_with_a_humidity_axis_is_rejected(self):
        mode = _tabulated_mode()
        mode["kappa"] = 0.25

        with pytest.raises(SchemaError, match="rh"):
            schema.validate(_tree("kappa", ("fine", mode)))

    def test_kappa_without_a_kappa_variable_is_rejected(self):
        mode = _dry_mode().drop_vars("kappa")

        with pytest.raises(SchemaError, match="kappa"):
            schema.validate(_tree("kappa", ("fine", mode)))

    def test_none_with_a_humidity_axis_is_rejected(self):
        with pytest.raises(SchemaError):
            schema.validate(_tree("none", ("fine", _tabulated_mode())))

    def test_unknown_growth_value_is_rejected(self):
        with pytest.raises(SchemaError, match="growth"):
            schema.validate(_tree("magic", ("fine", _dry_mode())))

    def test_kappa_growth_rejects_modified_gamma(self):
        """MOPSMAP refuses mod_gamma combined with hygroscopic growth."""
        mode = _dry_mode().drop_vars(["rm", "sigma"])
        mode["A"], mode["B"], mode["alpha"], mode["gamma"] = 1.0, 1.0, 1.0, 1.0
        mode.attrs["psd_type"] = "mod-gamma"

        with pytest.raises(SchemaError, match="mod-gamma"):
            schema.validate(_tree("kappa", ("fine", mode)))


# ---------------------------------------------------------------------------
# Conventions
# ---------------------------------------------------------------------------


class TestConventions:
    def test_negative_imaginary_index_is_rejected(self):
        mode = _dry_mode()
        mode["n_imag"] = (("wl",), np.full(3, -1e-3))

        with pytest.raises(SchemaError, match="n_imag"):
            schema.validate(_tree("kappa", ("fine", mode)))

    def test_unknown_psd_type_is_rejected(self):
        mode = _dry_mode()
        mode.attrs["psd_type"] = "gaussian"

        with pytest.raises(SchemaError, match="psd_type"):
            schema.validate(_tree("kappa", ("fine", mode)))

    def test_missing_psd_variable_is_rejected(self):
        mode = _dry_mode().drop_vars("sigma")

        with pytest.raises(SchemaError, match="sigma"):
            schema.validate(_tree("kappa", ("fine", mode)))

    def test_the_distribution_amplitude_is_a_required_variable(self):
        """The ratio between modes is part of what defines the species."""
        assert "n" in schema.expected_psd_variables("lognormal")
        assert "A" in schema.expected_psd_variables("mod-gamma")

    def test_a_mode_without_its_amplitude_is_rejected(self):
        mode = _dry_mode().drop_vars("n")

        with pytest.raises(SchemaError, match="'n'"):
            schema.validate(_tree("kappa", ("fine", mode)))

    def test_wavelength_in_nanometres_is_converted(self, tmp_path):
        mode = _dry_mode()
        mode = mode.assign_coords(wl=[440.0, 550.0, 670.0])
        mode["wl"].attrs["units"] = "nm"
        path = tmp_path / "specie.nc"
        schema.write(_tree("kappa", ("fine", mode)), path, check_units=False)

        back = schema.read(path)

        np.testing.assert_allclose(back["specie/fine"]["wl"].values, WL)

    def test_unknown_length_unit_is_rejected(self):
        mode = _dry_mode()
        mode["rm"].attrs["units"] = "furlong"

        with pytest.raises(SchemaError, match="furlong"):
            schema.validate(_tree("kappa", ("fine", mode)))
