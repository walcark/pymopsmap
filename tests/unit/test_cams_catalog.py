"""The built-in CAMS catalogue, and its agreement with the legacy adapter."""

from __future__ import annotations

import importlib.resources as resources
import json

import numpy as np
import pytest

from pymopsmap.species import schema

VERSIONS = ["47r1", "48r1", "49r1"]
SPECIES = [
    "ammonium",
    "black_carbon",
    "continen",
    "dust",
    "nitrate",
    "organic_matter",
    "sea_salt",
    "secondary_organic",
    "sulphate",
]


def _catalogue(version: str):
    path = resources.files("pymopsmap") / "data" / "cams" / f"{version}.nc"
    return schema.read(path)


@pytest.fixture(scope="module")
def tree():
    return _catalogue("49r1")


class TestCatalogueLayout:
    @pytest.mark.parametrize("version", VERSIONS)
    def test_every_version_ships_and_validates(self, version):
        tree = _catalogue(version)

        assert tree.attrs["source"] == "CAMS"
        assert tree.attrs["version"] == version

    def test_a_specie_absent_from_a_version_is_not_written(self):
        """secondary_organic did not exist in 47r1."""
        assert "secondary_organic" not in _catalogue("47r1").children
        assert "secondary_organic" in _catalogue("48r1").children

    def test_no_stored_value_is_missing(self):
        """One file per version removes the NaN padding of the source table."""
        for version in VERSIONS:
            tree = _catalogue(version)
            for node in tree.subtree:
                if not node.has_data:
                    continue
                for name, var in node.to_dataset().variables.items():
                    assert not var.isnull().any(), f"{node.path}/{name}"

    def test_every_specie_is_present_with_both_modes(self, tree):
        assert sorted(tree.children) == SPECIES  # 49r1 carries them all
        for name in SPECIES:
            assert list(tree[name].children) == ["fine", "coarse"]

    def test_growth_is_tabulated(self, tree):
        assert tree["sulphate"].attrs["growth"] == "tabulated"

    def test_discriminators_match_the_pydantic_literals(self, tree):
        mode = tree["sulphate/fine"]
        assert mode.attrs["psd_type"] == "lognormal"
        assert mode.attrs["shape_type"] == "sphere"


class TestConventions:
    def test_wavelength_is_in_micrometres(self, tree):
        wl = tree["sulphate/fine"]["wl"].values

        assert wl.min() == pytest.approx(0.2)
        assert wl.max() == pytest.approx(2.35)

    def test_width_is_sigma_not_its_logarithm(self, tree):
        sigma = tree["sulphate/fine"]["sigma"].values

        assert (sigma > 1.0).all()

    def test_imaginary_index_is_positive(self, tree):
        for name in SPECIES:
            for mode in ("fine", "coarse"):
                assert (tree[f"{name}/{mode}"]["n_imag"] >= 0).all()

    def test_truncation_radii_are_stored(self, tree):
        mode = tree["sulphate/fine"]

        assert float(mode["rmin"]) == pytest.approx(0.001)
        assert float(mode["rmax"]) == pytest.approx(40.0)

    def test_amplitudes_carry_the_modal_ratio(self, tree):
        raw = json.loads(
            (
                resources.files("pymopsmap") / "data" / "cams" / "modes.json"
            ).read_text()
        )
        for name in SPECIES:
            fine = float(tree[f"{name}/fine"]["n"])
            coarse = float(tree[f"{name}/coarse"]["n"])
            assert [fine, coarse] == pytest.approx(raw[name])


class TestAgreementWithTheLegacyAdapter:
    """The catalogue must carry the same physics as the adapter it replaces."""

    @pytest.mark.parametrize("specie", ["sulphate", "sea_salt", "dust"])
    @pytest.mark.parametrize("rh", [0.0, 50.0, 90.0])
    def test_same_size_distribution_and_refractive_index(self, specie, rh):
        from pymopsmap.adapters.input.cams import (
            CamsAerosol,
            CamsVersion,
            read_aerosol_microphysical_parameters,
        )

        wl = [0.44, 0.55, 0.67]
        legacy = read_aerosol_microphysical_parameters(
            aerosol=CamsAerosol(specie),
            version=CamsVersion.V49_R1,
            wl_microns=wl,
            rh=[rh],
        ).mixtures[0]

        tree = _catalogue("49r1")
        for legacy_mode, name in zip(legacy.modes, ("fine", "coarse")):
            node = tree[f"{specie}/{name}"].to_dataset().interp(rh=rh, wl=wl)

            assert float(node["rm"]) == pytest.approx(
                legacy_mode.psd.rm, abs=1e-6
            )
            assert float(node["sigma"]) == pytest.approx(
                legacy_mode.psd.sigma, abs=1e-6
            )
            assert float(node["n"]) == pytest.approx(legacy_mode.psd.n)
            np.testing.assert_allclose(
                node["n_real"].values, legacy_mode.n_real, rtol=1e-12
            )
            np.testing.assert_allclose(
                node["n_imag"].values, legacy_mode.n_imag, rtol=1e-12
            )
