"""The built-in OPAC catalogue, kappa flavour."""

from __future__ import annotations

import importlib.resources as resources

import pytest

import pymopsmap as pm
from pymopsmap.species import schema

SPECIES = [
    "inso",
    "miam",
    "micm",
    "minm",
    "mitr",
    "soot",
    "ssam",
    "sscm",
    "suso",
    "waso",
]

# Hess et al. 1998 Table 1c, as published on the GEISA OPAC page.
DRY = {
    "inso": (2.51, 0.471, 0.005, 20.0, 2.0),
    "waso": (2.24, 0.0212, 0.005, 20.0, 1.8),
    "soot": (2.00, 0.0118, 0.005, 20.0, 1.0),
    "ssam": (2.03, 0.209, 0.005, 20.0, 2.2),
    "sscm": (2.03, 1.75, 0.005, 60.0, 2.2),
    "minm": (1.95, 0.07, 0.005, 20.0, 2.6),
    "miam": (2.00, 0.39, 0.005, 20.0, 2.6),
    "micm": (2.15, 1.90, 0.005, 60.0, 2.6),
    "mitr": (2.20, 0.50, 0.02, 5.0, 2.5),
    "suso": (2.03, 0.0695, 0.005, 20.0, 1.7),
}


@pytest.fixture(scope="module")
def tree():
    path = resources.files("pymopsmap") / "data" / "opac" / "kappa.nc"
    return schema.read(path)


class TestLayout:
    def test_every_specie_ships(self, tree):
        assert sorted(tree.children) == SPECIES

    def test_growth_is_delegated_to_the_engine(self, tree):
        assert tree["waso"].attrs["growth"] == "kappa"

    def test_species_are_monomodal(self, tree):
        for name in SPECIES:
            assert list(tree[name].children) == ["only"]

    def test_no_humidity_axis(self, tree):
        """A kappa species holds dry values; the engine grows them."""
        assert "rh" not in tree["waso/only"].to_dataset().dims


class TestHessTable:
    @pytest.mark.parametrize("name", SPECIES)
    def test_the_dry_size_distribution_matches(self, tree, name):
        sigma, rm, rmin, rmax, density = DRY[name]
        ds = tree[f"{name}/only"].to_dataset()

        assert float(ds["sigma"]) == pytest.approx(sigma)
        assert float(ds["rm"]) == pytest.approx(rm)
        assert float(ds["rmin"]) == pytest.approx(rmin)
        assert float(ds["rmax"]) == pytest.approx(rmax)
        assert float(ds["density_dry"]) == pytest.approx(density)

    def test_it_agrees_with_the_mopsmap_user_guide(self, tree):
        """The desert example of the guide pins four of these species."""
        guide = {
            "waso": (0.0212, 2.24, 0.005, 20.0),
            "minm": (0.07, 1.95, 0.005, 20.0),
            "miam": (0.39, 2.00, 0.005, 20.0),
            "micm": (1.90, 2.15, 0.005, 60.0),
        }
        for name, (rm, sigma, rmin, rmax) in guide.items():
            ds = tree[f"{name}/only"].to_dataset()
            assert (
                float(ds["rm"]),
                float(ds["sigma"]),
                float(ds["rmin"]),
                float(ds["rmax"]),
            ) == pytest.approx((rm, sigma, rmin, rmax))


class TestHygroscopicity:
    def test_soluble_species_carry_a_positive_kappa(self, tree):
        for name in ("waso", "ssam", "sscm", "suso"):
            assert float(tree[f"{name}/only"].to_dataset()["kappa"]) > 0

    def test_insoluble_species_carry_a_zero_kappa(self, tree):
        """MOPSMAP requires a kappa on every mode as soon as rH is set."""
        for name in ("inso", "soot", "minm", "miam", "micm", "mitr"):
            assert float(tree[f"{name}/only"].to_dataset()["kappa"]) == 0.0


class TestRefractiveIndex:
    def test_the_grid_is_the_mopsmap_one(self, tree):
        wl = tree["waso/only"].to_dataset()["wl"].values

        assert len(wl) == 61
        assert wl[0] == pytest.approx(0.25)
        assert wl[-1] == pytest.approx(40.0)

    def test_the_imaginary_part_is_positive(self, tree):
        for name in SPECIES:
            assert (tree[f"{name}/only"].to_dataset()["n_imag"] >= 0).all()

    def test_soot_is_the_most_absorbing(self, tree):
        at_550 = {
            name: float(
                tree[f"{name}/only"].to_dataset()["n_imag"].interp(wl=0.55)
            )
            for name in SPECIES
        }

        assert max(at_550, key=at_550.get) == "soot"


class TestLoading:
    def test_it_loads_through_the_public_api(self):
        specie = pm.load(pm.OPAC.WASO)

        assert specie.name == "waso"
        assert specie.growth is schema.Growth.KAPPA
        assert specie.rh_range is None

    def test_the_engine_receives_the_humidity(self):
        point = pm.load(pm.OPAC.WASO).at(wl=[0.55], rh=80.0)

        assert point.engine_rh == 80.0
        assert point.modes[0].kappa == pytest.approx(0.249)
        assert point.modes[0].density == pytest.approx(1.8)

    def test_a_dry_request_needs_no_humidity(self):
        point = pm.load(pm.OPAC.SOOT).at(wl=[0.55])

        assert point.engine_rh is None
        assert point.modes[0].psd.rm == pytest.approx(0.0118)


class TestMixes:
    def test_the_climatologies_are_constructors(self):
        mix = pm.opac_mix("continental_average")

        assert sorted(mix.weights) == ["inso", "soot", "waso"]

    def test_the_concentrations_are_in_si(self):
        """The published table is in cm-3; the catalogue speaks m-3."""
        mix = pm.opac_mix("continental_average")

        assert mix.weights["waso"] == pytest.approx(7000.0 * 1e6)

    @pytest.mark.parametrize(
        "name,total",
        [
            ("continental_clean", 2600.15),
            ("continental_average", 15300.4),
            ("continental_polluted", 50000.6),
            ("urban", 158001.5),
            ("desert", 2300.142),
            ("maritime_clean", 1520.0032),
            ("maritime_polluted", 9000.0032),
            ("maritime_tropical", 600.0013),
            ("arctic", 6601.91),
            ("antarctic", 42.9523),
        ],
    )
    def test_each_mixture_sums_to_its_published_total(self, name, total):
        mix = pm.opac_mix(name)

        assert sum(mix.weights.values()) / 1e6 == pytest.approx(total)
