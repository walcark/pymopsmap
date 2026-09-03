"""The output combination registry."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from pymopsmap.engine.outputs import COMBINE, combine

WL = np.array([0.44, 0.55, 0.67])


def _angstrom(values: np.ndarray) -> np.ndarray:
    """MOPSMAP's per-interval exponent; the first wavelength is NaN."""
    out = np.full(len(WL), np.nan)
    # A non-absorbing species gives a 0/0 here; NaN is the answer.
    with np.errstate(divide="ignore", invalid="ignore"):
        out[1:] = -np.log(values[1:] / values[:-1]) / np.log(WL[1:] / WL[:-1])
    return out


def _species(kext: float, ssa: float, spheres: bool = True) -> xr.Dataset:
    """A physically self-consistent single-species result."""
    kext_a = np.full(3, kext) * np.array([1.2, 1.0, 0.85])
    ksca_a = kext_a * ssa
    backscatter = ksca_a / 25.0
    vol_dens, cross_dens = 3.0e-12, 2.0e-6
    ds = xr.Dataset(
        {
            "kext": (("wl",), kext_a),
            "ssa": (("wl",), np.full(3, ssa)),
            "ksca": (("wl",), ksca_a),
            "g": (("wl",), np.full(3, 0.7)),
            "n": 1.0e9,
            "mass_conc": 4.0e-9,
            "vol_dens": vol_dens,
            "cross_dens": cross_dens,
            "reff": 0.75 * vol_dens / cross_dens,
            "backscatter": (("wl",), backscatter),
            "lidar_ratio": (("wl",), kext_a / backscatter),
            "depol_ratio": (("wl",), np.full(3, 0.05)),
            "ext_to_mass": (("wl",), np.full(3, 3.0)),
            "back_to_mass": (("wl",), np.full(3, 0.1)),
            "angstrom_ext": (("wl",), _angstrom(kext_a)),
            "angstrom_sca": (("wl",), _angstrom(ksca_a)),
            "angstrom_abs": (("wl",), _angstrom(kext_a - ksca_a)),
            "angstrom_back": (("wl",), _angstrom(backscatter)),
        },
        coords={"wl": WL},
    )
    ds.attrs["shape_types"] = ["sphere"] if spheres else ["spheroid"]
    return ds


class TestRegistryCoverage:
    def test_every_parsed_variable_has_a_rule(self):
        parsed = set(_species(1e-4, 0.9).data_vars)

        assert parsed <= set(COMBINE)

    def test_an_unregistered_variable_is_refused(self):
        one = _species(1e-4, 0.9).assign(mystery=("wl", np.ones(3)))

        with pytest.raises(KeyError, match="mystery"):
            combine([one, one], weights=[1.0, 1.0])


class TestIdentity:
    """A mixture of one species must reproduce that species exactly."""

    @pytest.mark.parametrize("name", sorted(_species(1e-4, 0.9).data_vars))
    def test_every_variable_round_trips(self, name):
        one = _species(1e-4, 0.9)

        out = combine([one], weights=[1.0])

        np.testing.assert_allclose(
            out[name].values, one[name].values, rtol=1e-12, err_msg=name
        )


class TestTwoSpecies:
    @pytest.fixture
    def mixed(self):
        return combine(
            [_species(1e-4, 1.0), _species(1e-5, 0.2)], weights=[1.0, 1.0]
        )

    def test_extinction_is_additive(self, mixed):
        expected = _species(1e-4, 1.0)["kext"] + _species(1e-5, 0.2)["kext"]

        np.testing.assert_allclose(mixed["kext"].values, expected.values)

    def test_albedo_is_a_ratio_of_sums_not_a_mean(self, mixed):
        a, b = _species(1e-4, 1.0), _species(1e-5, 0.2)
        expected = (a["ksca"] + b["ksca"]) / (a["kext"] + b["kext"])

        np.testing.assert_allclose(mixed["ssa"].values, expected.values)
        assert not np.allclose(mixed["ssa"].values, 0.6)

    def test_asymmetry_is_weighted_by_scattering(self, mixed):
        a, b = _species(1e-4, 1.0), _species(1e-5, 0.2)
        expected = (a["ksca"] * a["g"] + b["ksca"] * b["g"]) / (
            a["ksca"] + b["ksca"]
        )

        np.testing.assert_allclose(mixed["g"].values, expected.values)

    def test_angstrom_is_recomputed_from_the_combined_spectrum(self, mixed):
        expected = _angstrom(mixed["kext"].values)

        np.testing.assert_allclose(
            mixed["angstrom_ext"].values, expected, equal_nan=True
        )

    def test_lidar_ratio_is_the_ratio_of_the_sums(self, mixed):
        np.testing.assert_allclose(
            mixed["lidar_ratio"].values,
            (mixed["kext"] / mixed["backscatter"]).values,
        )


class TestSweptDimensions:
    """A swept humidity puts the wavelength axis last, not first."""

    @pytest.fixture
    def swept(self):
        one = _species(1e-4, 0.9)
        return xr.concat(
            [one, _species(2e-4, 0.8)],
            dim=xr.DataArray([0.0, 90.0], dims="rh", name="rh"),
        )

    def test_it_combines_without_broadcasting_on_the_wrong_axis(self, swept):
        out = combine([swept, swept], weights=[1.0, 1.0])

        assert out["kext"].dims == ("rh", "wl")

    def test_the_angstrom_follows_the_wavelength_axis(self, swept):
        out = combine([swept, swept], weights=[1.0, 1.0])

        expected = _angstrom(out["kext"].sel(rh=90.0).values)
        np.testing.assert_allclose(
            out["angstrom_ext"].sel(rh=90.0).values, expected, equal_nan=True
        )

    def test_each_swept_point_gets_its_own_exponent(self, swept):
        out = combine([swept, swept], weights=[1.0, 1.0])

        assert out["angstrom_ext"].dims == ("rh", "wl")
        assert np.isnan(out["angstrom_ext"].isel(wl=0)).all()


class TestWeights:
    def test_a_weight_scales_the_additive_variables(self):
        one = _species(1e-4, 0.9)

        out = combine([one], weights=[3.0])

        np.testing.assert_allclose(
            out["kext"].values, 3.0 * one["kext"].values
        )

    def test_a_weight_leaves_the_intensive_variables_alone(self):
        one = _species(1e-4, 0.9)

        out = combine([one], weights=[3.0])

        np.testing.assert_allclose(out["ssa"].values, one["ssa"].values)
        np.testing.assert_allclose(out["g"].values, one["g"].values)


class TestEffectiveRadius:
    def test_it_is_recomputed_from_the_summed_moments(self):
        a, b = _species(1e-4, 1.0), _species(1e-5, 0.2)

        out = combine([a, b], weights=[1.0, 1.0])

        expected = (
            0.75
            * (a["vol_dens"] + b["vol_dens"])
            / (a["cross_dens"] + b["cross_dens"])
        )
        np.testing.assert_allclose(out["reff"].values, expected.values)

    def test_it_is_not_available_for_non_spherical_modes(self):
        """cross_dens and vol_dens then carry different powers of the
        aspect ratio, so the moments cannot be recovered."""
        with pytest.warns(UserWarning, match="reff"):
            out = combine(
                [_species(1e-4, 1.0), _species(1e-5, 0.2, spheres=False)],
                weights=[1.0, 1.0],
            )

        assert np.isnan(out["reff"].values).all()
