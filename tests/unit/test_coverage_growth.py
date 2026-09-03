"""Size-parameter coverage is checked on the radii MOPSMAP will use."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from pymopsmap.engine.coverage import clip_modes_to_coverage
from pymopsmap.models.microparams import LognormalPSD, MicroParameters, Sphere

# A coarse sea-salt-like mode, kappa 0.916. Dry it fits the coverage at every
# wavelength; at 90 percent humidity the growth factor is 2.10 and it does not.
# The real OPAC coarse mode goes to 60 um, which already leaves the coverage at
# 355 nm while dry, so it would not isolate growth as the cause.
WL = [0.355, 0.532, 1.064]


def _sscm(kappa: float | None = 0.916) -> MicroParameters:
    return MicroParameters(
        wavelength=WL,
        n_real=1.50,
        n_imag=1e-4,
        shape=Sphere(),
        psd=LognormalPSD(rm=1.75, sigma=2.03, n=1e6, rmin=0.005, rmax=50.0),
        kappa=kappa,
    )


class TestGrownRadii:
    def test_a_dry_mode_keeps_every_wavelength(self):
        _, mask = clip_modes_to_coverage([_sscm()], rh=None)

        assert mask.all()

    def test_growth_pushes_a_wavelength_out_of_coverage(self):
        """The dry radius fits; the grown one does not."""
        with pytest.warns(UserWarning, match="clipped"):
            _, mask = clip_modes_to_coverage([_sscm()], rh=90.0)

        assert not mask.all()

    def test_a_zero_kappa_mode_is_unaffected(self):
        _, dry = clip_modes_to_coverage([_sscm(kappa=0.0)], rh=None)
        _, wet = clip_modes_to_coverage([_sscm(kappa=0.0)], rh=90.0)

        np.testing.assert_array_equal(dry, wet)

    def test_the_clipped_modes_keep_their_own_radii(self):
        """Growth is predicted for the check, not applied to what is sent."""
        with pytest.warns(UserWarning):
            modes, _ = clip_modes_to_coverage([_sscm()], rh=90.0)

        assert modes[0].psd.rmax == pytest.approx(50.0)


class TestShapeProvenance:
    def test_a_result_records_the_shapes_it_came_from(self):
        from pymopsmap.engine.output_format import shape_types

        assert shape_types([_sscm()]) == ["sphere"]

    def test_a_mixed_result_records_every_shape(self):
        from pymopsmap.engine.output_format import shape_types
        from pymopsmap.models.microparams import Spheroid

        spheroid = _sscm().model_copy(
            update={"shape": Spheroid(mode="oblate", aspect_ratio=1.5)}
        )

        assert shape_types([_sscm(), spheroid]) == ["sphere", "spheroid"]

    def test_the_effective_radius_survives_a_spherical_mixture(self):
        """It was always NaN: the provenance was never recorded."""
        from pymopsmap.engine.output_format import combine

        one = xr.Dataset(
            {
                "kext": (("wl",), np.full(3, 1e-4)),
                "ksca": (("wl",), np.full(3, 9e-5)),
                "vol_dens": 3.0e-12,
                "cross_dens": 2.0e-6,
                "reff": 0.75 * 3.0e-12 / 2.0e-6,
            },
            coords={"wl": WL},
            attrs={"shape_types": ["sphere"]},
        )

        assert np.isfinite(combine([one, one], weights=[1.0, 1.0])["reff"])
