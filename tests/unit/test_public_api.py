"""The public surface, and that it cannot grow back unnoticed."""

from __future__ import annotations

import pymopsmap as pm

# Every name pymopsmap exports, and why it earns its place. Shapes and size
# distributions live under pm.shapes and pm.psd rather than at the top level:
# there are twelve of them, and they are only needed to build a custom species.
EXPECTED = {
    # Loading a species, and mixing several
    "load",
    "Specie",
    "Mix",
    "Mode",
    "CAMS",
    "OPAC",
    "opac_mix",
    # Building one by hand
    "shapes",
    "psd",
    "MicroParameters",
    # Asking for outputs
    "OutputType",
    # Errors a caller acts on
    "CoverageError",
    "DomainError",
    "MopsmapError",
    "SchemaError",
}


class TestSurface:
    def test_it_is_exactly_what_is_declared(self):
        assert set(pm.__all__) == EXPECTED

    def test_it_stays_small(self):
        assert len(pm.__all__) <= 16

    def test_every_exported_name_resolves(self):
        for name in pm.__all__:
            assert getattr(pm, name, None) is not None, name


class TestNamespaces:
    def test_shapes_are_reachable(self):
        assert pm.shapes.Sphere().type == "sphere"

    def test_size_distributions_are_reachable(self):
        assert (
            pm.psd.LognormalPSD(
                rm=0.1, sigma=1.5, n=1e6, rmin=0.001, rmax=40.0
            ).type
            == "lognormal"
        )

    def test_the_accessor_is_registered_on_import(self):
        import xarray as xr

        assert hasattr(xr.Dataset(), "mopsmap")
