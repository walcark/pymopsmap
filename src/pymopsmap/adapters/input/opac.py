"""OPAC aerosol mix adapter — builds MicroParameters from OPAC .nc files."""

from dataclasses import dataclass, field
from enum import StrEnum

import xarray as xr

from pymopsmap.models import (
    LognormalPSD,
    MicroParameters,
    OptiProps,
    ParametricSweep,
    ParticleMixture,
    Sphere,
)
from pymopsmap.models.output_request import DEFAULT_OUTPUT, OutputRequest
from pymopsmap.utils import DATA_PATH

from .opac_download import download_opac_microparams

# Number concentrations [cm⁻³] for standard OPAC mixes — Hess et al. 1998, Table 3
OPAC_MIXES: dict[str, dict[str, float]] = {
    "continental_clean": {"waso": 2600.0, "inso": 0.15},
    "continental_average": {"waso": 7000.0, "inso": 0.4, "soot": 8300.0},
    "continental_polluted": {"waso": 15700.0, "inso": 0.6, "soot": 34300.0},
    "urban": {"waso": 15000.0, "inso": 1.5, "soot": 130000.0},
    "maritime_clean": {"waso": 1500.0, "ssam": 20.0, "sscm": 0.1},
    "maritime_polluted": {"waso": 15000.0, "ssam": 20.0, "sscm": 0.0032},
    "maritime_tropical": {"waso": 590.0, "ssam": 10.0, "sscm": 0.0013},
    "desert": {"waso": 2000.0, "minm": 269.5, "miam": 30.5, "micm": 0.142},
    "arctic": {"inso": 0.01, "waso": 1300.0, "soot": 5300.0, "ssam": 1.90},
    "antarctic": {"ssam": 0.047, "mitr": 0.0053, "suso": 42.9},
}


class OpacMixName(StrEnum):
    CONTINENTAL_CLEAN = "continental_clean"
    CONTINENTAL_AVERAGE = "continental_average"
    CONTINENTAL_POLLUTED = "continental_polluted"
    URBAN = "urban"
    MARITIME_CLEAN = "maritime_clean"
    MARITIME_POLLUTED = "maritime_polluted"
    MARITIME_TROPICAL = "maritime_tropical"
    DESERT = "desert"
    ARCTIC = "arctic"
    ANTARCTIC = "antarctic"


@dataclass
class OpacMix:
    """An OPAC aerosol mix defined by species and their number concentrations."""

    components: dict[str, float] | OpacMixName = field(default_factory=dict)

    def __post_init__(self):
        if isinstance(self.components, OpacMixName):
            self.components = OPAC_MIXES[self.components]

    def compute(
        self,
        wavelengths: list[float],
        rhs: list[float],
        output_types: OutputRequest = DEFAULT_OUTPUT,
        quiet: bool = False,
    ) -> OptiProps:
        """Compute optical properties for this mix over a wavelength × RH grid."""
        from pymopsmap.engine import compute_optical_properties

        sweep = self._build_sweep(wavelengths, rhs)
        return compute_optical_properties(sweep, output_types=output_types, quiet=quiet)

    def _build_sweep(
        self, wavelengths: list[float], rhs: list[float]
    ) -> ParametricSweep:
        datasets = {
            sp: xr.open_dataset(DATA_PATH / f"opac/{sp}.nc") for sp in self.components
        }
        sweep = ParametricSweep()
        for rh in rhs:
            modes = [
                _build_mode(sp, n, datasets[sp], wavelengths, rh)
                for sp, n in self.components.items()
            ]
            sweep.add(ParticleMixture(modes), {"rh": rh})
        return sweep


def _build_mode(
    species: str,
    n_cm3: float,
    ds: xr.Dataset,
    wavelengths: list[float],
    rh: float,
) -> MicroParameters:
    sel = ds.interp(rh=rh)
    refr = sel.interp(wavelength=wavelengths)
    return MicroParameters(
        wavelength=wavelengths,
        n_real=refr["n_real"].values.tolist(),
        n_imag=(-refr["n_imag"].values).tolist(),
        shape=Sphere(),
        psd=LognormalPSD(
            rm=float(sel["rmod"]),
            sigma=float(sel["sigma"]),
            n=n_cm3 * 1e6,  # cm⁻³ → m⁻³
            rmin=float(sel["rmin"]),
            rmax=float(sel["rmax"]),
        ),
    )


__all__ = [
    "OPAC_MIXES",
    "OpacMixName",
    "OpacMix",
    "download_opac_microparams",
]
