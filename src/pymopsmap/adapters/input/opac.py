"""OPAC aerosol mix adapter — builds MicroParameters from OPAC .nc files."""

from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np
import xarray as xr
from scipy.optimize import curve_fit

from pymopsmap.models import (
    LognormalPSD,
    MicroParameters,
    OptiProps,
    ParametricSweep,
    ParticleMixture,
    Sphere,
)
from pymopsmap.models.output_request import DEFAULT_OUTPUT, OutputRequest
from pymopsmap.utils import DATA_PATH, MOPSMAP_PATH, check_within_grid

from .opac_download import download_opac_microparams

# Number concentrations [cm⁻³] for  OPAC mixes — Hess et al. 1998, Table 3
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

# κ values from Zieger et al. 2013 (ACP 13:10609), Table 4 modified g(RH).
# Non-hygroscopic species (inso, soot, mineral modes) have κ = 0.
OPAC_KAPPA_ZIEGER2013: dict[str, float] = {
    "waso": 0.249,
    "ssam": 0.916,
    "sscm": 0.916,
    "suso": 0.804,
    "inso": 0.0,
    "soot": 0.0,
    "minm": 0.0,
    "miam": 0.0,
    "micm": 0.0,
    "mitr": 0.0,
}

# Dry material densities [g/cm³] — Hess et al. 1998
OPAC_DRY_DENSITY: dict[str, float] = {
    "inso": 2.0,
    "waso": 1.8,
    "soot": 1.0,
    "ssam": 2.2,
    "sscm": 2.2,
    "suso": 1.7,
    "minm": 2.6,
    "miam": 2.6,
    "micm": 2.6,
    "mitr": 2.6,
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


class OpacHumidityMode(StrEnum):
    GEISA = "geisa"
    KAPPA = "kappa"


@dataclass
class OpacMix:
    """An OPAC aerosol mix defined by species and number concentrations."""

    components: dict[str, float] | OpacMixName = field(default_factory=dict)

    def __post_init__(self):
        if isinstance(self.components, OpacMixName):
            self.components = OPAC_MIXES[self.components]

    def compute(
        self,
        wavelengths: list[float],
        rhs: list[float],
        output_types: OutputRequest = DEFAULT_OUTPUT,
        mode: OpacHumidityMode = OpacHumidityMode.GEISA,
        quiet: bool = False,
    ) -> OptiProps:
        """
        Compute optical properties for this mix over a wavelength × RH grid.

        Parameters
        ----------
        mode : OpacHumidityMode
            GEISA — wet PSD and refractive indices interpolated directly from
            GEISA tables (more faithful to OPAC data).
            KAPPA — hygroscopic growth from the κ parameterisation fitted to
            the GEISA data; refractive indices mixed with water using the
            volume-weighted rule (reproduces Gasteiger & Wiegner 2018 Fig. 5).
        """
        from pymopsmap.engine import compute_optical_properties

        datasets = {}
        for sp in self.components:
            datasets[sp] = xr.open_dataset(DATA_PATH / f"opac/{sp}.nc")

        for sp, ds in datasets.items():
            where = f"OPAC {sp}"
            check_within_grid(
                "wavelength", wavelengths, ds["wavelength"], where
            )
            if mode == OpacHumidityMode.GEISA:
                check_within_grid("rh", rhs, ds["rh"], where)

        sweep = ParametricSweep()
        for rh in rhs:
            if mode == OpacHumidityMode.GEISA:
                modes = [
                    _build_mode_geisa(sp, n, datasets[sp], wavelengths, rh)
                    for sp, n in self.components.items()  # type: ignore[union-attr]
                ]
            else:
                modes = [
                    _build_mode_kappa(sp, n, datasets[sp], wavelengths, rh)
                    for sp, n in self.components.items()  # type: ignore[union-attr]
                ]
            sweep.add(ParticleMixture(modes), {"rh": rh})
        return compute_optical_properties(
            sweep,
            output_types=output_types,
            quiet=quiet,
        )


# ---------------------------------------------------------------------------
# GEISA mode
# ---------------------------------------------------------------------------


def _build_mode_geisa(
    species: str,
    n_cm3: float,
    ds: xr.Dataset,
    wavelengths: list[float],
    rh: float,
) -> MicroParameters:
    """Wet MicroParameters from RH-interpolated GEISA tables."""
    sel = ds.interp(rh=rh)
    refr = sel.interp(wavelength=wavelengths)

    rmod_dry = float(ds.sel(rh=0, method="nearest")["rmod"])
    gf3 = (float(sel["rmod"]) / rmod_dry) ** 3
    rho_dry = OPAC_DRY_DENSITY[species]
    rho_wet = (rho_dry + 1.0 * (gf3 - 1.0)) / gf3

    return MicroParameters(
        wavelength=wavelengths,
        n_real=refr["n_real"].values.tolist(),
        n_imag=(-refr["n_imag"].values).tolist(),
        shape=Sphere(),
        psd=LognormalPSD(
            rm=float(sel["rmod"]),
            sigma=float(sel["sigma"]),
            n=n_cm3 * 1e6,
            rmin=float(sel["rmin"]),
            rmax=float(sel["rmax"]),
        ),
        density=rho_wet,
    )


# ---------------------------------------------------------------------------
# KAPPA mode
# ---------------------------------------------------------------------------


def _fit_kappa(ds: xr.Dataset) -> float:
    """Fit κ from GEISA rmod vs RH (radius growth factor model)."""
    rh_vals = ds.rh.values.astype(float)
    rmod_vals = ds["rmod"].values.astype(float)

    if rmod_vals.std() < 1e-12:
        return 0.0  # non-hygroscopic

    gf = rmod_vals / rmod_vals[0]
    mask = rh_vals > 0
    rh_fit, gf_fit = rh_vals[mask], gf[mask]

    def model(rh, kappa):
        return (1.0 + kappa * rh / (100.0 - rh)) ** (1.0 / 3.0)

    popt, _ = curve_fit(model, rh_fit, gf_fit, p0=[0.5], bounds=(0.0, 5.0))
    return float(popt[0])


def _water_refrac(wavelengths: list[float]) -> tuple[list[float], list[float]]:
    """
    Interpolate Segelstein water refractive index to target wavelengths [µm].
    Returns (n_real, n_imag) with n_imag > 0 (abs coefficient convention).
    """
    segelstein = MOPSMAP_PATH.parent / "data" / "refr_water_segelstein"
    data = np.loadtxt(segelstein, comments="#")
    nr = np.interp(wavelengths, data[:, 0], data[:, 1])
    ni = np.interp(wavelengths, data[:, 0], data[:, 2])
    return nr.tolist(), ni.tolist()


def _build_mode_kappa(
    species: str,
    n_cm3: float,
    ds: xr.Dataset,
    wavelengths: list[float],
    rh: float,
) -> MicroParameters:
    """
    Wet MicroParameters via manual κ parameterisation (Zieger et al. 2013).

    Applies the κ growth formula externally and mixes the dry refractive index
    with water (volume-weighted rule), so that MOPSMAP receives standard wet
    parameters without needing kappa/rH commands.
    """
    kappa = OPAC_KAPPA_ZIEGER2013[species]

    gf3 = 1.0 + kappa * rh / (100.0 - rh) if rh > 0 else 1.0
    gf = gf3 ** (1.0 / 3.0)
    f_w = (gf3 - 1.0) / gf3  # water volume fraction

    sel_dry = ds.sel(rh=0, method="nearest")
    refr_dry = sel_dry.interp(wavelength=wavelengths)

    nr_dry = refr_dry["n_real"].values.tolist()
    ni_dry = (-refr_dry["n_imag"].values).tolist()  # positive convention

    nr_water, ni_water = _water_refrac(wavelengths)
    nr_wet = [(1.0 - f_w) * nr + f_w * nw for nr, nw in zip(nr_dry, nr_water)]
    ni_wet = [(1.0 - f_w) * ni + f_w * nw for ni, nw in zip(ni_dry, ni_water)]

    rho_dry = OPAC_DRY_DENSITY[species]
    rho_wet = (1.0 - f_w) * rho_dry + f_w * 1.0

    return MicroParameters(
        wavelength=wavelengths,
        n_real=nr_wet,
        n_imag=ni_wet,
        shape=Sphere(),
        psd=LognormalPSD(
            rm=float(sel_dry["rmod"]) * gf,
            sigma=float(sel_dry["sigma"]),
            n=n_cm3 * 1e6,
            rmin=float(sel_dry["rmin"]) * gf,
            rmax=float(sel_dry["rmax"]),
        ),
        density=rho_wet,
    )


__all__ = [
    "OPAC_MIXES",
    "OPAC_DRY_DENSITY",
    "OPAC_KAPPA_ZIEGER2013",
    "OpacMixName",
    "OpacHumidityMode",
    "OpacMix",
    "download_opac_microparams",
]
