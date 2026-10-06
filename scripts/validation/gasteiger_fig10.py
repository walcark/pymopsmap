"""
Reproduce Figure 10 of Gasteiger and Wiegner (2018), GMD 11:2739.

    PYMOPSMAP_DATASET_SOURCE=... \
        pixi run -e dev python -m scripts.validation.gasteiger_fig10

"Modeled correction factors Cts for total scattering (a) and Cbs for
hemispheric backscattering (b) of an Aurora 3000 nephelometer as a function of
particle size."

A nephelometer does not see the whole sphere: it collects between about 10 and
171 degrees, with a sensitivity of its own. The correction factor is the ratio
between the coefficient a perfect instrument would measure and the one this
geometry gives, normalised on the value for the smallest particles, which
scatter like Rayleigh.

Every input is in the authors' own script,
``bin/mopsmap/misc/paper_examples/sect_57_neph_truncation/calc_trunc_corr.py``:
the angular sensitivity of Mueller et al. (2011), the OPAC mineral refractive
index, a non-absorbing fraction of 0.5, and a lognormal mode of width 1.6 cut
at 5 um, which is a PM10 inlet.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap import MicroParameters
from pymopsmap.engine import run_point
from pymopsmap.engine.outputs import OutputType
from pymopsmap.psd import LognormalPSD
from pymopsmap.shapes import Sphere, SpheroidDistrFile

FIGURES = Path(__file__).resolve().parents[2] / "docs" / "figures"
AR_KANDLER = "ar_kandler"

# The three operating wavelengths of the instrument, and the OPAC mineral
# index at each, as the shipped refr_mineral of the example gives it.
WAVELENGTHS_UM = [0.45, 0.525, 0.635]
N_REAL = [1.53, 1.53, 1.53]
N_IMAG = [0.00830, 0.00665, 0.00450]

NONABS_FRACTION = 0.5
SIGMA = 1.6
R_MIN_UM = 0.001
R_MAX_UM = 5.0
N_ANGLES = 1801

OUTPUTS = frozenset({OutputType.PHASE_FUNCTION})
COLOURS = ["b", "g", "r"]


def mode_radii(count: int = 150) -> np.ndarray:
    """The mode radii of the reference script: 0.001 to 2 um, log spaced."""
    return np.logspace(np.log10(0.001), np.log10(2.0), count)


def weights(theta_deg: np.ndarray) -> dict[str, np.ndarray]:
    """
    The four angular weights the correction factors are ratios of.

    The ideal instrument weights by sin(theta) over the whole sphere, or over
    the backward hemisphere. The real one sees 10 to 171 degrees with the
    sensitivity Mueller et al. (2011) measured, and its backscatter channel
    opens gradually between 70.25 and 110.24 degrees.
    """
    theta = np.radians(theta_deg)
    # sin(pi) comes out a hair below zero, and a fractional power of that is
    # not a number.
    sine = np.clip(np.sin(theta), 0.0, None)
    sensitivity = 1.01 * sine**1.19
    inside = (theta_deg > 10.0) & (theta_deg < 171.0)
    opening = np.clip((theta_deg - 70.25) / 39.99, 0.0, 1.0)
    return {
        "total_ideal": sine,
        "total_real": np.where(inside, sensitivity, 0.0),
        "back_ideal": np.where(theta_deg >= 90.0, sine, 0.0),
        "back_real": np.where(inside, sensitivity * opening, 0.0),
    }


def compute() -> xr.Dataset:
    """
    The two correction factors, over the mode radii, for both shapes.

    Returns
    -------
    xr.Dataset
        Dimensions ``(shape, rmod, wl)``, variables ``cts`` and ``cbs``,
        each already normalised on its Rayleigh value.
    """
    radii = mode_radii()
    per_shape = []
    for spherical in (True, False):
        raw = np.full((len(radii), len(WAVELENGTHS_UM), 2), np.nan)
        for index, rmod in enumerate(radii):
            raw[index] = _factors(rmod, spherical)
            if index % 30 == 0:
                print(
                    f"  {'spheres' if spherical else 'spheroids'},"
                    f" rmod = {rmod:7.4f} um"
                )
        # "Normalization to the smallest particles": the first radius is deep
        # in the Rayleigh regime, where the correction is a property of the
        # geometry alone.
        per_shape.append(
            xr.Dataset(
                {
                    "cts": (("rmod", "wl"), raw[:, :, 0] / raw[0, :, 0]),
                    "cbs": (("rmod", "wl"), raw[:, :, 1] / raw[0, :, 1]),
                },
                coords={"rmod": radii, "wl": WAVELENGTHS_UM},
            )
        )
    return xr.concat(
        per_shape,
        dim=xr.DataArray(["spheres", "spheroids"], dims="shape", name="shape"),
    )


def _factors(rmod: float, spherical: bool) -> np.ndarray:
    """The two uncorrected ratios at one mode radius, per wavelength."""
    mode = MicroParameters(
        wavelength=WAVELENGTHS_UM,
        n_real=N_REAL,
        n_imag=N_IMAG,
        shape=(
            Sphere()
            if spherical
            else SpheroidDistrFile(distr_filename=AR_KANDLER)
        ),
        psd=LognormalPSD(
            rm=rmod, sigma=SIGMA, n=1e6, rmin=R_MIN_UM, rmax=R_MAX_UM
        ),
        nonabs_fraction=NONABS_FRACTION,
    )
    point = run_point([mode], OUTPUTS, quiet=True, n_angles=N_ANGLES)

    theta = point["theta"].values
    weight = weights(theta)
    radians = np.radians(theta)
    phase = point["phase"].transpose("wl", "theta").values

    def integrate(name: str) -> np.ndarray:
        return np.trapezoid(weight[name] * phase, x=radians, axis=1)

    return np.stack(
        [
            integrate("total_ideal") / integrate("total_real"),
            integrate("back_ideal") / integrate("back_real"),
        ],
        axis=1,
    )


def report(data: xr.Dataset) -> None:
    """Print the three statements section 5.7 makes about the figure."""
    # The published panels stop at rmod = 1 um; the grid runs to 2.
    shown = data.sel(rmod=slice(None, 1.0))

    print("\n=== statements of section 5.7 ===")
    coarse = shown["cts"].sel(shape="spheres").isel(rmod=-1)
    print(f"Cts at rmod = 1 um, spheres: {coarse.values.round(2)}")
    print("  article: 'underestimates total scattering by a factor of ~ 2'")

    effect = (
        shown["cts"].sel(shape="spheroids") / shown["cts"].sel(shape="spheres")
        - 1.0
    )
    print(f"\nshape effect on Cts, at most: {float(abs(effect).max()):.1%}")
    print("  article: 'less than 3 %'")

    effect = (
        shown["cbs"].sel(shape="spheroids") / shown["cbs"].sel(shape="spheres")
        - 1.0
    )
    print(f"shape effect on Cbs, at most: {float(abs(effect).max()):.1%}")
    print("  article: 'The maximum shape effect on Cbs is 7 %'")


def plot(data: xr.Dataset, path: Path) -> None:
    """Draw the two panels in the layout of the published figure."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for wavelength, colour in zip(WAVELENGTHS_UM, COLOURS):
        for shape_name, style in (("spheres", "-"), ("spheroids", ":")):
            label = f"$\\lambda$ = {wavelength * 1e3:3.0f} nm, {shape_name}"
            axes[0].semilogx(
                data["rmod"],
                data["cts"].sel(shape=shape_name, wl=wavelength),
                style,
                color=colour,
                label=label,
            )
            axes[1].semilogx(
                data["rmod"],
                data["cbs"].sel(shape=shape_name, wl=wavelength),
                style,
                color=colour,
            )
    axes[0].set_yticks(np.arange(0.0, 3.1, 0.2))
    axes[0].set_ylim(1.0, 2.0)
    axes[0].set_ylabel("$C_{ts}$")
    axes[0].legend(fontsize=8)
    axes[0].text(0.011, 1.02, "(a)")
    axes[1].set_ylim(0.8, 1.2)
    axes[1].set_xlim(0.01, 1.0)
    axes[1].set_ylabel("$C_{bs}$")
    axes[1].set_xlabel("mode radius $r_{mod}$ [$\\mu$m]")
    axes[1].text(0.011, 0.81, "(b)")
    for axis in axes:
        axis.grid(which="both", alpha=0.3)
    fig.suptitle("Gasteiger and Wiegner (2018), Figure 10, recomputed")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print(f"\nwrote {path}")


def main() -> None:
    print("Computing the nephelometer correction factors:")
    data = compute()
    report(data)
    plot(data, FIGURES / "gasteiger_fig10.png")


if __name__ == "__main__":
    main()
