"""
Reproduce Figure 8 of Gasteiger and Wiegner (2018), GMD 11:2739.

    pixi run -e dev python -m scripts.validation.gasteiger_fig8

"Lognormal size distributions (SD) with same rmod, sigma, N0, and rmax
assuming different size equivalences for aggregate particles (shape D,
xi_vc = 0.8708)."

Nothing here calls MOPSMAP. The figure shows what the ``size_equ`` option
means before any optics are computed: the same four numbers describe three
different populations of particles, depending on whether the radius given is
read as the sphere of equal cross section, of equal volume, or of equal
volume-to-cross-section ratio.

Equations 3 and 4 of the article relate the three:

    xi_vc = rv / rc        and        rvcr = xi_vc^3 * rc

so a radius given as rv sits at rc = rv / xi_vc, and one given as rvcr sits at
rc = rvcr / xi_vc^3. The article states the second: "rmod = xi_vc^-3 * 0.5 um
= 0.757 um", which this script checks.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

FIGURES = Path(__file__).resolve().parents[2] / "docs" / "figures"

# Shape D of Gasteiger et al. (2011b), the aggregate of Table 5.
XI_VC = 0.8708

# The mode of section 5.4, read as a cross-section-equivalent radius.
R_MOD_UM = 0.5
SIGMA = 2.0
R_MIN_UM = 0.001
R_MAX_UM = 1.75

# The amplitude Table 5 is computed with: the article writes N0 = 10^3.66 cm-3
# beside "a concentration of N = 100 cm-3 in the range from rmin to rmax", and
# the two disagree. The second is the one MOPSMAP echoes back, so it is the
# one used here, converted to the total the lognormal is normalised on.
IN_RANGE_FRACTION = 0.9646465211511944
N0_CM3 = 100.0 / IN_RANGE_FRACTION

# The three readings, and how each one maps onto a cross-section-equivalent
# radius. The colours are the ones the caption gives.
EQUIVALENCES = {
    "cs": ("SD assuming $r_c$", 1.0, "black"),
    "vol": ("SD assuming $r_v$", XI_VC, "red"),
    "vol_cs_ratio": ("SD assuming $r_{vcr}$", XI_VC**3, "green"),
}


def dn_drc(radii: np.ndarray, r_mod: float) -> np.ndarray:
    """
    dN/dr of a lognormal, as log_distr.f90 writes it.

    Parameters
    ----------
    radii : ndarray
        Cross-section-equivalent radii in micrometres.
    r_mod : float
        Modal radius in the same units.

    Returns
    -------
    ndarray
        Number density per micrometre, in cm-3 um-1.
    """
    log_sigma = np.log(SIGMA)
    return (
        N0_CM3
        / (np.sqrt(2.0 * np.pi) * radii * log_sigma)
        * np.exp(-(np.log(radii / r_mod) ** 2) / (2.0 * log_sigma**2))
    )


# The peaks of the three published curves, read off the figure at 250 dpi.
PUBLISHED_PEAKS = {"cs": 380.0, "vol": 332.0, "vol_cs_ratio": 252.0}


def report() -> None:
    """Check the conversions, and the one thing that does not agree."""
    print("=== section 5.4 ===")
    converted = R_MOD_UM / XI_VC**3
    print(f"rmod read as rvcr sits at rc = {converted:.3f} um")
    print("  article: 'xi_vc^-3 * 0.5 um = 0.757 um'")

    # Mass follows the cube of the radius at fixed number, which is what
    # Table 5 shows between its first and last columns.
    ratio = (converted / R_MOD_UM) ** 3
    print(f"\nmass ratio between the rvcr and the rc reading: {ratio:.2f}")
    print("  Table 5: 1103 / 318 ug m-3 = 3.47")

    print("\n=== peak of each curve ===")
    print(f"{'':16s}{'published':>11s}{'computed':>11s}{'ratio':>9s}")
    for key, (label, factor, _) in EQUIVALENCES.items():
        radii = np.logspace(-2, 1, 20000)
        peak = float(dn_drc(radii, R_MOD_UM / factor).max())
        published = PUBLISHED_PEAKS[key]
        name = label.split()[-1].strip("$")
        print(
            f"{name:16s}{published:11.0f}{peak:11.1f}{published / peak:9.3f}"
        )
    print(f"\nsqrt(2 pi) = {np.sqrt(2.0 * np.pi):.3f}")


def plot(path: Path) -> None:
    """Draw the three distributions on the cross-section-equivalent axis."""
    import matplotlib.pyplot as plt

    radii = np.linspace(0.02, 3.0, 600)

    fig, axis = plt.subplots(figsize=(8, 6))
    axis.set_ylim(0.0, 500.0)
    for label, factor, colour in EQUIVALENCES.values():
        r_mod = R_MOD_UM / factor
        values = dn_drc(radii, r_mod)
        # The cutoff travels with the reading, exactly as rmod does.
        values[radii > R_MAX_UM / factor] = np.nan
        axis.plot(radii, values, color=colour, label=label)

    axis.set_xlim(0.0, 3.0)
    axis.annotate(
        "published peaks: 380, 332, 252,\n"
        "a factor $\\sqrt{2\\pi}$ above the\nlognormal normalisation",
        xy=(1.55, 420.0),
        fontsize=8,
    )
    axis.set_xlabel("$r_c$ [$\\mu$m]")
    axis.set_ylabel("$dN/dr_c$ [cm$^{-3}$ $\\mu$m$^{-1}$]")
    axis.legend()
    axis.grid(alpha=0.3)

    # The red and green axes of the published figure: the same curves read in
    # terms of the other two size definitions.
    for offset, (factor, colour) in enumerate(
        ((XI_VC, "red"), (XI_VC**3, "green")), start=1
    ):
        twin = axis.secondary_xaxis(
            1.0 + 0.09 * (offset - 1),
            functions=(lambda x, f=factor: x * f, lambda x, f=factor: x / f),
        )
        twin.tick_params(colors=colour)
        twin.set_xlabel(
            "$r_v$" if colour == "red" else "$r_{vcr}$", color=colour
        )

    fig.suptitle("Gasteiger and Wiegner (2018), Figure 8, recomputed")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print(f"\nwrote {path}")


def main() -> None:
    report()
    plot(FIGURES / "gasteiger_fig8.png")


if __name__ == "__main__":
    main()
