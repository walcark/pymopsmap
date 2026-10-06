"""Mode: one aerosol mode described in Python rather than read from a file."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import xarray as xr

from pymopsmap.psd import PSD, FixedPSD
from pymopsmap.shapes import Shape, Spheroid

from .schema import AMPLITUDE_FIELD

# The aspect ratios the optical data set holds, Table 1 of Gasteiger and
# Wiegner (2018). A ratio between two of them is split over both, which is
# what init_shape.f90 does with one that reaches MOPSMAP unbinned.
ASPECT_RATIO_GRID = (
    1.0,
    1.2,
    1.4,
    1.6,
    1.8,
    2.0,
    2.2,
    2.4,
    2.6,
    2.8,
    3.0,
    3.4,
    3.8,
    4.2,
    4.6,
    5.0,
)


@dataclass(frozen=True)
class Mode:
    """
    One mode of a hand-built species.

    Any numeric parameter accepts a ``DataArray`` where a float would go,
    which sweeps it: distinct dimensions multiply, and two parameters sharing
    a dimension vary together. They are declared on the species rather than on
    the compute call, so a parameter name is never ambiguous between modes.

    Parameters
    ----------
    shape : Shape
        Particle shape.
    psd : PSD
        Size distribution, amplitude included.
    n_real, n_imag : float or DataArray
        Refractive index, on the wavelengths given to ``compute``.
    density_dry : float, optional
        Dry material density in g cm-3.
    kappa : float, optional
        Hygroscopicity. Its presence makes the species grow with humidity.
    nonabs_fraction : float
        The share of the particles that does not absorb at all, the rest
        absorbing the more for it so the average index holds. Section 3.1 of
        Gasteiger and Wiegner (2018).
    """

    shape: Shape
    psd: PSD
    n_real: Any
    n_imag: Any
    density_dry: float | None = None
    kappa: float | None = None
    nonabs_fraction: float = 0.0
    name: str = field(default="only")

    @classmethod
    def from_particles(
        cls,
        radii: Sequence[float],
        aspect_ratios: Sequence[float],
        *,
        n_real: Any,
        n_imag: Any,
        weights: Sequence[float] | None = None,
        spheroid: str = "prolate",
        r_max: float | None = None,
        density_dry: float | None = None,
        kappa: float | None = None,
        nonabs_fraction: float = 0.0,
        grid: Sequence[float] = ASPECT_RATIO_GRID,
    ) -> list[Mode]:
        """
        Turn a table of measured particles into the modes MOPSMAP runs.

        One mode per distinct size and aspect ratio, the amplitude counting
        how many particles fell there. A ratio between two grid points is
        split over both, which is what MOPSMAP does with an unbinned one
        (``init_shape.f90``), so the contributions are the same and there are
        fewer of them: a microscopy table of forty thousand particles comes
        out as a few thousand modes.

        Parameters
        ----------
        radii : sequence of float
            Particle radii in micrometres, one per particle. Microscopy
            tables usually give a diameter, which is twice this.
        aspect_ratios : sequence of float
            Ratio of the long axis to the short one, at least one.
        n_real, n_imag : float or sequence of float
            Refractive index, shared by every particle.
        weights : sequence of float, optional
            What each row counts for. One each by default, so the amplitudes
            come out as particle counts.
        spheroid : {'prolate', 'oblate'}
            Which way the spheroids are stretched.
        r_max : float, optional
            Radius above which particles are clipped rather than dropped, as
            section 5.8 of the article does at 47.5 um. Ratios are always
            clipped at the largest grid point.
        density_dry, kappa, nonabs_fraction
            Carried onto every mode.
        grid : sequence of float
            The aspect ratio grid to bin onto.

        Returns
        -------
        list of Mode
            Ready for ``Specie.custom``.

        Raises
        ------
        ValueError
            If the inputs disagree in length, or an aspect ratio is below one.
        """
        counts = _bin_particles(radii, aspect_ratios, weights, grid, r_max)
        return [
            cls(
                shape=Spheroid(mode=spheroid, aspect_ratio=ratio),  # type: ignore[arg-type]
                psd=FixedPSD(radius=radius, n=weight),
                n_real=n_real,
                n_imag=n_imag,
                density_dry=density_dry,
                kappa=kappa,
                nonabs_fraction=nonabs_fraction,
                name=f"mode_{index + 1}",
            )
            for index, ((radius, ratio), weight) in enumerate(
                sorted(counts.items())
            )
        ]

    @classmethod
    def from_index_distribution(
        cls,
        n_imag: Sequence[float],
        weights: Sequence[float],
        *,
        n_real: Any,
        shape: Shape,
        psd: PSD,
        density_dry: float | None = None,
        kappa: float | None = None,
    ) -> list[Mode]:
        """
        Spread one mode over a distribution of imaginary refractive indices.

        Particles of one aerosol are not equally absorbing, and section 5.6 of
        the article shows that averaging the index is not the same as
        averaging what it produces. One mode per bin, its amplitude cut by the
        share of particles in that bin, is how MOPSMAP is told the difference.

        Parameters
        ----------
        n_imag : sequence of float
            The imaginary index of each bin, its midpoint for a histogram.
        weights : sequence of float
            How many particles fall in each bin. Normalised here, so counts
            work as well as fractions.
        n_real : float or sequence of float
            The real index, shared by every bin.
        shape : Shape
            Particle shape, shared by every bin.
        psd : PSD
            The size distribution of the whole mode. Its amplitude is what
            gets divided between the bins.
        density_dry, kappa
            Carried onto every mode.

        Returns
        -------
        list of Mode
            One per non-empty bin, ready for ``Specie.custom``.

        Raises
        ------
        ValueError
            If the two sequences disagree in length, or every weight is zero.
        """
        if len(n_imag) != len(weights):
            raise ValueError(
                f"n_imag and weights must have the same length, got "
                f"{len(n_imag)} and {len(weights)}."
            )
        total = float(sum(weights))
        if total <= 0.0:
            raise ValueError("the weights of a distribution cannot sum to 0.")

        amplitude = AMPLITUDE_FIELD[psd.type]
        if amplitude is None:
            raise ValueError(
                f"a {psd.type} size distribution keeps its amplitude in a "
                "file, so it cannot be split between index bins."
            )
        whole = float(getattr(psd, amplitude))

        modes = []
        for index, (imaginary, weight) in enumerate(zip(n_imag, weights)):
            if weight <= 0.0:
                continue
            modes.append(
                cls(
                    shape=shape,
                    psd=psd.model_copy(
                        update={amplitude: whole * float(weight) / total}
                    ),
                    n_real=n_real,
                    n_imag=float(imaginary),
                    density_dry=density_dry,
                    kappa=kappa,
                    name=f"mode_{index + 1}",
                )
            )
        return modes

    def to_dataset(self, wl: list[float] | None = None) -> xr.Dataset:
        """Render the mode in the canonical schema."""
        variables: dict[str, Any] = {
            "n_real": _as_variable(self.n_real),
            "n_imag": _as_variable(self.n_imag),
        }
        for name in _fields_of(self.psd):
            variables[name] = _as_variable(getattr(self.psd, name))
        for name in _fields_of(self.shape):
            variables[name] = _as_variable(getattr(self.shape, name))
        if self.density_dry is not None:
            variables["density_dry"] = self.density_dry
        if self.kappa is not None:
            variables["kappa"] = self.kappa
        if self.nonabs_fraction:
            variables["nonabs_fraction"] = self.nonabs_fraction

        coords = {"wl": wl} if wl is not None else {}
        ds = xr.Dataset(variables, coords=coords)
        ds.attrs.update(psd_type=self.psd.type, shape_type=self.shape.type)
        _stamp_units(ds)
        return ds


def _bin_particles(
    radii: Sequence[float],
    aspect_ratios: Sequence[float],
    weights: Sequence[float] | None,
    grid: Sequence[float],
    r_max: float | None,
) -> dict[tuple[float, float], float]:
    """How much weight lands on each size and grid aspect ratio."""
    if len(radii) != len(aspect_ratios):
        raise ValueError(
            f"radii and aspect_ratios must have the same length, got "
            f"{len(radii)} and {len(aspect_ratios)}."
        )
    if weights is not None and len(weights) != len(radii):
        raise ValueError(
            f"weights must have the same length as radii, got "
            f"{len(weights)} and {len(radii)}."
        )

    size = np.asarray(radii, dtype=float)
    ratio = np.asarray(aspect_ratios, dtype=float)
    share = (
        np.ones(len(size))
        if weights is None
        else np.asarray(weights, dtype=float)
    )
    if (ratio < 1.0).any():
        raise ValueError(
            "an aspect ratio is the long axis over the short one, so it is "
            f"at least one; got {float(ratio.min()):g}."
        )

    points = np.asarray(grid, dtype=float)
    if r_max is not None:
        size = np.minimum(size, r_max)
    ratio = np.clip(ratio, points[0], points[-1])

    upper = np.clip(np.searchsorted(points, ratio, side="left"), 1, None)
    lower = upper - 1
    below = (points[upper] - ratio) / (points[upper] - points[lower])

    counts: dict[tuple[float, float], float] = defaultdict(float)
    for radius, low, high, fraction, weight in zip(
        size, lower, upper, below, share
    ):
        if fraction > 0.0:
            counts[(float(radius), float(points[low]))] += weight * fraction
        if fraction < 1.0:
            counts[(float(radius), float(points[high]))] += weight * (
                1.0 - fraction
            )
    return counts


def _fields_of(model) -> list[str]:
    return [name for name in type(model).model_fields if name != "type"]


def _as_variable(value: Any) -> Any:
    """Keep a DataArray as it is; anything else becomes a scalar or a list."""
    if isinstance(value, xr.DataArray):
        return value
    if isinstance(value, (list, tuple, np.ndarray)):
        return ("wl", np.asarray(value, dtype=float))
    return value


LENGTHS = {"rm", "rmin", "rmax", "radius"}


def _stamp_units(ds: xr.Dataset) -> None:
    for name in ds.variables:
        key = str(name)
        if key in LENGTHS or key == "wl":
            ds[key].attrs.update(units="um")
        elif key in AMPLITUDE_FIELD.values():
            ds[key].attrs.update(units="m-3")
        elif key == "density_dry":
            ds[key].attrs.update(units="g cm-3")
        elif key in ("n_real", "n_imag", "sigma", "kappa"):
            ds[key].attrs.update(units="1")
        elif key == "nonabs_fraction":
            ds[key].attrs.update(units="1")
    if "n_imag" in ds:
        ds["n_imag"].attrs.update(sign="positive")
