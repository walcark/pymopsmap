"""Mix: an external mixture of species."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

import xarray as xr

from pymopsmap.engine.outputs import DEFAULT_OUTPUT, OutputRequest

from .catalog import CatalogSpecie
from .specie import Specie, load

Entry = CatalogSpecie | Specie | str | Path


class Currency(StrEnum):
    """How the weights of a mixture are expressed."""

    NUMBER = "number"
    MASS = "mass"
    OPTICAL_DEPTH = "optical_depth"


@dataclass(frozen=True)
class Mix:
    """
    Several species mixed externally, each with its own concentration.

    Optical properties are linear in the number concentration, so each species
    is computed once and the results are combined afterwards. Changing a
    concentration then costs nothing: it never triggers a run.

    Parameters
    ----------
    entries : dict
        Catalogue entries, paths or already loaded species, mapped to their
        number concentration in m-3.
    """

    entries: dict[Entry, float]
    currency: Currency = Currency.NUMBER
    wl_ref: float | None = None
    rh_ref: float | None = None
    species: list[Specie] = field(init=False)

    def __post_init__(self) -> None:
        if not self.entries:
            raise ValueError("Cannot build an empty mixture.")
        object.__setattr__(
            self,
            "species",
            [
                entry if isinstance(entry, Specie) else load(entry)
                for entry in self.entries
            ],
        )

    @classmethod
    def from_mass(cls, entries: dict[Entry, float], *, rh_ref: float) -> Mix:
        """
        Build a mixture from mass concentrations, in kg m-3.

        Parameters
        ----------
        entries : dict
            Species mapped to their mass concentration.
        rh_ref : float
            The humidity at which those masses hold. It is required: a wet
            particle carries water, so a mass without a humidity does not
            identify a state. Dry masses mean ``rh_ref=0``.
        """
        return cls(entries, currency=Currency.MASS, rh_ref=rh_ref)

    @classmethod
    def from_optical_depth(
        cls, entries: dict[Entry, float], *, wl_ref: float, rh_ref: float
    ) -> Mix:
        """
        Build a mixture from fractions of the total optical depth.

        The fractions are normalised, and the resulting extinction is
        normalised too: it sums to one at ``wl_ref``, to be multiplied by the
        optical depth of interest.

        Parameters
        ----------
        entries : dict
            Species mapped to their share of the total optical depth.
        wl_ref : float
            The wavelength the fractions refer to, in micrometres.
        rh_ref : float
            The humidity at which the fractions were observed. It is required
            and has no default: an optical depth fraction is a derived
            quantity that depends on humidity, because species do not grow
            alike. What stays fixed is the number concentration, which this
            inverts the fractions into.
        """
        return cls(
            entries,
            currency=Currency.OPTICAL_DEPTH,
            wl_ref=wl_ref,
            rh_ref=rh_ref,
        )

    @property
    def weights(self) -> dict[str, float]:
        """The requested weight per species, in the currency it was given."""
        return dict(zip((s.name for s in self.species), self.entries.values()))

    def contributions(self, result: xr.Dataset) -> dict[str, float]:
        """
        Share of the total extinction held by each species, at ``wl_ref``.

        Parameters
        ----------
        result : xr.Dataset
            A result produced by this mixture.

        Returns
        -------
        dict
            One share per species, summing to one.
        """
        wl = self.wl_ref if self.wl_ref is not None else float(result["wl"][0])
        total = float(result["kext"].sel(wl=wl, method="nearest"))
        scales = result.attrs["scales"]
        return {
            specie.name: scale
            * float(specie.compute(wl=[wl], rh=self.rh_ref)["kext"].isel(wl=0))
            / total
            for specie, scale in zip(self.species, scales)
        }

    def _scales(self, quiet: bool) -> list[float]:
        """
        Factor applied to each species so it reaches its requested weight.

        The species file already holds the ratio between its modes, so the
        mixture scales the whole species rather than replacing an amplitude
        that a single number could not express.

        Number concentrations resolve straight away. Masses and optical depth
        fractions are derived quantities, so they are inverted against a
        reference run of each species. That inversion is only possible because
        the species are computed separately: it needs their individual results
        before the weights are known.
        """
        requested = list(self.entries.values())
        if self.currency is Currency.NUMBER:
            return [
                value / specie.amplitude
                for specie, value in zip(self.species, requested)
            ]

        references = [
            specie.compute(
                wl=[self._reference_wavelength(specie)],
                rh=self.rh_ref,
                quiet=quiet,
            )
            for specie in self.species
        ]

        if self.currency is Currency.MASS:
            return [
                value / float(reference["mass_conc"])
                for value, reference in zip(requested, references)
            ]

        total = sum(requested)
        return [
            (value / total) / float(reference["kext"].isel(wl=0))
            for value, reference in zip(requested, references)
        ]

    def _reference_wavelength(self, specie: Specie) -> float:
        """Where to evaluate a species when inverting a derived weight."""
        if self.wl_ref is not None:
            return self.wl_ref
        # Mass is wavelength independent, so any covered wavelength works.
        low, high = specie.wl_range
        return 0.5 * (low + high)

    def compute(
        self,
        wl: list[float],
        rh: float | list[float] | None = None,
        outputs: OutputRequest = DEFAULT_OUTPUT,
        quiet: bool = False,
    ) -> xr.Dataset:
        """
        Compute the optical properties of the mixture.

        Parameters
        ----------
        wl : list of float
            Wavelengths in micrometres, shared by every species.
        rh : float or list of float, optional
            Relative humidity in percent.
        outputs : OutputRequest
            Which families of optical properties to compute.
        quiet : bool
            Silence the dataset download progress bars.

        Returns
        -------
        xr.Dataset
        """
        from pymopsmap.engine.outputs import combine

        scales = self._scales(quiet)
        results = [
            specie.compute(wl=wl, rh=rh, outputs=outputs, quiet=quiet)
            for specie in self.species
        ]
        mixed = combine(results, weights=scales)
        mixed.attrs["species"] = [specie.name for specie in self.species]
        mixed.attrs["concentrations"] = [
            scale * specie.amplitude
            for scale, specie in zip(scales, self.species)
        ]
        mixed.attrs["scales"] = scales
        return mixed
