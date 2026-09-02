"""Mix: an external mixture of species."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import xarray as xr

from pymopsmap.models.output_request import DEFAULT_OUTPUT, OutputRequest

from .catalog import CamsSpecie
from .specie import Specie, load

Entry = CamsSpecie | Specie | str | Path


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

    @property
    def weights(self) -> dict[str, float]:
        """Requested concentration per species, in m-3."""
        return dict(zip((s.name for s in self.species), self.entries.values()))

    def _scales(self) -> list[float]:
        """
        Factor applied to each species so it reaches its requested total.

        The species file already holds the ratio between its modes, so the
        mixture scales the whole species rather than replacing an amplitude
        that a single number could not express.
        """
        return [
            requested / specie.amplitude
            for specie, requested in zip(self.species, self.entries.values())
        ]

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
        from pymopsmap.engine.output_format import combine

        results = [
            specie.compute(wl=wl, rh=rh, outputs=outputs, quiet=quiet)
            for specie in self.species
        ]
        return combine(results, weights=self._scales())
