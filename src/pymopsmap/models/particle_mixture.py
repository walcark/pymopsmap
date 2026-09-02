"""ParticleMixture and ParametricSweep: particle system models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from .microparams import MicroParameters

if TYPE_CHECKING:
    from pymopsmap.models.optiprops import OptiProps
    from pymopsmap.models.output_request import OutputRequest


@dataclass
class ParticleMixture:
    """A particle population described by one or more set of micro-parameters.

    For instance, the Sea Salt CAMS specie is modelled as a bi-modal particle
    size distribution, each mode being a log-normal Particle Size Distribution.
    """

    modes: list[MicroParameters]

    def __post_init__(self):
        """Ensure broadcasting of micro-parameters modes into a list."""
        if isinstance(self.modes, MicroParameters):
            self.modes = [self.modes]

    def compute(
        self,
        output_types: OutputRequest | None = None,
        rh: float | None = None,
        quiet: bool = False,
    ) -> OptiProps:
        """Compute the optical properties of the mixture."""
        from pymopsmap.engine import compute_optical_properties

        return compute_optical_properties(
            self,
            output_types=output_types,
            rh=rh,
            quiet=quiet,
        )


@dataclass
class ParametricSweep:
    """A parametric sweep over multiple ParticleMixture configurations."""

    mixtures: list[ParticleMixture] = field(default_factory=list)
    params: list[dict[str, float]] = field(default_factory=list)

    def add(self, mixture: ParticleMixture, params: dict[str, float]) -> None:
        if self.mixtures:
            ref_wl = self.mixtures[0].modes[0].wavelength
            for mp in mixture.modes:
                if list(mp.wavelength) != list(ref_wl):
                    raise ValueError(
                        "All modes in a ParametricSweep must share the same "
                        f"wavelength array. Got {mp.wavelength} vs {ref_wl}."
                    )
        self.mixtures.append(mixture)
        self.params.append(params)

    def compute(
        self,
        output_types: OutputRequest | None = None,
        rh: float | None = None,
        quiet: bool = False,
    ) -> OptiProps:
        from pymopsmap.engine import compute_optical_properties

        return compute_optical_properties(
            self, output_types=output_types, rh=rh, quiet=quiet
        )
