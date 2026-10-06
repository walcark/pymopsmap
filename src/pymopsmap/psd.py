"""Particle size distributions, and the command each one writes."""

from enum import StrEnum
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from pymopsmap.utils import SortedPosFloat64List
from pymopsmap.varying import AboveOne, NonNegative, Positive, Varying


class FixedPSD(BaseModel):
    """
    Defines the radius r of a single particle (given in μm).
    The total particle number density n needs to be given in units of m^-3.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    type: Literal["fixed"] = "fixed"
    radius: Positive
    n: NonNegative

    @property
    def command(self) -> str:
        return f"size {self.radius} {self.n}"


class LognormalPSD(BaseModel):
    """
    Defines a log-normal size distribution.

    Particles in the radius range from rmin to rmax are covered.
    The total particle number density n0 needs to be given in units of m^-3.
    Note that the actual particle number density of the modeled ensemble may
    be lower than n0 because of clipping at rmin and rmax.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    type: Literal["lognormal"] = "lognormal"
    rm: Positive
    sigma: AboveOne
    n: NonNegative
    rmin: Positive
    rmax: Positive

    @model_validator(mode="after")
    def check_rmin_rmax(self):
        if self.rmax <= self.rmin:
            raise ValueError(
                f"rmax > rmin required, "
                f"got rmax={self.rmax} and rmin={self.rmin})."
            )
        return self

    @property
    def command(self) -> str:
        return (
            f"size log_normal {self.rm} {self.sigma}"
            f" {self.n} {self.rmin} {self.rmax}"
        )


class ModifiedGammaPSD(BaseModel):
    """
    Defines a modified gamma distribution:

            n(r) = A * r^alpha * exp[-B * r^gamma]

    Particles in the radius range from rmin to rmax are covered.
    The parameter A needs to be given in units of m−3.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    type: Literal["mod-gamma"] = "mod-gamma"
    A: Positive
    B: Positive
    alpha: Varying
    gamma: Varying
    rmin: Positive
    rmax: Positive

    @model_validator(mode="after")
    def check_bounds(self):
        if self.rmax <= self.rmin:
            raise ValueError(
                f"rmax > rmin required, "
                f"got rmax={self.rmax} and rmin={self.rmin})."
            )
        return self

    @property
    def command(self) -> str:
        return (
            f"size mod_gamma {self.A} {self.rmin} {self.rmax}"
            f" {self.alpha} {self.B} {self.gamma}"
        )


class FileDefinedPSD(BaseModel):
    type: Literal["bin-file"] = "bin-file"
    filename: str

    @property
    def command(self) -> str:
        return f"size bin_file '{self.filename}'"


class DistrType(StrEnum):
    DNDR = "dndr"
    DNDLNR = "dndlnr"
    DNDLOGR = "dndlogr"
    DADR = "dadr"
    DADLNR = "dadlnr"
    DADLOGR = "dadlogr"
    DVDR = "dvdr"
    DVDLNR = "dvdlnr"
    DVDLOGR = "dvdlogr"


class DistrListPSD(BaseModel):
    type: Literal["distr-list"] = "distr-list"
    radii: SortedPosFloat64List
    concentrations: list[float]
    distr_type: DistrType

    @model_validator(mode="after")
    def check_lengths(self):
        if len(self.radii) != len(self.concentrations):
            raise ValueError(
                f"radii and concentrations must have the same length, "
                f"got {len(self.radii)} vs {len(self.concentrations)}"
            )
        if len(self.radii) < 2:
            raise ValueError(
                "DistrListPSD requires at least 2 radius/concentration pairs"
            )
        return self

    @property
    def command(self) -> str:
        pairs = " ".join(
            f"{r} {c}" for r, c in zip(self.radii, self.concentrations)
        )
        return f"size distr_list {self.distr_type.value} {pairs}"


PSD = Annotated[
    FixedPSD | LognormalPSD | FileDefinedPSD | ModifiedGammaPSD | DistrListPSD,
    Field(discriminator="type"),
]
