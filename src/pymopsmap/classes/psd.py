"""
psd.py

Author  : Kévin Walcarius
Date    : 2025-11-19
Version : 1.0
License : MIT
Summary : Defines the different input Particle Size Distribution (PSD)
          available in MOPSMAP. Also defined a type Size to simplify
          future loading with hydra.
"""

from typing import Literal, Annotated, Union
from pydantic import BaseModel, Field, model_validator


# =================================================================================================
# 1) Fixed size distribution
# =================================================================================================
class FixedPSD(BaseModel):
    """
    Defines the radius r of a single particle (given in μm).
    The total particle number density n needs to be given in units of m^-3.
    """

    type: Literal["fixed"] = "fixed"
    radius: float = Field(gt=0)
    n: float = Field(gt=0)

    @property
    def command(self) -> str:
        return f"size {self.radius} {self.n}"


# =================================================================================================
# 2) Log-normal distribution
# =================================================================================================
class LognormalPSD(BaseModel):
    """
    Defines a log-normal size distribution (particle number density per particle radius
    interval) according to:

        n(r) = 1/sqrt(2) * n/ln(sigma) * 1/r * exp[-0.5 * (ln(r/rm) / ln(sigma))²)

    Particles in the radius range from rmin to rmax are covered.
    The total particle number density n0 needs to be given in units of m^-3. Note that
    the actual particle number density of the modeled ensemble may be lower than n0
    because of clipping at rmin and rmax.
    """

    type: Literal["lognormal"] = "lognormal"
    rm: float = Field(gt=0)
    sigma: float = Field(gt=1)
    n: float = Field(gt=0)
    rmin: float = Field(gt=0)
    rmax: float = Field(gt=0)

    @model_validator(mode="after")
    def check_rmin_rmax(self):
        if self.rmax <= self.rmin:
            raise ValueError(
                f"rmax > rmin required, got rmax={self.rmax} and rmin={self.rmin})."
            )
        return self

    @property
    def command(self) -> str:
        return (
            f"size log_normal {self.rm} {self.sigma} {self.n} {self.rmin} {self.rmax}"
        )


# =================================================================================================
# 2) Modified Gamma distribution
# =================================================================================================
class ModifiedGammaPSD(BaseModel):
    """
    Defines a modified gamma distribution (particle number density per particle radius
    interval) according to:

            n(r) = A * r^alpha * exp[-B * r^gamma]

    Particles in the radius range from rmin to rmax are covered.
    The parameter A needs to be given in units of m−3.
    """

    type: Literal["mod-gamma"] = "mod-gamma"
    A: float = Field(gt=0)
    B: float = Field(gt=0)
    alpha: float
    gamma: float
    rmin: float = Field(gt=0)
    rmax: float = Field(gt=0)

    @model_validator(mode="after")
    def check_bounds(self):
        if self.rmax <= self.rmin:
            raise ValueError(
                f"rmax > rmin required, got rmax={self.rmax} and rmin={self.rmin})."
            )
        return self

    @property
    def command(self) -> str:
        return f"size mod_gamma {self.A} {self.rmin} {self.rmax} {self.alpha} {self.B} {self.gamma}"


# =================================================================================================
# 4) File-defined distribution
# =================================================================================================
class FileDefinedPSD(BaseModel):
    type: Literal["bin-file"] = "bin-file"
    filename: str

    @property
    def command(self) -> str:
        return f"size bin_file '{self.filename}'"


# =================================================================================================
# Annotated union
# =================================================================================================
Size = Annotated[
    Union[FixedPSD, LognormalPSD, FileDefinedPSD],
    Field(discriminator="type"),
]
