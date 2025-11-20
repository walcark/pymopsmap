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
    type: Literal["fixed"] = "fixed"
    radius: float = Field(gt=0)

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} size fixed {self.radius}"


# =================================================================================================
# 2) Log-normal distribution
# =================================================================================================
class LognormalPSD(BaseModel):
    type: Literal["lognormal"] = "lognormal"
    r_eff: float = Field(gt=0)
    sigma: float = Field(gt=1)
    n: float = Field(gt=0)
    rmin: float = Field(gt=0)
    rmax: float = Field(gt=0)

    @model_validator(mode="after")
    def check_rmin_rmax(self):
        if self.rmax <= self.rmin:
            raise ValueError(f"rmax ({self.rmax}) must be > rmin ({self.rmin}).")
        return self

    def to_section(self, num: int | None = None) -> str:
        return (
            f"mode {num} size log_normal {self.r_eff} "
            f"{self.sigma} {self.n} {self.rmin} {self.rmax}"
        )


# =================================================================================================
# 2) Modified Gamma distribution
# =================================================================================================
class ModifiedGammaPSD(BaseModel):
    type: Literal["mod-gamma"] = "mod-gamma"
    r_eff: float = Field(gt=0)
    v_eff: float = Field(gt=0)
    rmin: float = Field(gt=0)
    rmax: float = Field(gt=0)

    @model_validator(mode="after")
    def check_bounds(self):
        if self.rmax <= self.rmin:
            raise ValueError(f"rmax ({self.rmax}) must be > rmin ({self.rmin}).")
        return self

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} size mod_gamma {self.r_eff} {self.v_eff} {self.rmin} {self.rmax}"


# =================================================================================================
# 4) File-defined distribution
# =================================================================================================
class FileDefinedPSD(BaseModel):
    type: Literal["distr-file"] = "distr-file"
    distr_filename: str

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} size distr file {self.distr_filename}"


# =================================================================================================
# Annotated union
# =================================================================================================
Size = Annotated[
    Union[FixedPSD, LognormalPSD, FileDefinedPSD],
    Field(discriminator="type"),
]
