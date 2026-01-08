"""
microparams.py

Author  : Kévin Walcarius
Date    : 2025-01-08
Version : 1.0
License : MIT
Summary : Class to encapsulate the microphysical parameters of an
          atmospheric constituant.
"""

from pydantic import BaseModel, PositiveFloat, model_validator

from pymopsmap.utils import Float64List, PosFloat64List, SortedPosFloat64List
from typing import Annotated, Literal, Union

from pydantic import Field


# --------------------------------------------------------------------------
# Shapes
# --------------------------------------------------------------------------
class Sphere(BaseModel):
    type: Literal["sphere"] = "sphere"

    @property
    def command(self) -> str:
        return "shape sphere"


class Spheroid(BaseModel):
    type: Literal["spheroid"] = "spheroid"
    mode: Literal["oblate", "prolate"]
    aspect_ratio: float = Field(ge=1)

    @property
    def command(self) -> str:
        return f"shape spheroid {self.mode} {self.aspect_ratio}"


class SpheroidLognormal(BaseModel):
    type: Literal["spheroid-lognormal"] = "spheroid-lognormal"
    zeta1: float = Field(ge=0, le=1)
    zeta2: float = Field(ge=0, le=1)
    aspect_ratio: float = Field(ge=1.2, le=5.0)
    sigma_ar: float = Field(gt=0)

    @property
    def command(self) -> str:
        return (
            f"shape spheroid log_normal "
            f"{self.zeta1} {self.zeta2} {self.aspect_ratio} {self.sigma_ar}"
        )


class SpheroidDistrFile(BaseModel):
    type: Literal["spheroid-distr-file"] = "spheroid-distr-file"
    distr_filename: str

    @property
    def command(self) -> str:
        return f"shape spheroid distr_file {self.distr_filename}"


class Irregular(BaseModel):
    type: Literal["irregular"] = "irregular"
    shape_id: Literal["A", "B", "C", "D", "E", "F"]

    @property
    def command(self) -> str:
        return f"shape irregular {self.shape_id}"


class IrregularDistrFile(BaseModel):
    type: Literal["irregular-distr-file"] = "irregular-distr-file"
    distr_filename: str

    @property
    def command(self) -> str:
        return f"shape irregular distr_file {self.distr_filename}"


class IrregularOverlay(BaseModel):
    type: Literal["irregular-overlay"] = "irregular-overlay"
    distr_filename: str
    xmin: float
    xmax: float

    @property
    def command(self) -> str:
        return f"shape irregular_overlay {self.distr_filename} {self.xmin} {self.xmax}"


Shape = Annotated[
    Union[
        Sphere,
        Spheroid,
        SpheroidLognormal,
        SpheroidDistrFile,
        Irregular,
        IrregularDistrFile,
        IrregularOverlay,
    ],
    Field(discriminator="type"),
]


# --------------------------------------------------------------------------
# Particle size distributions
# --------------------------------------------------------------------------
class FixedPSD(BaseModel):
    """
    Defines the radius r of a single particle (given in μm).
    The total particle number density n needs to be given in units of m^-3.
    """

    type: Literal["fixed"] = "fixed"
    radius: float = Field(gt=0)
    n: float = Field(ge=0)

    @property
    def command(self) -> str:
        return f"size {self.radius} {self.n}"


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
    n: float = Field(ge=0)
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
        return f"size log_normal {self.rm} {self.sigma} {self.n} {self.rmin} {self.rmax}"


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


class FileDefinedPSD(BaseModel):
    type: Literal["bin-file"] = "bin-file"
    filename: str

    @property
    def command(self) -> str:
        return f"size bin_file '{self.filename}'"


PSD = Annotated[
    Union[FixedPSD, LognormalPSD, FileDefinedPSD, ModifiedGammaPSD],
    Field(discriminator="type"),
]


# =================================================================================================
# MicroParameters
# =================================================================================================
class MicroParameters(BaseModel):
    wavelength: SortedPosFloat64List
    n_real: PosFloat64List
    n_imag: Float64List
    shape: Shape
    psd: PSD
    kappa: PositiveFloat | None = None
    density: PositiveFloat | None = None

    def command(self, num: int | None = None) -> str:
        mode: str = f"mode {num} "
        string = (
            mode
            + self.shape.command
            + "\n"
            + mode
            + self.psd.command
            + "\n"
            + mode
            + self.refr_index.command
        )
        if self.kappa is not None:
            string += "\n" + mode + f"kappa {self.kappa}"

        if self.density is not None:
            string += "\n" + mode + f"density {self.density}"

        return string
