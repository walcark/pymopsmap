"""Particle shapes, and the MOPSMAP command each one writes."""

from typing import Annotated, Literal

from pydantic import BaseModel, Field


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
        return (
            f"shape irregular_overlay {self.distr_filename}"
            f" {{self.xmin}} {{self.xmax}}"
        )


Shape = Annotated[
    Sphere
    | Spheroid
    | SpheroidLognormal
    | SpheroidDistrFile
    | Irregular
    | IrregularDistrFile
    | IrregularOverlay,
    Field(discriminator="type"),
]
