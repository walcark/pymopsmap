"""Particle shapes, and the MOPSMAP command each one writes."""

from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field

from pymopsmap.varying import (
    AtLeastOne,
    Positive,
    UnitInterval,
    Varying,
    aspect_ratio_range,
)


class Sphere(BaseModel):
    type: Literal["sphere"] = "sphere"

    @property
    def command(self) -> str:
        return "shape sphere"


class Spheroid(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    type: Literal["spheroid"] = "spheroid"
    mode: Literal["oblate", "prolate"]
    aspect_ratio: AtLeastOne

    @property
    def command(self) -> str:
        return f"shape spheroid {self.mode} {self.aspect_ratio}"


class SpheroidLognormal(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    type: Literal["spheroid-lognormal"] = "spheroid-lognormal"
    zeta1: UnitInterval
    zeta2: UnitInterval
    aspect_ratio: aspect_ratio_range(1.2, 5.0)  # type: ignore[valid-type]
    sigma_ar: Positive

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
    model_config = ConfigDict(arbitrary_types_allowed=True)
    type: Literal["irregular-overlay"] = "irregular-overlay"
    distr_filename: str
    xmin: Varying
    xmax: Varying

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
