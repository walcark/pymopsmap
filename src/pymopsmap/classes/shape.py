# pymopsmap/classes/shape.py
from typing import Literal
from pydantic import BaseModel, Field

# =================================================================================================
# 1) SPHERE
# =================================================================================================
class Sphere(BaseModel):
    type: Literal["sphere"] = "sphere"

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} shape sphere"


# =================================================================================================
# 2) SPHEROID
# =================================================================================================
class Spheroid(BaseModel):
    type: Literal["spheroid"] = "spheroid"
    mode: Literal["oblate", "prolate"]
    aspect_ratio: float = Field(gt=1)

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} shape spheroid {self.mode} {self.aspect_ratio}"


# =================================================================================================
# 3) LOG-NORMAL SPHEROID
# =================================================================================================
class SpheroidLognormal(BaseModel):
    type: Literal["spheroid-lognormal"] = "spheroid-lognormal"
    zeta1: float = Field(ge=0, le=1)
    zeta2: float = Field(ge=0, le=1)
    aspect_ratio: float = Field(ge=1.2, le=5.0)
    sigma_ar: float = Field(gt=0)

    def to_section(self, num: int | None = None) -> str:
        return (
            f"mode {num} shape spheroid log_normal "
            f"{self.zeta1} {self.zeta2} {self.aspect_ratio} {self.sigma_ar}"
        )


# =================================================================================================
# 4) GENERALIZED SPHEROID
# =================================================================================================
class SpheroidDistrFile(BaseModel):
    type: Literal["spheroid-distr-file"] = "spheroid-distr-file"
    distr_filename: str

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} shape spheroid distr_file {self.distr_filename}"


# =================================================================================================
# 5) PRE-DEFINED SHAPES (from Gasteiger at al. (2011))
# =================================================================================================
class Irregular(BaseModel):
    type: Literal["irregular"] = "irregular"
    shape_id: Literal["A", "B", "C", "D", "E", "F"]

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} shape irregular {self.shape_id}"


# =================================================================================================
# 5) MULTIPLE PRE-DEFINED SHAPES FROM FILE
# =================================================================================================
class IrregularDistrFile(BaseModel):
    type: Literal["irregular-distr-file"] = "irregular-distr-file"
    distr_filename: str

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} shape irregular distr_file {self.distr_filename}"


# =================================================================================================
# 5) IRREGULATE OVERLAY
# =================================================================================================
class IrregularOverlay(BaseModel):
    type: Literal["irregular-overlay"] = "irregular-overlay"
    distr_filename: str
    xmin: float
    xmax: float

    def to_section(self, num: int | None = None) -> str:
        return f"mode {num} shape irregular_overlay {self.distr_filename} {self.xmin} {self.xmax}"