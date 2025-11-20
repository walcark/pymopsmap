from pydantic import BaseModel, field_validator

from .refractive_index import RefractiveIndex
from .shape import Shape
from .psd import Size


class Mode(BaseModel):
    shape: Shape
    size: Size
    refr_index: RefractiveIndex
    kappa: float | None = None
    density: float | None = None

    @field_validator("kappa", "density")
    @classmethod
    def check_positive(cls, v, info):
        if v is None:
            return v
        if v <= 0:
            raise ValueError(f"{info.field_name} must be > 0 (got {v})")
        return v

    def command(self, num: int | None = None) -> str:
        """
        Convert a Mode to MOPSMAP section text.
        """
        mode: str = f"mode {num}"
        string = (
            mode
            + self.shape.command
            + "\n"
            + mode
            + self.size.command
            + "\n"
            + mode
            + self.refr_index.command
        )
        if self.kappa is not None:
            string += f"mode {num} kappa {self.kappa}"

        if self.density is not None:
            string += f"mode {num} kappa {self.density}"

        return string
