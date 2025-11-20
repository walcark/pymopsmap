from pydantic import BaseModel, field_validator

from .refractive_index import RefractiveIndex
from .shape import Shape
from .psd import PSD


class Mode(BaseModel):
    shape: Shape
    psd: PSD
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
