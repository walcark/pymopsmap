"""MicroParameters: one aerosol mode, validated."""

from pydantic import (
    BaseModel,
    NonNegativeFloat,
    PositiveFloat,
    model_validator,
)

from pymopsmap.psd import PSD
from pymopsmap.shapes import Shape
from pymopsmap.utils import (
    Float64List,
    PosFloat64List,
    SortedPosFloat64List,
)


class MicroParameters(BaseModel):
    wavelength: SortedPosFloat64List
    n_real: PosFloat64List | float
    n_imag: Float64List | float
    shape: Shape
    psd: PSD
    kappa: NonNegativeFloat | None = None
    density: PositiveFloat | None = None

    @model_validator(mode="after")
    def broadcast_refractive_index(self):
        n = len(self.wavelength)
        if isinstance(self.n_real, float):
            self.n_real = [self.n_real] * n
        elif len(self.n_real) == 1 and n > 1:
            self.n_real = self.n_real * n
        elif len(self.n_real) != n:
            raise ValueError(
                f"n_real length ({len(self.n_real)}) must match"
                f" wavelength length ({n})"
            )
        if isinstance(self.n_imag, float):
            self.n_imag = [self.n_imag] * n
        elif len(self.n_imag) == 1 and n > 1:
            self.n_imag = self.n_imag * n
        elif len(self.n_imag) != n:
            raise ValueError(
                f"n_imag length ({len(self.n_imag)}) must match"
                f" wavelength length ({n})"
            )
        return self
