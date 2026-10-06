"""MicroParameters: one aerosol mode, validated."""

from typing import Literal

from pydantic import (
    BaseModel,
    Field,
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

# How MOPSMAP reads a size given for a nonspherical particle (read_input.f90,
# keyword size_equ): as the radius of the sphere of equal cross section, of
# equal volume, or of equal volume-to-cross-section ratio. Section 2.1 of
# Gasteiger and Wiegner (2018) defines the three, and its Table 5 shows that
# the choice moves the mass-to-backscatter factor by a factor of two.
SizeEquivalence = Literal["cs", "vol", "vol_cs_ratio"]


class MicroParameters(BaseModel):
    wavelength: SortedPosFloat64List
    n_real: PosFloat64List | float
    n_imag: Float64List | float
    shape: Shape
    psd: PSD
    kappa: NonNegativeFloat | None = None
    density: PositiveFloat | None = None
    # The fraction of the mode made of particles that do not absorb. MOPSMAP
    # splits the mode in two, giving the rest an imaginary index divided by
    # 1 - X so the average is unchanged (init_wavelength_refr.f90, line 223).
    # Section 3.1 of Gasteiger and Wiegner (2018) introduces it.
    nonabs_fraction: float = Field(default=0.0, ge=0.0, lt=1.0)
    size_equ: SizeEquivalence = "cs"

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
