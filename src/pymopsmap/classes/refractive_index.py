"""
refractive_index.py

Author  : Kévin Walcarius
Date    : 2025-11-19
Version : 1.0
License : MIT
Summary : Defines the Refractive Index to be used in Mopsmap. Two
          types of refractive index are accepted: constant and
          wavelength-dependant.
"""

from pydantic import BaseModel, model_validator

from pymopsmap.utils import PosFloat64List, SortedPosFloat64List, get_tempfile


class RefractiveIndex(BaseModel):
    """
    Simple RefractiveIndex handler for MOPSMAP.
    """

    wl: SortedPosFloat64List | None = None
    n_real: PosFloat64List | None = None
    n_imag: PosFloat64List | None = None
    filename: str | None = None

    @model_validator(mode="after")
    def validate_logic(self):
        # No need for wl, n_real and n_imag
        if self.filename is not None:
            return self

        # Test shapes of wl, n_real and n_imag
        if any(v is None for v in [self.wl, self.n_real, self.n_imag]):
            raise ValueError("wl, n_real and n_imag should be specified.")

        Lw, Lr, Li = len(self.wl), len(self.n_real), len(self.n_imag)
        if not (Lw == Lr == Li):
            raise ValueError(
                f"Sizes must match: wl({Lw}) n_real({Lr}) n_imag({Li})"
            )
        self._write_temp_file()

        return self

    def _write_temp_file(self) -> str:
        self.filename = get_tempfile(filename="ri.txt")
        with open(self.filename, "w") as f:
            for wl, nr, ni in zip(self.wl, self.n_real, self.n_imag):
                f.write(f"{wl:.6f} {nr:.6f} {ni:.6f}\n")

    @property
    def command(self) -> str:
        return f"refrac file '{self.filename}'"
