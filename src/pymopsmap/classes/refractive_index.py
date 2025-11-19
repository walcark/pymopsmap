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
from typing import List, Any
from pydantic import BaseModel, field_validator, model_validator
from pathlib import Path
import tempfile, uuid


class RefractiveIndex(BaseModel):
    """
    Unified representation of refractive index for MOPSMAP.

    Accepted forms
    --------------
    1) Constant:
        RefractiveIndex(n_real=1.5, n_imag=0.01)

    2) Wavelength-dependent lists:
        RefractiveIndex(
            wl=[0.44, 0.55],
            n_real=[1.50, 1.48],
            n_imag=[0.01, 0.015]
        )

    3) File-based:
        RefractiveIndex(filename="my_ri.txt")

    Internally:
        - scalar inputs are auto-wrapped into lists
        - wavelength length == n_real length == n_imag length
        - all wavelengths must be sorted and > 0
        - all imaginary parts must be >= 0
    """
    wl: List[float] | None = None
    n_real: List[float] | float | None = None
    n_imag: List[float] | float | None = None
    filename: str | None = None

    @field_validator("wl", "n_real", "n_imag", mode="before")
    @classmethod
    def auto_wrap_lists(cls, v: Any):
        if v is None:
            return None
        if isinstance(v, (float, int)):
            return [float(v)]
        if isinstance(v, list):
            return v
        raise TypeError(f"Invalid type {type(v)}: expected float or list.")

    @model_validator(mode="after")
    def validate_logic(self):
        # Case 3: file-based
        if self.filename is not None:
            return self

        # Case 1 & 2: lists or scalars → now always lists
        if self.wl is None:
            if (isinstance(self.n_real, list) or isinstance(self.n_imag, list)):
                raise ValueError("If wl is None, n_real/n_imag must be scalar.")
            return self

        # Case 2: wavelength-dependent
        Lw, Lr, Li = len(self.wl), len(self.n_real),  len(self.n_imag)
        if not (Lw == Lr == Li):
            raise ValueError(f"Sizes must match: wl({Lw}) n_real({Lr}) n_imag({Li})")
        if any(w <= 0 for w in self.wl):
            raise ValueError("Wavelengths must be > 0")    
        if self.wl != sorted(self.wl):
            raise ValueError("Wavelengths must be sorted ascending.")
        if any(ni < 0 for ni in self.n_imag):
            raise ValueError("Imaginary part must be >= 0")
        return self

    def _write_temp_file(self) -> str:
        fname = f"ri_{uuid.uuid4().hex}.txt"
        path = Path(tempfile.gettempdir()) / fname

        with open(path, "w") as f:
            for wl, nr, ni in zip(self.wl, self.n_real, self.n_imag):
                f.write(f"{wl:.6f} {nr:.6f} {ni:.6f}\n")

        return str(path)

    def to_section(self, num: int | None = None) -> str:
        # Case 3: file-based
        if self.filename:
            return f"mode {num} refrac file '{self.filename}'"

        # Case 1: constant
        if self.wl is None:
            return f"mode {num} refrac constant {self.n_real} {self.n_imag}"

        # Case 2: wavelength-dependent
        tmpfile = self._write_temp_file()
        return f"mode {num} refrac file '{tmpfile}'"