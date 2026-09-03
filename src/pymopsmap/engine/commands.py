"""MOPSMAP input file command builders."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pymopsmap.microparams import MicroParameters
from pymopsmap.psd import PSD
from pymopsmap.shapes import Shape
from pymopsmap.utils import PosFloat64List, SortedPosFloat64List, get_tempfile

if TYPE_CHECKING:
    pass

# Scientific notation, matching the format MOPSMAP uses in its own data files
# and reads with list-directed input. Fixed point with six decimals wrote every
# imaginary index below 5e-7 as zero, and collapsed nearby wavelengths onto the
# same value, which MOPSMAP rejects as a non-ascending grid.
_FLOAT = ".10e"


def microparams_command(
    mp: MicroParameters | list[MicroParameters],
) -> str:
    """Returns the full command of a MicroParams instance or list."""
    mpli = [mp] if isinstance(mp, MicroParameters) else mp

    mp_command = "\n".join(
        [
            _single_microparams_command(m, i + 1)
            for m, i in zip(mpli, range(len(mpli)))
        ]
    )

    return mp_command


def _single_microparams_command(mp: MicroParameters, num: int = 1) -> str:
    """Returns the full command of a MicroParams instance."""
    mode: str = f"mode {num} "
    string = (
        mode
        + shape_command(mp.shape)
        + "\n"
        + mode
        + psd_command(mp.psd)
        + "\n"
        + mode
        + refr_command(
            wl=mp.wavelength,
            nr=mp.n_real,  # type: ignore[arg-type]
            ni=mp.n_imag,  # type: ignore[arg-type]
            mode_index=num,
        )
    )
    if mp.kappa is not None:
        string += "\n" + mode + f"kappa {mp.kappa}"

    if mp.density is not None:
        string += "\n" + mode + f"density {mp.density}"

    return string


def wl_command(wavelengths: SortedPosFloat64List | None = None) -> str:
    if wavelengths is not None and len(wavelengths) == 1:
        return f"wavelength {wavelengths[0]:{_FLOAT}}"
    return "wavelength from_refrac_file"


def write_refr_file(
    wl: SortedPosFloat64List,
    nr: PosFloat64List,
    ni: PosFloat64List,
    mode_index: int = 1,
) -> Path:
    """
    Write refractive index and wavelenth in a file. Returns
    the file location.

    The file name carries the mode index: every mode of a mixture needs its
    own file, otherwise the modes overwrite each other and MOPSMAP reads the
    same refractive index for all of them.
    """
    filename = get_tempfile(filename=f"ri_{mode_index}.txt")

    with open(filename, "w") as f:
        for w, r, i in zip(wl, nr, ni):
            f.write(f"{w:{_FLOAT}} {r:{_FLOAT}} {i:{_FLOAT}}\n")

    return filename


def refr_command(
    wl: SortedPosFloat64List,
    nr: PosFloat64List,
    ni: PosFloat64List,
    mode_index: int = 1,
) -> str:
    # MOPSMAP bug: interpolate_linear returns weight_upper=weight_lower=1.0 for
    # single-element arrays, doubling the refractive index and causing an
    # out-of-range error. Use constant refrac command to bypass the file path.
    if len(wl) == 1:
        return f"refrac {nr[0]:{_FLOAT}} {ni[0]:{_FLOAT}}"
    filename = str(write_refr_file(wl=wl, nr=nr, ni=ni, mode_index=mode_index))
    return f"refrac file '{filename}'"


def shape_command(shape: Shape) -> str:
    return shape.command


def psd_command(psd: PSD) -> str:
    return psd.command
