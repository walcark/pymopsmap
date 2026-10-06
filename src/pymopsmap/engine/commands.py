"""MOPSMAP input file command builders."""

from __future__ import annotations

import hashlib
from pathlib import Path

from pymopsmap.engine.workspace import Workspace
from pymopsmap.microparams import MicroParameters
from pymopsmap.utils import PosFloat64List, SortedPosFloat64List

# Scientific notation, matching the format MOPSMAP uses in its own data files
# and reads with list-directed input. Fixed point with six decimals wrote every
# imaginary index below 5e-7 as zero, and collapsed nearby wavelengths onto the
# same value, which MOPSMAP rejects as a non-ascending grid.
_FLOAT = ".10e"


def microparams_command(
    modes: list[MicroParameters], workspace: Workspace
) -> str:
    """The mode block of the launch file, one stanza per mode."""
    return "\n".join(
        _single_microparams_command(mode, workspace, index)
        for index, mode in enumerate(modes, start=1)
    )


def _single_microparams_command(
    mp: MicroParameters, workspace: Workspace, num: int = 1
) -> str:
    """Returns the full command of a MicroParams instance."""
    mode: str = f"mode {num} "
    string = (
        mode
        + mp.shape.command
        + "\n"
        + mode
        + mp.psd.command
        + "\n"
        + mode
        + refr_command(
            workspace=workspace,
            wl=mp.wavelength,
            nr=mp.n_real,  # type: ignore[arg-type]
            ni=mp.n_imag,  # type: ignore[arg-type]
        )
    )
    if mp.nonabs_fraction:
        string += "\n" + mode + f"refrac nonabs_fraction {mp.nonabs_fraction}"

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
    workspace: Workspace,
    wl: SortedPosFloat64List,
    nr: PosFloat64List,
    ni: PosFloat64List,
) -> Path:
    """
    Write the refractive index of one mode, and return the file it went to.

    The file is named after what it holds rather than after the mode: two
    modes with the same refractive index share one file, and two with
    different ones cannot overwrite each other. An ensemble built one measured
    particle at a time runs to thousands of modes on a single spectrum, which
    is a single file here and was a thousand before.
    """
    content = "".join(
        f"{w:{_FLOAT}} {r:{_FLOAT}} {i:{_FLOAT}}\n"
        for w, r, i in zip(wl, nr, ni)
    )
    digest = hashlib.blake2b(content.encode(), digest_size=8).hexdigest()
    filename = workspace.file(f"ri_{digest}.txt")
    if not filename.exists():
        filename.write_text(content)
    return filename


def refr_command(
    workspace: Workspace,
    wl: SortedPosFloat64List,
    nr: PosFloat64List,
    ni: PosFloat64List,
) -> str:
    # MOPSMAP bug: interpolate_linear returns weight_upper=weight_lower=1.0 for
    # single-element arrays, doubling the refractive index and causing an
    # out-of-range error. Use constant refrac command to bypass the file path.
    if len(wl) == 1:
        return f"refrac {nr[0]:{_FLOAT}} {ni[0]:{_FLOAT}}"
    filename = str(write_refr_file(workspace, wl=wl, nr=nr, ni=ni))
    return f"refrac file '{filename}'"
