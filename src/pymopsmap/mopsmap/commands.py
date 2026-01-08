"""
mospmap_commands.py

Author  : Kévin Walcarius
Date    : 2025-11-19
Version : 1.0
License : MIT
Summary : Module used to define all the Mopsmap commands, as defined
          in the guide: https://mopsmap.net/mopsmap_userguide.pdf.
"""

from pathlib import Path
from pymopsmap.utils import SortedPosFloat64List, PosFloat64List, get_tempfile
from pymopsmap.classes import MicroParameters, Shape, PSD


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
        + refr_command(wl=mp.wavelength, nr=mp.n_real, ni=mp.n_imag)
    )
    if mp.kappa is not None:
        string += "\n" + mode + f"kappa {mp.kappa}"

    if mp.density is not None:
        string += "\n" + mode + f"density {mp.density}"

    return string


def wl_command() -> str:
    return "wavelength from_refrac_file"


def write_refr_file(
    wl: SortedPosFloat64List, nr: PosFloat64List, ni: PosFloat64List
) -> Path:
    """
    Write refractive index and wavelenth in a file. Returns
    the file location.
    """
    filename = get_tempfile(filename="ri.txt")

    with open(filename, "w") as f:
        for w, r, i in zip(wl, nr, ni):
            f.write(f"{w:.6f} {r:.6f} {i:.6f}\n")

    return filename


def refr_command(
    wl: SortedPosFloat64List, nr: PosFloat64List, ni: PosFloat64List
) -> str:
    filename = str(write_refr_file(wl=wl, nr=nr, ni=ni))
    return f"refrac file '{filename}'"


def shape_command(shape: Shape) -> str:
    return shape.command


def psd_command(psd: PSD) -> str:
    return psd.command


if __name__ == "__main__":
    from .microparams import Sphere, Spheroid, FixedPSD

    mp: MicroParameters = MicroParameters(
        wavelength=[1.5],
        n_real=[1.0],
        n_imag=[1e-4],
        shape=Sphere(),
        psd=FixedPSD(radius=1.0, n=1.0),
    )
    print(microparams_command(mp))
    print(
        refr_command(
            wl=[1.0, 3.0, 5.0], nr=[1.0, 1.0, 1.0], ni=[0.0001, 0.0001, 0.0001]
        )
    )
    print(shape_command(Sphere()))
    print(shape_command(Spheroid(aspect_ratio=1.1, mode="oblate")))
