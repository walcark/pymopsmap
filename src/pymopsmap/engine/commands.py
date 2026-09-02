"""MOPSMAP input file command builders."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pymopsmap.models import PSD, MicroParameters, Shape
from pymopsmap.utils import PosFloat64List, SortedPosFloat64List, get_tempfile

if TYPE_CHECKING:
    from pymopsmap.models.output_request import OutputRequest


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
        return f"wavelength {wavelengths[0]:.6f}"
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
            f.write(f"{w:.6f} {r:.6f} {i:.6f}\n")

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
        return f"refrac {nr[0]:.6f} {ni[0]:.6f}"
    filename = str(write_refr_file(wl=wl, nr=nr, ni=ni, mode_index=mode_index))
    return f"refrac file '{filename}'"


def shape_command(shape: Shape) -> str:
    return shape.command


def psd_command(psd: PSD) -> str:
    return psd.command


def output_request_commands(
    output_types: OutputRequest,
    ascii_base: Path,
    nc_path: Path,
    n_angles: int = 2000,
    rh: float | None = None,
) -> str:
    """Build the output section of a MOPSMAP launch file."""
    from pymopsmap.models.output_request import OutputType

    _ASCII_TYPES = {
        OutputType.PHASE_FUNCTION,
        OutputType.SCATTERING_MATRIX,
        OutputType.VOLUME_SCATTERING_FUNCTION,
        OutputType.LIDAR,
        OutputType.COEFF,
    }
    lines = [f"output num_theta {n_angles}", wl_command()]
    if rh is not None:
        lines.append(f"rH {rh}")
    lines.append("output integrated")
    lines.append(f"output netcdf '{nc_path}'")
    ascii_needed = output_types & _ASCII_TYPES
    if ascii_needed:
        lines.append(f"output ascii_file '{ascii_base}'")
        for otype in sorted(ascii_needed, key=lambda x: x.value):
            lines.append(f"output {otype.value}")
    return "\n".join(lines)


if __name__ == "__main__":
    from .microparams import FixedPSD, Sphere, Spheroid

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
