from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pymopsmap.classes import MicroParameters, OptiProps


def compute_optical_properties(
    mp: MicroParameters | list[MicroParameters],
) -> OptiProps:
    """
    Compute optical properties from one or many MicroParameters.
    """
    from .launch_file_format import write_launching_file
    from .launcher import launch_mopsmap
    from .output_format import format_mopsmap_outputs

    files = write_launching_file(mp=mp)
    dico = launch_mopsmap(input_filename=files["mopsmap"])
    op = format_mopsmap_outputs(dico)
    return op


__all__ = ["compute_optical_properties"]
