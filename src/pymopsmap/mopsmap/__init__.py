from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from pymopsmap.classes import (
        MicroParameters,
        MicroParametersDispatch,
        OptiProps,
    )

    Modes = MicroParameters | list[MicroParameters]


def compute_optical_properties(
    dispatch: Modes | MicroParametersDispatch,
) -> OptiProps:
    """
    Compute optical properties for a dispatch of microphysical parameters.
    Concatenates all dispatched optical properties in a single OptiProps
    instance.

    Parameters
    ----------
    dispatch : Modes | MicroParametersDispatch
        Either a dispatch of microphysical parameters of a list of modes
        of microphysical parameters.

    Returns
    -------
    OptiProps
        The resulting optical properties.
    """
    from pymopsmap.classes import MicroParametersDispatch, extend_optiprops

    if not isinstance(dispatch, MicroParametersDispatch):
        return _compute_optical_properties_single_dispatch(dispatch)

    return extend_optiprops(
        index=dispatch.params,
        optiprops_li=[
            _compute_optical_properties_single_dispatch(mode)
            for mode in dispatch
        ],
    )


def _compute_optical_properties_single_dispatch(modes: Modes) -> OptiProps:
    """
    Compute optical properties for a single dispatch of microphysical
    parameters.
    """
    from .launch_file_format import write_launching_file
    from .launcher import launch_mopsmap
    from .output_format import format_mopsmap_outputs

    files = write_launching_file(mp=modes)
    dico = launch_mopsmap(input_filename=files["mopsmap"])
    op = format_mopsmap_outputs(dico)
    return op


__all__ = ["compute_optical_properties"]
