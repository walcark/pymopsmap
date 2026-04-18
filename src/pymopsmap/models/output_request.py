"""OutputType enum and OutputRequest frozen-set type."""

from enum import StrEnum


class OutputType(StrEnum):
    INTEGRATED = "integrated"
    LIDAR = "lidar"
    PHASE_FUNCTION = "phase_function"
    SCATTERING_MATRIX = "scattering_matrix"
    VOLUME_SCATTERING_FUNCTION = "volume_scattering_function"
    COEFF = "coeff"


OutputRequest = frozenset[OutputType]
DEFAULT_OUTPUT: OutputRequest = frozenset({OutputType.INTEGRATED})
