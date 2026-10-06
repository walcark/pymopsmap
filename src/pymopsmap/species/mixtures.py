"""The OPAC climatological mixtures."""

from __future__ import annotations

from .catalog import OpacSpecie
from .mix import Mix

# Hess, Koepke and Schult (1998), BAMS 79:831, as published on
# https://geisa.aeris-data.fr/opac/ . Number concentrations in cm-3; the
# published per-type total is kept beside each entry, and a test checks that
# the components sum to it.
OPAC_COMPOSITIONS: dict[str, dict[OpacSpecie, float]] = {
    "continental_clean": {
        OpacSpecie.WASO: 2600.0,
        OpacSpecie.INSO: 0.15,
    },
    "continental_average": {
        OpacSpecie.WASO: 7000.0,
        OpacSpecie.INSO: 0.4,
        OpacSpecie.SOOT: 8300.0,
    },
    "continental_polluted": {
        OpacSpecie.WASO: 15700.0,
        OpacSpecie.INSO: 0.6,
        OpacSpecie.SOOT: 34300.0,
    },
    "urban": {
        OpacSpecie.WASO: 28000.0,
        OpacSpecie.INSO: 1.5,
        OpacSpecie.SOOT: 130000.0,
    },
    "desert": {
        OpacSpecie.WASO: 2000.0,
        OpacSpecie.MINM: 269.5,
        OpacSpecie.MIAM: 30.5,
        OpacSpecie.MICM: 0.142,
    },
    "maritime_clean": {
        OpacSpecie.WASO: 1500.0,
        OpacSpecie.SSAM: 20.0,
        OpacSpecie.SSCM: 3.2e-3,
    },
    "maritime_polluted": {
        OpacSpecie.WASO: 3800.0,
        OpacSpecie.SSAM: 20.0,
        OpacSpecie.SSCM: 3.2e-3,
        OpacSpecie.SOOT: 5180.0,
    },
    "maritime_tropical": {
        OpacSpecie.WASO: 590.0,
        OpacSpecie.SSAM: 10.0,
        OpacSpecie.SSCM: 1.3e-3,
    },
    "arctic": {
        OpacSpecie.WASO: 1300.0,
        OpacSpecie.INSO: 0.01,
        OpacSpecie.SSAM: 1.9,
        OpacSpecie.SOOT: 5300.0,
    },
    "antarctic": {
        OpacSpecie.SUSO: 42.9,
        OpacSpecie.SSAM: 0.047,
        OpacSpecie.MITR: 0.0053,
    },
}

PER_CM3_TO_PER_M3 = 1e6


def opac_mix(name: str) -> Mix:
    """
    Build one of the OPAC climatological mixtures.

    Parameters
    ----------
    name : str
        A key of ``OPAC_COMPOSITIONS``, for instance
        ``"continental_average"``.

    Returns
    -------
    Mix
        The components at their published concentrations, converted to m-3.
    """
    if name not in OPAC_COMPOSITIONS:
        raise KeyError(
            f"Unknown OPAC mixture '{name}'. "
            f"Available: {sorted(OPAC_COMPOSITIONS)}."
        )
    return Mix(
        {
            specie: value * PER_CM3_TO_PER_M3
            for specie, value in OPAC_COMPOSITIONS[name].items()
        }
    )
