"""
microparams_dispatch.py

Author  : Kévin Walcarius
Date    : 2025-01-08
Version : 1.0
License : MIT
Summary : Class to encapsulate the evolution of microphysical
          parameters of an aerosol specie with external parameters.
"""

from dataclasses import dataclass, field
from typing import Iterable

from .microparams import MicroParameters

MPs = MicroParameters | list[MicroParameters]


@dataclass
class MicroParametersDispatch(Iterable):
    modes: list[list[MicroParameters]] = field(default_factory=list)
    params: list[dict[str, float]] = field(default_factory=list)

    def append(self, modes: MPs, params: dict[str, float]):
        modes = [modes] if isinstance(modes, MicroParameters) else modes
        self.modes.append(modes)
        self.params.append(params)
        assert len(self.modes) == len(self.params)

    def __iter__(self):
        for mode in self.modes:
            yield mode
