"""MicroParametersDispatch: batch aerosol configurations over external parameters."""

from collections.abc import Iterable
from dataclasses import dataclass, field

from .microparams import MicroParameters

MPs = MicroParameters | list[MicroParameters]


@dataclass
class MicroParametersDispatch(Iterable):
    modes: list[list[MicroParameters]] = field(default_factory=list)
    params: list[dict[str, float]] = field(default_factory=list)

    def append(self, modes: MPs, params: dict[str, float]):
        modes = [modes] if isinstance(modes, MicroParameters) else modes
        if self.modes:
            ref_wl = self.modes[0][0].wavelength
            for mp in modes:
                if list(mp.wavelength) != list(ref_wl):
                    raise ValueError(
                        "All MicroParameters in a dispatch must share the same "
                        f"wavelength array. Got {mp.wavelength} vs {ref_wl}."
                    )
        self.modes.append(modes)
        self.params.append(params)
        assert len(self.modes) == len(self.params)

    def __iter__(self):
        yield from self.modes
