
from pydantic import BaseModel

from .refractive_index import RefractiveIndex
from .shape import Shape
from .psd import Size

class Mode(BaseModel):
    shape: Shape
    size: Size
    refr_index: RefractiveIndex
    kappa: float | None = None
    density: float | None = None

    def to_section(self, num: int | None = None) -> str:
        """
        Convert a Mode to MOPSMAP section text.
        """
        string = (
            self.shape.to_section(num) + "\n" +
            self.size.to_section(num) + "\n" +
            self.refr_index.to_section(num)
        )
        if self.kappa is not None:
            string += f"mode {num} kappa {self.kappa}"

        if self.density is not None:
            string += f"mode {num} kappa {self.density}"

        return string
        
