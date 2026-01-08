import numpy as np
from pymopsmap.classes import (
    Sphere,
    FixedPSD,
    extend_optiprops,
    MicroParameters,
)


from pymopsmap.mopsmap import compute_optical_properties


def main():
    mp: MicroParameters = MicroParameters(
        wavelength=np.linspace(0.330, 1.5, 100),
        n_real=[1.0] * 100,
        n_imag=[1e-4] * 100,
        shape=Sphere(),
        psd=FixedPSD(radius=1.0, n=1.0),
    )

    op = compute_optical_properties(mp=mp)

    index = [
        {"a": 1, "b": 1},
        {"a": 2, "b": 1},
        {"a": 1, "b": 2},
        {"a": 2, "b": 2},
    ]
    opli = [op] * 4

    ops = extend_optiprops(index, opli)
    print(ops)


if __name__ == "__main__":
    main()
