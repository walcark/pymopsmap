import numpy as np

from pymopsmap.adapters import (
    CamsAerosol,
    CamsVersion,
    cams_to_kext,
)


def main():
    import matplotlib.pyplot as plt

    version = CamsVersion.V49_R1
    wls = np.linspace(0.330, 2.1, 50)
    rhs = np.linspace(0.0, 95.0, 10)
    for specie in CamsAerosol:
        df = cams_to_kext(specie, version, wls, rhs)
        print(df)
        df.plot()
        plt.show()


if __name__ == "__main__":
    main()
