import numpy as np

from pymopsmap.adapt import cams_to_smartg, CamsVersion, CamsAerosol
from pymopsmap.utils import DATA_PATH


def main():
    version = CamsVersion.V49_R1
    wls = np.linspace(0.330, 2.1, 50)
    rhs = np.linspace(0.0, 95.0, 10)
    for specie in CamsAerosol:
        cams_to_smartg(specie, version, wls, rhs, DATA_PATH / "smartg")


if __name__ == "__main__":
    main()
