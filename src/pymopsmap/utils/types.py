from typing import Annotated

import numpy as np
from pydantic import AfterValidator, BeforeValidator


# Methods
def coerce_as_list(value: float | list | np.ndarray):
    value_arr = np.atleast_1d(np.asarray(value))
    return value_arr.tolist()


def assert_strictly_positive(value: list[float]) -> list[float]:
    if any(v <= 0 for v in value):
        raise ValueError("Input value should be strictly positive.")
    return value


def assert_sorted(value: list[float]) -> list[float]:
    value_sorted = sorted(value)
    if any(v != vs for v, vs in zip(value, value_sorted)):
        raise ValueError("Input list should be sorted.")
    return value


# Types
type Float64List = Annotated[list[float], BeforeValidator(coerce_as_list)]

type PosFloat64List = Annotated[
    Float64List, AfterValidator(assert_strictly_positive)
]

type SortedFloat64List = Annotated[Float64List, AfterValidator(assert_sorted)]

type SortedPosFloat64List = Annotated[
    SortedFloat64List, AfterValidator(assert_strictly_positive)
]
