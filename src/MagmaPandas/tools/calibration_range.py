import warnings as w
from typing import Dict

import numpy as np
import pandas as pd

from MagmaPandas.core.magma_protocol import Magma


def _check_calibration_range_Series(melt: pd.Series, calibration_range: Dict) -> None:

    checks = np.array([], dtype=bool)
    for element, *range in calibration_range:
        row = [element] if isinstance(element, str) else element
        np.append(checks, not (range[0] < melt[row].sum() < range[1]))

    if sum(checks) < 1:
        return

    w.warn(f"sample has composition outside the thermometer calibration range")


def _check_calibration_range_dataframe(melt: Magma, calibration_range: Dict) -> None:

    checks = pd.DataFrame(dtype=bool)

    for element, *range in calibration_range:
        col = [element] if isinstance(element, str) else element
        checks[str(element)] = ~melt[col].sum(axis=1).between(range[0], range[1])

    outside_range = checks.any(axis=1)

    if outside_range.sum() < 1:
        return

    w.warn(
        f"samples {*melt.index[outside_range],} have compositions outside the thermometer calibration range"
    )


def _check_calibration_range(melt: pd.Series | Magma, calibration_range) -> None:

    if isinstance(melt, pd.Series):
        _check_calibration_range_Series(melt, calibration_range)
    elif isinstance(melt, pd.DataFrame):
        _check_calibration_range_dataframe(melt, calibration_range)
