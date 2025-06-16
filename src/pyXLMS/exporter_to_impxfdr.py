#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .data import create_crosslink_from_csm
from .transform_util import get_available_keys
from .transform_filter import filter_target_decoy
from .exporter_to_msannika import to_msannika

from typing import Optional
from typing import Dict
from typing import Any
from typing import List


def to_impxfdr(
    data: List[Dict[str, Any]],
    filename: Optional[str],
    targets_only: bool = True,
) -> pd.DataFrame:
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    _ok = check_input(targets_only, "targets_only", bool)
    if targets_only:
        data = filter_target_decoy(data)["Target-Target"]
    if "data_type" not in data[0] or data[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
        )
    available_keys = get_available_keys(data)
    if (
        not available_keys["alpha_proteins"]
        or not available_keys["alpha_proteins_crosslink_positions"]
        or not available_keys["beta_proteins"]
        or not available_keys["beta_proteins_crosslink_positions"]
    ):
        raise RuntimeError(
            "Can't export to IMP-X-FDR because not all necessary information is available!"
        )
    if data[0]["data_type"] == "crosslink":
        return to_msannika(data, filename, format="xlsx")
    return to_msannika(
        [create_crosslink_from_csm(csm) for csm in data], filename, format="xlsx"
    )
