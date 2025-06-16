#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
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
    _ok = check_input(targets_only, "targets_only", bool)
    if targets_only:
        data = filter_target_decoy(data)["Target-Target"]
    return to_msannika(data, filename, format="xlsx")
