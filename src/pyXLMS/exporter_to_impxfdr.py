#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .exporter_to_msannika import to_msannika

from typing import Optional
from typing import Dict
from typing import Any
from typing import List


def to_impxfdr(
    data: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    # assert csms
    return to_msannika(data, filename, format="xlsx")
