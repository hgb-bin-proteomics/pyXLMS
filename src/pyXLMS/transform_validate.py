#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .data import check_input_multi
from .data import create_parser_result
from .transform_util import get_available_keys

from typing import Optional
from typing import Dict
from typing import List
from typing import Tuple
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def validate(
    data: List[Dict[str, Any]] | Dict[str, Any],
    fdr: float = 0.01,
    formula: Literal["D/T", "(TD+DD)/TT", "(TD-DD)/TT"] = "D/T",
    separate_intra_inter: bool = False,
    ignore_missing_labels: bool = False,
) -> List[Dict[str, Any]] | Dict[str, Any]:
    return
