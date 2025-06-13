#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import sqlite3
import warnings
import pandas as pd
from tqdm import tqdm
from os.path import splitext

from .data import check_input
from .data import create_crosslink
from .data import create_csm
from .data import create_parser_result
from .constants import MODIFICATIONS
from .parser_util import format_sequence
from .parser_util import get_bool_from_value

from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __csms_to_msannika(
    data: List[Dict[str, Any]],
    filename: Optional[str],
    format: Literal["csv", "xlsx"],
) -> pd.DataFrame:
    msannika_df = pd.DataFrame({})
    return msannika_df


def __xls_to_msannika(
    data: List[Dict[str, Any]],
    filename: Optional[str],
    format: Literal["csv", "xlsx"],
) -> pd.DataFrame:
    msannika_df = pd.DataFrame({})
    return msannika_df


def to_msannika(
    data: List[Dict[str, Any]],
    filename: Optional[str] = None,
    format: Literal["csv", "xlsx"] = "csv",
) -> pd.DataFrame:
    msannika_df = pd.DataFrame({})
    return msannika_df
