#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pandas as pd
from .data import check_input

from typing import List
from typing import Dict
from typing import Any

def __crosslinks_to_dataframe(data: List[Dict[str, Any]]) -> pd.DataFrame:
    return

def __csms_to_dataframe(data: List[Dict[str, Any]]) -> pd.DataFrame:
    return

def to_dataframe(data: List[Dict[str, Any]]) -> pd.DataFrame:
    """Returns a pandas DataFrame of the given crosslinks or crosslink-spectrum-matches.

    Parameters
    ----------
    data : list
        A list of crosslinks or crosslink-spectrum-matches.

    Returns
    -------
    pandas.DataFrame
        The pandas DataFrame created from the list of input crosslinks or crosslink-spectrum-matches.

    Raises
    ------
    TypeError
        If the list does not contain crosslinks or crosslink-spectrum-matches.
    ValueError
        If the list does not contain any objects.
    """
    ## input checks
    check_input(data, "data", list, dict)
    ## function calls
    if len(data) > 0:
        if "data_type" in data[0] and data[0]["data_type"] == "crosslink":
            return __crosslinks_to_dataframe(data)
        elif "data_type" in data[0] and data[0]["data_type"] == "crosslink-spectrum-match":
            return __csms_to_dataframe(data)
        else:
            raise TypeError("The given data object is not supported!")
    else:
        raise ValueError("Parameter data has to be at least of length one!")
