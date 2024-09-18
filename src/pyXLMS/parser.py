#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pandas as pd
from os.path import splitext

from .data import create_crosslink
from .data import create_csm
from .transform import to_dataframe

from typing import BinaryIO
from typing import List
from typing import Dict
from typing import Union
from typing import Any

def create_parser_result(search_engine: str, csms: List[Dict[str, any]], crosslinks: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Creates a parser result data structure.

    Contains all necessary data elements that should be contained in a result returned by a crosslink search engine result parser.

    Parameters
    ----------
    search_engine : str
        Name of the identifying crosslink search engine.
    csms : list of dict
        List of crosslink-spectrum-matches as created by `data.create_csm`.
    crosslinks : list of dict
        List of crosslinks as created by `data.create_crosslink`.

    Returns
    -------
    dict
        The parser result data structure which is a dictionary with keys data_type, search_engine, crosslink-spectrum-matches and
        crosslinks.
    """
    return {"data_type": "parser_result",
            "search_engine": search_engine,
            "crosslink-spectrum-matches": csms,
            "crosslinks": crosslinks}


def read_custom():
    return


def read_msannika(input: Union[str, BinaryIO], format: str = "auto", sep: str = "\t") -> Dict[str, Any]:
    data = None
    if format == "auto" and type(input) is not str:
        raise ValueError("Can't detect format for file-like objects. Please specify format manually!")
    if format == "auto":
        file_extension = splitext(input)
        if file_extension == ".tsv" or file_extension == ".csv":
            data = pd.read_csv(input, sep = sep)
        elif file_extension == ".xlsx":
            data = pd.read_excel(input, engine = "openpyxl")
        else:
            raise ValueError("")
    return


def read():
    return
