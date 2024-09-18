#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pandas as pd
from os.path import splitext

#from .data import create_crosslink
#from .data import create_csm
#from .transform import to_dataframe

from typing import BinaryIO
from typing import Union
#from typing import Any

def read_custom():
    return


def read_msannika(input: Union[str, BinaryIO], format: str = "auto", sep: str = "\t"):
    _data = None
    if format == "auto" and type(input) is not str:
        raise ValueError("Can't detect format for file-like objects. Please specify format manually!")
    if format == "auto":
        file_extension = splitext(input)
        if file_extension == ".tsv" or file_extension == ".csv":
            _data = pd.read_csv(input, sep = sep)
        elif file_extension == ".xlsx":
            _data = pd.read_excel(input, engine = "openpyxl")
        else:
            raise ValueError("")
    return


def read():
    return
