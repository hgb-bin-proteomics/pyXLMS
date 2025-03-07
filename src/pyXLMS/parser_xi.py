#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from os.path import splitext

from .data import check_input
from .data import create_crosslink
from .data import create_csm
from .data import create_parser_result
from .constants import XI_MODIFICATION_MAPPING
from .parser_util import format_sequence
from .parser_util import get_bool_from_value

from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List


def detect_xi_filetype():
    return

def read_xisearch(data: pd.DataFrame, modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING) -> List[Dict[str, Any]]:
    return

def read_xifdr_csms(data: pd.DataFrame, modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING) -> List[Dict[str, Any]]:
    return

def read_xifdr_crosslinks(data: pd.DataFrame, modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING) -> List[Dict[str, Any]]:
    return

def read_xi(
    files: str | List[str] | BinaryIO,
    modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING
) -> Dict[str, Any]:
    """Read an MS Annika result file.

    Reads an MS Annika crosslink-spectrum-matches result file or crosslink result file in ``.csv`` or ``.xlsx`` format and returns a
    ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MS Annika result file(s) or a file-like object/stream.
    modifications: dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of modification names to modification masses.
    format : "auto", "csv", "tsv", "txt", or "xlsx", default = "auto"
        The format of the result file. ``"auto"`` is only available if the name/path to the MS Annika result file is given.
    sep : str, default = "\t"
        Seperator used in the ``.csv`` or ``.tsv`` file. Parameter is ignored if the file is in ``.xlsx`` format.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    ValueError
        If the input format is not supported or cannot be inferred.
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslinks or crosslink-spectrum-matches.
    KeyError
        If one of the found post-translational-modifications could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_msannika
    >>> csms_from_xlsx = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")

    >>> from pyXLMS.parser import read_msannika
    >>> crosslinks_from_xlsx = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx")

    >>> from pyXLMS.parser import read_msannika
    >>> csms_from_tsv = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt")

    >>> from pyXLMS.parser import read_msannika
    >>> crosslinks_from_tsv = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt")
    """
    ## check input
    _ok = check_input(modifications, "modifications", dict, tuple)

    ## helper functions


    ## data structures
    crosslinks = list()
    csms = list()

    ## handle input
    if not isinstance(files, list):
        inputs = [files]
    else:
        inputs = files

    for input in inputs:
        ## reading data
        data = None
        if format == "auto" and not isinstance(input, str):
            raise ValueError(
                "Can't detect format for file-like objects. Please specify format manually!"
            )
        # and isinstance specified for type checking
        if format == "auto" and isinstance(input, str):
            file_extension = splitext(input)[1]
            if (
                file_extension == ".txt"
                or file_extension == ".tsv"
                or file_extension == ".csv"
            ):
                data = pd.read_csv(input, sep=sep)
            elif file_extension == ".xlsx":
                data = pd.read_excel(input, engine="openpyxl")
            else:
                raise ValueError(
                    f"Detected file extension {file_extension} is not supported! Input file has to be a valid file with extension '.csv', '.tsv' or '.xlsx'!"
                )
        elif format in ["csv", "tsv", "txt", "xlsx"]:
            if format == "xlsx":
                data = pd.read_excel(input, engine="openpyxl")
            else:
                data = pd.read_csv(input, sep=sep)
        else:
            raise ValueError(
                f"Provided input format {format} is not supported! Input format has to be of type 'csv', 'tsv' or 'xlsx'!"
            )
        if data is None:
            raise RuntimeError(
                "Something went wrong while reading the file! Please file a bug report!"
            )
        # this should be impossible, but check here for pyright
        if not isinstance(data, pd.DataFrame):
            raise RuntimeError(
                "Something went wrong while reading the file! Please file a bug report!"
            )
        ## detect input file type

        ## process data

    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="xiSearch/xiFDR",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )
