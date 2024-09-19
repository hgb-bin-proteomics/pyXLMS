#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from os.path import splitext

from .data import create_crosslink
from .data import create_csm

from typing import BinaryIO
from typing import List
from typing import Dict
from typing import Literal
from typing import Any

# todo: move to data
def create_parser_result(search_engine: str, csms: List[Dict[str, any]] | None, crosslinks: List[Dict[str, Any]] | None) -> Dict[str, Any]:
    """Creates a parser result data structure.

    Contains all necessary data elements that should be contained in a result returned by a crosslink search engine result parser.

    Parameters
    ----------
    search_engine : str
        Name of the identifying crosslink search engine.
    csms : list of dict, or None
        List of crosslink-spectrum-matches as created by ``data.create_csm()``.
    crosslinks : list of dict, or None
        List of crosslinks as created by ``data.create_crosslink()``.

    Returns
    -------
    dict
        The parser result data structure which is a dictionary with keys ``data_type``, ``search_engine``, ``crosslink-spectrum-matches`` and
        ``crosslinks``.

    Examples
    --------
    >>> from pyXLMS import parser.create_parser_result
    >>> result = create_parser_result("MS Annika", None, None)
    >>> result["data_type"]
    'parser_result'
    >>> result["search_engine"]
    'MS Annika'
    """
    return {"data_type": "parser_result",
            "search_engine": search_engine,
            "crosslink-spectrum-matches": csms,
            "crosslinks": crosslinks}


def read_custom():
    return


def read_msannika(input: str | BinaryIO, format: Literal["auto", "csv", "tsv", "xlsx"] = "auto", sep: str = "\t") -> Dict[str, Any]:
    """Read an MS Annika result file.

    Reads an MS Annika crosslink-spectrum-matches result file or crosslink result file in ``.csv`` or ``.xlsx`` format and returns a
    ``parser_result``.

    Parameters
    ----------
    input : str or file stream
        The name/path of the MS Annika result file or a file-like object/stream.
    format : "auto", "csv", "tsv", or "xlsx", default = "auto"
        The format of the result file. ``"auto"`` is only available if the name/path to the MS Annika result file is given.
    """
    ## reading data
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
            raise ValueError(f"Detected file extension {file_extension} is not supported! Input file has to be a valid file with extension '.csv', '.tsv' or '.xlsx'!")
    elif format in ["csv", "tsv", "xlsx"]:
        if format == "xlsx":
            data = pd.read_excel(input, engine = "openpyxl")
        else:
            data = pd.read_csv(input, sep = sep)
    else:
        raise ValueError(f"Provided input format {format} is not supported! Input format has to be of type 'csv', 'tsv' or 'xlsx'!")
    if data is None:
        raise RuntimeError("Something went wrong while reading the file! Please file a bug report!")
    ## detect input file type
    col_names = data.columns.values.tolist()
    is_crosslink_dataframe = "# CSMs" in col_names
    ## process data
    crosslinks = list()
    csms = list()
    if is_crosslink_dataframe:
        for i, row in data.iterrows():
            # create crosslink
            crosslink = create_crosslink(peptide_a = "",
                                         xl_position_peptide_a = 0,
                                         proteins_a = [""],
                                         xl_position_proteins_a = [0],
                                         decoy_a = False,
                                         peptide_b = "",
                                         xl_position_peptide_b = 0,
                                         proteins_b = [""],
                                         xl_position_proteins_b = [0],
                                         decoy_b = False,
                                         score = 0.0)
            crosslinks.append(crosslink)
    else:
        for i, row in data.iterrows():
            # create csm
            csm = create_csm(peptide_a = "",
                             modifications_a = {},
                             xl_position_peptide_a = 0,
                             proteins_a = [""],
                             xl_position_proteins_a = [0],
                             pep_position_proteins_a = [0],
                             score_a = 0.0,
                             decoy_a = False,
                             peptide_b = "",
                             modifications_b = {},
                             xl_position_peptide_b = 0,
                             proteins_b = [""],
                             xl_position_proteins_b = [0],
                             pep_position_proteins_b = [0],
                             score_b = 0.0,
                             decoy_b = False,
                             score = 0.0,
                             spectrum_file = "",
                             scan_nr = 0,
                             charge = 0,
                             rt = 0.0,
                             im_cv = None)
            csms.append(csm)
    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError("No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!")
    ## return parser result
    return create_parser_result(search_engine = "MS Annika",
                                csms = csms if len(csms) > 0 else None,
                                crosslinks = crosslinks if len(crosslinks) > 0 else None)


def read():
    return
