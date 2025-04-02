#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from tqdm import tqdm

from .data import check_input
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


def read_maxquant(
    files: str | List[str] | BinaryIO,
    decoy_prefix: str = "REV__",
    modifications: Dict[str, float] = MODIFICATIONS,
    sep: str = "\t",
) -> Dict[str, Any]:
    """Read a MaxQuant result file.

    Reads a MaxQuant crosslink-spectrum-matches result file "crosslinkMsms.txt" in ``.txt`` (tab delimited) format
    and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MaxQuant result file(s) or a file-like object/stream.
    decoy_prefix : str, default = "REV__"
        The prefix that indicates that a protein is from the decoy database.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
    sep : str, default = "\t"
        Seperator used in the ``.txt`` file.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches.
    KeyError
        If one of the found post-translational-modifications could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_maxquant
    >>> csms_from_xlsx = read_maxquant("data/maxquant/run1/crosslinkMsms.txt")
    """
    ## check input
    _ok = check_input(decoy_prefix, "decoy_prefix", str)
    _ok = check_input(modifications, "modifications", dict, float)
    _ok = check_input(sep, "sep", str)

    ## data structures
    csms = list()

    ## handle input
    if not isinstance(files, list):
        inputs = [files]
    else:
        inputs = files

    ## process data
    for input in inputs:
        data = pd.read_csv(input, sep=sep, low_memory=False)
        xl = data.dropna(axis=0, subset=["Proteins2"])
        for i, row in tqdm(xl.iterrows(), total=xl.shape[0], desc="Reading MaxQuant CSMs..."):
            # create csm
            csm = create_csm()
            csms.append(csm)
    ## check results
    if len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="MaxQuant",
        csms=csms,
        crosslinks=None,
    )


def read_maxlynx(
    files: str | List[str] | BinaryIO,
    decoy_prefix: str = "REV__",
    modifications: Dict[str, float] = MODIFICATIONS,
    sep: str = "\t",
) -> Dict[str, Any]:
    """Read a MaxLynx result file.

    Reads a MaxLynx crosslink-spectrum-matches result file "crosslinkMsms.txt" in ``.txt`` (tab delimited) format
    and returns a ``parser_result``. This is an alias for the MaxQuant reader.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MaxLynx result file(s) or a file-like object/stream.
    decoy_prefix : str, default = "REV__"
        The prefix that indicates that a protein is from the decoy database.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
    sep : str, default = "\t"
        Seperator used in the ``.txt`` file.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches.
    KeyError
        If one of the found post-translational-modifications could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_maxlynx
    >>> csms_from_xlsx = read_maxlynx("data/maxquant/run1/crosslinkMsms.txt")
    """
    return read_maxquant(files, decoy_prefix, modifications, sep)
