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

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

def detect_xi_filetype(data: pd.DataFrame) -> Literal["xisearch", "xifdr_csms", "xifdr_crosslinks"]:
    """Detects the source application of the data.

    Detects whether the input data is originating from xiSearch or xiFDR, and if xiFDR which type of data is
    being read (CSMs or crosslinks).

    Parameters
    ----------
    data : pd.DataFrame
        The input data originating from xiSearch or xiFDR.

    Returns
    -------
    str
        "xisearch" if a xiSearch result file was read, "xifdr_csms" if CSMs from xiFDR were read,
        "xifdr_crosslinks" if crosslinks from xiFDR were read.

    Raises
    ------
    ValueError
        If the data source could not be determined.

    Examples
    --------
    >>> from pyXLMS.parser_xi import detect_xi_filetype
    >>> import pandas as pd
    >>> df1 = pd.read_csv("data/xi/r1_Xi1.7.6.7.csv")
    >>> detect_xi_filetype(df1)
    'xisearch'

    >>> from pyXLMS.parser_xi import detect_xi_filetype
    >>> import pandas as pd
    >>> df2 = pd.read_csv("data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv")
    >>> detect_xi_filetype(df2)
    'xifdr_csms'

    >>> from pyXLMS.parser_xi import detect_xi_filetype
    >>> import pandas as pd
    >>> df3 = pd.read_csv("data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv")
    >>> detect_xi_filetype(df3)
    'xifdr_crosslinks'
    """
    col_names = data.columns.values.tolist()
    if "AllScore" in col_names:
        return "xisearch"
    if "LinkPos1" in col_names:
        return "xifdr_csms"
    if "ToSite" in col_names:
        return "xifdr_crosslinks"

    raise ValueError("Could not infer data source, are you sure you read a xi result file?")

    return "err"

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
    """Read a xiSearch/xiFDR result file.

    Reads a xiSearch crosslink-spectrum-matches result file or a xiFDR crosslink-spectrum-matches
    result file or crosslink result file in ``.csv`` format and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the xiSearch/xiFDR result file(s) or a file-like object/stream.
    modifications: dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"Ccm"``) to their modifications (e.g. ``("C", "Carbamidomethyl", 57.021464)``).

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) contain no crosslinks or crosslink-spectrum-matches.

    Examples
    --------
    >>> from pyXLMS.parser import read_xi
    >>> csms_from_xiSearch = read_xi("data/xi/r1_Xi1.7.6.7.csv")

    >>> from pyXLMS.parser import read_xi
    >>> csms_from_xiFDR = read_xi("data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv")

    >>> from pyXLMS.parser import read_xi
    >>> crosslinks_from_xiFDR = read_xi("data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv")
    """
    ## check input
    _ok = check_input(modifications, "modifications", dict, tuple)

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
        data = pd.read_csv(input)
        ## detect input file type
        xi_file_type = detect_xi_filetype(data)
        ## process data
        if xi_file_type == "xifdr_csms":
            csms += read_xifdr_csms(data, modifications)
        elif xi_file_type ==  "xifdr_crosslinks":
            crosslinks += read_xifdr_crosslinks(data, modifications)
        else:
            csms += read_xisearch(data, modifications)

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
