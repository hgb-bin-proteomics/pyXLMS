#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import sqlite3
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


def detect_scout_filetype(
    data: pd.DataFrame,
) -> Literal["scout_csms_unfiltered", "scout_csms_filtered", "scout_xl"]:
    r"""Detects the Scout-related source of the data.

    Detects whether the input data is unfiltered crosslink-spectrum-matches, filtered crosslink-spectrum-matches,
    or crosslinks from Scout.

    Parameters
    ----------
    data : pd.DataFrame
        The input data originating from Scout.

    Returns
    -------
    str
        "scout_csms_unfiltered" if a Scout unfiltered CSMs file was read, "scout_csms_filtered" if a Scout filtered CSMs file was read,
        "scout_xl" if a Scout crosslink/residue pair result file was read.


    Raises
    ------
    ValueError
        If the data source could not be determined.

    Examples
    --------
    >>> from pyXLMS.parser import detect_scout_filetype
    >>> import pandas as pd
    >>> df1 = pd.read_csv("data/scout/Cas9_Unfiltered_CSMs.csv")
    >>> detect_scout_filetype(df1)
    'scout_csms_unfiltered'

    >>> from pyXLMS.parser import detect_scout_filetype
    >>> import pandas as pd
    >>> df2 = pd.read_csv("data/scout/Cas9_Filtered_CSMs.csv")
    >>> detect_scout_filetype(df2)
    'scout_csms_filtered'

    >>> from pyXLMS.parser import detect_scout_filetype
    >>> import pandas as pd
    >>> df3 = pd.read_csv("data/scout/Cas9_Residue_Pairs.csv")
    >>> detect_scout_filetype(df3)
    'scout_xl'
    """
    ## check input
    _ok = check_input(data, "data", pd.DataFrame)

    raise ValueError(
        "Could not infer data source, are you sure you read a Scout result file?"
    )

    return "err"


def __read_scout_csms_unfiltered(
    data: pd.DataFrame,
    modifications: Dict[str, Tuple[str, float]] = SCOUT_MODIFICATION_MAPPING
) -> List[Dict[str, Any]]:
    r"""Reads crosslink-spectrum-matches from a Scout unfiltered CSMs result.

    Parameters
    ----------
    data : pandas.DataFrame
        The Scout unfiltered CSMs result data.
    modifications : dict of str, tuple, default = ``constants.SCOUT_MODIFICATION_MAPPING``
        Mapping of Scout sequence elements (e.g. ``"+15.994900"``) and modifications (e.g ``"Oxidation of Methionine"``)
        to their modifications (e.g. ``("Oxidation", 15.994915)``).

    Returns
    -------
    list of dict
        The read crosslink-spectrum-matches.

    Notes
    -----
    This function should not be called directly, it is called from ``read_scout()``.
    """
    return


def __read_scout_csms_filtered(
    data: pd.DataFrame,
    modifications: Dict[str, Tuple[str, float]] = SCOUT_MODIFICATION_MAPPING
) -> List[Dict[str, Any]]:
    r"""Reads crosslink-spectrum-matches from a Scout filtered CSMs result.

    Parameters
    ----------
    data : pandas.DataFrame
        The Scout filtered CSMs result data.
    modifications : dict of str, tuple, default = ``constants.SCOUT_MODIFICATION_MAPPING``
        Mapping of Scout sequence elements (e.g. ``"+15.994900"``) and modifications (e.g ``"Oxidation of Methionine"``)
        to their modifications (e.g. ``("Oxidation", 15.994915)``).

    Returns
    -------
    list of dict
        The read crosslink-spectrum-matches.

    Notes
    -----
    This function should not be called directly, it is called from ``read_scout()``.
    """
    return


def __read_scout_crosslinks(
    data: pd.DataFrame
) -> List[Dict[str, Any]]:
    r"""Reads crosslinks from a Scout crosslink/residue pair result.

    Parameters
    ----------
    data : pandas.DataFrame
        The Scout crosslink/residue pair result data.

    Returns
    -------
    list of dict
        The read crosslinks.

    Notes
    -----
    This function should not be called directly, it is called from ``read_scout()``.
    """
    return


def read_scout(
    files: str | List[str] | BinaryIO,
    modifications: Dict[str, Tuple[str, float]] = SCOUT_MODIFICATION_MAPPING
) -> Dict[str, Any]:
    r"""Read a Scout result file.

    Reads a Scout filtered or unfiltered crosslink-spectrum-matches result file or crosslink/residue pair result file in ``.csv`` format
    and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the Scout result file(s) or a file-like object/stream.
    modifications : dict of str, tuple, default = ``constants.SCOUT_MODIFICATION_MAPPING``
        Mapping of Scout sequence elements (e.g. ``"+15.994900"``) and modifications (e.g ``"Oxidation of Methionine"``)
        to their modifications (e.g. ``("Oxidation", 15.994915)``).

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslinks or crosslink-spectrum-matches.
    KeyError
        If one of the found post-translational-modifications could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_scout
    >>> csms_unfiltered = read_scout("data/scout/Cas9_Unfiltered_CSMs.csv")

    >>> from pyXLMS.parser import read_scout
    >>> csms_filtered = read_scout("data/scout/Cas9_Filtered_CSMs.csv")

    >>> from pyXLMS.parser import read_scout
    >>> crosslinks = read_scout("data/scout/Cas9_Residue_Pairs.csv")
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
        data = pd.read_csv(input, low_memory=False)
        ## detect input file type
        scout_file_type = detect_scout_filetype(data)
        ## process data
        if scout_file_type == "scout_csms_unfiltered":
            csms += __read_scout_csms_unfiltered(
                data, modifications
            )
        elif xi_file_type == "scout_csms_filtered":
            csms += __read_scout_csms_filtered(
                data, modifications
            )
        else:
            crosslinks += __read_scout_crosslinks(
                data
            )

    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="Scout",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )
