#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .transform_util import get_available_keys
from .exporter_util import __get_filename

from typing import Optional
from typing import Dict
from typing import Any
from typing import List


def __xls_to_xinet(
    xls: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports crosslinks to xiNET format.

    Parameters
    ----------
    xls : list of dict of str, any
        A list of crosslinks.
    filename : str, or None
        If not None, the data will be written to a file with the specified filename.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame in xiNET format.

    Notes
    -----
    This function should not be called directly, it is called from ``to_xinet()``.
    """
    protein1 = list()
    pepseq1 = list()
    linkpos1 = list()
    protein2 = list()
    pepseq2 = list()
    linkpos2 = list()
    score = list()
    id = list()
    has_scores = True
    for i, xl in enumerate(xls):
        protein1.append(";".join(xl["alpha_proteins"]))
        pepseq1.append(xl["alpha_peptide"])
        linkpos1.append(
            ";".join([str(pos) for pos in xl["alpha_proteins_crosslink_positions"]])
        )
        protein2.append(";".join(xl["beta_proteins"]))
        pepseq2.append(xl["beta_peptide"])
        linkpos2.append(
            ";".join([str(pos) for pos in xl["beta_proteins_crosslink_positions"]])
        )
        if xl["score"] is not None:
            score.append(xl["score"])
        else:
            has_scores = False
        id.append(i)
    xinet_df = pd.DataFrame()
    if has_scores:
        xinet_df = pd.DataFrame(
            {
                "Protein1": protein1,
                "PepSeq1": pepseq1,
                "LinkPos1": linkpos1,
                "Protein2": protein2,
                "PepSeq2": pepseq2,
                "LinkPos2": linkpos2,
                "Score": score,
                "Id": id,
            }
        )
    else:
        xinet_df = pd.DataFrame(
            {
                "Protein1": protein1,
                "PepSeq1": pepseq1,
                "LinkPos1": linkpos1,
                "Protein2": protein2,
                "PepSeq2": pepseq2,
                "LinkPos2": linkpos2,
                "Id": id,
            }
        )
    if filename is not None:
        xinet_df.to_csv(__get_filename(filename, "csv"), index=False)
    return xinet_df


def to_xinet(
    crosslinks: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports a list of crosslinks to xiNET format.

    Exports a list of crosslinks to xiNET format. The tool xiNET is accessible
    via the link
    `crosslinkviewer.org <https://crosslinkviewer.org/>`_.
    Requires that ``alpha_proteins``, ``beta_proteins``, ``alpha_proteins_crosslink_positions`` and
    ``beta_proteins_crosslink_positions`` fields are set for all crosslinks.

    Parameters
    ----------
    crosslinks : list of dict of str, any
        A list of crosslinks.
    filename : str, or None
        If not None, the exported data will be written to a file with the specified filename.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing crosslinks in xiNET format.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If 'crosslinks' parameter contains elements of mixed data type.
    ValueError
        If the provided 'crosslinks' parameter contains no elements.
    RuntimeError
        If not all of the required information is present in the input data.

    Notes
    -----
    The optional ``Score`` column in the xiNET table will only be available if all crosslinks have assigned scores.

    Examples
    --------
    >>> from pyXLMS.exporter import to_xinet
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> crosslinks = targets_only(pr)["crosslinks"]
    >>> to_xinet(crosslinks, filename="crosslinks_xiNET.csv")
        Protein1           PepSeq1 LinkPos1 Protein2         PepSeq2 LinkPos2   Score   Id
    0       Cas9            GQKNSR      779     Cas9          GQKNSR      779  119.83    0
    1       Cas9             SDKNR      866     Cas9           SDKNR      866  114.43    1
    2       Cas9            DKQSGK      677     Cas9          DKQSGK      677  200.98    2
    3       Cas9            DKQSGK      677     Cas9           HSIKK       48   94.47    3
    4       Cas9             VPSKK       34     Cas9           VPSKK       34  110.48    4
    ..       ...               ...      ...      ...             ...      ...     ...  ...
    260     Cas9     MDGTEELLVKLNR      396     Cas9   MDGTEELLVKLNR      396  305.63  260
    261     Cas9    TILDFLKSDGFANR      688     Cas9       YDENDKLIR      952  110.46  261
    262     Cas9    IEEGIKELGSQILK      793     Cas9  SSFEKNPIDFLEAK     1180  288.36  262
    263     Cas9  KIECFDSVEISGVEDR      575     Cas9  TILDFLKSDGFANR      688  376.15  263
    264     Cas9    SSFEKNPIDFLEAK     1180     Cas9  SSFEKNPIDFLEAK     1180  437.10  264
    [265 rows x 8 columns]

    >>> from pyXLMS.exporter import to_xinet
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> crosslinks = targets_only(pr)["crosslinks"]
    >>> df = to_xinet(crosslinks, filename=None)
    """
    _ok = check_input(crosslinks, "crosslinks", list, dict)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    if len(crosslinks) == 0:
        raise ValueError("Provided crosslinks contain no elements!")
    if "data_type" not in crosslinks[0] or crosslinks[0]["data_type"] != "crosslink":
        raise TypeError(
            "Unsupported data type for input crosslinks! Parameter crosslinks has to be a list of crosslinks!"
        )
    available_keys = get_available_keys(crosslinks)
    if (
        not available_keys["alpha_proteins"]
        or not available_keys["beta_proteins"]
        or not available_keys["alpha_proteins_crosslink_positions"]
        or not available_keys["beta_proteins_crosslink_positions"]
    ):
        raise RuntimeError(
            "Can't export to xiNET because not all necessary information is available!"
        )
    return __xls_to_xinet(crosslinks, filename)
