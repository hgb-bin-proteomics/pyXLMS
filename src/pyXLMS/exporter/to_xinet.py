#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from ..data import check_input
from ..transform.util import get_available_keys
from .util import __get_filename

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
    peppos1 = list()
    pepseq1 = list()
    linkpos1 = list()
    protein2 = list()
    peppos2 = list()
    pepseq2 = list()
    linkpos2 = list()
    score = list()
    id = list()
    has_scores = True
    for i, xl in enumerate(xls):
        pos1 = xl["alpha_peptide_crosslink_position"]
        protein1.append(";".join(xl["alpha_proteins"]))
        peppos1.append(
            ";".join(
                [
                    str(pos - pos1 + 1)
                    for pos in xl["alpha_proteins_crosslink_positions"]
                ]
            )
        )
        pepseq1.append(xl["alpha_peptide"])
        linkpos1.append(pos1)
        pos2 = xl["beta_peptide_crosslink_position"]
        protein2.append(";".join(xl["beta_proteins"]))
        peppos2.append(
            ";".join(
                [str(pos - pos2 + 1) for pos in xl["beta_proteins_crosslink_positions"]]
            )
        )
        pepseq2.append(xl["beta_peptide"])
        linkpos2.append(pos2)
        if xl["score"] is not None:
            score.append(xl["score"])
        else:
            has_scores = False
        id.append(i + 1)
    xinet_df = pd.DataFrame()
    if has_scores:
        xinet_df = pd.DataFrame(
            {
                "Protein1": protein1,
                "PepPos1": peppos1,
                "PepSeq1": pepseq1,
                "LinkPos1": linkpos1,
                "Protein2": protein2,
                "PepPos2": peppos2,
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
                "PepPos1": peppos1,
                "PepSeq1": pepseq1,
                "LinkPos1": linkpos1,
                "Protein2": protein2,
                "PepPos2": peppos2,
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
    >>> from pyXLMS.transform import filter_proteins
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> crosslinks = targets_only(pr)["crosslinks"]
    >>> cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    >>> to_xinet(cas9, filename="crosslinks_xiNET.csv")
        Protein1 PepPos1           PepSeq1  LinkPos1 Protein2 PepPos2         PepSeq2  LinkPos2   Score   Id
    0       Cas9     777            GQKNSR         3     Cas9     777          GQKNSR         3  119.83    1
    1       Cas9     864             SDKNR         3     Cas9     864           SDKNR         3  114.43    2
    2       Cas9     676            DKQSGK         2     Cas9     676          DKQSGK         2  200.98    3
    3       Cas9     676            DKQSGK         2     Cas9      45           HSIKK         4   94.47    4
    4       Cas9      31             VPSKK         4     Cas9      31           VPSKK         4  110.48    5
    ..       ...     ...               ...       ...      ...     ...             ...       ...     ...  ...
    248     Cas9     387     MDGTEELLVKLNR        10     Cas9     387   MDGTEELLVKLNR        10  305.63  249
    249     Cas9     682    TILDFLKSDGFANR         7     Cas9     947       YDENDKLIR         6  110.46  250
    250     Cas9     788    IEEGIKELGSQILK         6     Cas9    1176  SSFEKNPIDFLEAK         5  288.36  251
    251     Cas9     575  KIECFDSVEISGVEDR         1     Cas9     682  TILDFLKSDGFANR         7  376.15  252
    252     Cas9    1176    SSFEKNPIDFLEAK         5     Cas9    1176  SSFEKNPIDFLEAK         5  437.10  253
    [253 rows x 10 columns]

    >>> from pyXLMS.exporter import to_xinet
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> from pyXLMS.transform import filter_proteins
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> crosslinks = targets_only(pr)["crosslinks"]
    >>> cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    >>> df = to_xinet(cas9, filename=None)
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
