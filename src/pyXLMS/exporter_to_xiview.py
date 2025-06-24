#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .transform_util import get_available_keys
from .exporter_to_xinet import to_xinet
from .exporter_util import __get_filename

from typing import Optional
from typing import Dict
from typing import Any
from typing import List


def __xls_to_xiview_minimal(
    xls: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports crosslinks to xiVIEW format.

    Parameters
    ----------
    xls : list of dict of str, any
        A list of crosslinks.
    filename : str, or None
        If not None, the data will be written to a file with the specified filename.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame in xiVIEW format.

    Notes
    -----
    This function should not be called directly, it is called from ``to_xiview()``.
    """
    abspos1 = list()
    abspos2 = list()
    protein1 = list()
    protein2 = list()
    decoy1 = list()
    decoy2 = list()
    score = list()
    has_decoys = True
    has_scores = True

    for xl in xls:
        abspos1.append(
            ";".join([str(pos) for pos in xl["alpha_proteins_crosslink_positions"]])
        )
        abspos2.append(
            ";".join([str(pos) for pos in xl["beta_proteins_crosslink_positions"]])
        )
        protein1.append(";".join(xl["alpha_proteins"]))
        protein2.append(";".join(xl["beta_proteins"]))
        if xl["alpha_decoy"] is not None and xl["beta_decoy"] is not None:
            if xl["alpha_decoy"]:
                decoy1.append("TRUE")
            else:
                decoy1.append("FALSE")
            if xl["beta_decoy"]:
                decoy2.append("TRUE")
            else:
                decoy2.append("FALSE")
        else:
            has_decoys = False
        if xl["score"] is not None:
            score.append(xl["score"])
        else:
            has_scores = False

    xiview_df = pd.DataFrame(
        {
            "AbsPos1": abspos1,
            "AbsPos2": abspos2,
            "Protein1": protein1,
            "Protein2": protein2,
        }
    )

    if has_decoys:
        xiview_df["Decoy1"] = decoy1
        xiview_df["Decoy2"] = decoy2
    if has_scores:
        xiview_df["Score"] = score

    if filename is not None:
        xiview_df.to_csv(__get_filename(filename, "csv"), index=False)
    return xiview_df


def to_xiview(
    crosslinks: List[Dict[str, Any]],
    filename: Optional[str],
    minimal: bool = True,
) -> pd.DataFrame:
    r"""Exports a list of crosslinks to xiVIEW format.

    Exports a list of crosslinks to xiVIEW format. The tool xiVIEW is accessible
    via the link
    `xiview.org/ <https://xiview.org/>`_.
    Requires that ``alpha_proteins``, ``beta_proteins``, ``alpha_proteins_crosslink_positions`` and
    ``beta_proteins_crosslink_positions`` fields are set for all crosslinks.

    Parameters
    ----------
    crosslinks : list of dict of str, any
        A list of crosslinks.
    filename : str, or None
        If not None, the exported data will be written to a file with the specified filename.
    minimal : bool, default = True
        Which xiVIEW format to return, if ``minimal = True`` the minimal xiVIEW format is returned. Otherwise
        the "CSV without peak lists" format is returned (internally this just calls ``exporter.to_xinet()``).
        For more information on the xiVIEW formats please refer to the `xiVIEW specification <https://xiview.org/csv-formats.php>`_.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing crosslinks in xiVIEW format.

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
    The optional ``Score`` column in the xiVIEW table will only be available if all crosslinks have assigned scores,
    the optional ``Decoy*`` columns will only be available if all crosslinks have assigned target and decoy labels.

    Examples
    --------
    >>> from pyXLMS.exporter import to_xiview
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> from pyXLMS.transform import filter_proteins
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> crosslinks = targets_only(pr)["crosslinks"]
    >>> cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    >>> to_xiview(cas9, filename="crosslinks_xiVIEW.csv")
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
            "Can't export to xiVIEW because not all necessary information is available!"
        )
    if minimal:
        return __xls_to_xiview_minimal(crosslinks, filename)
    return to_xinet(crosslinks, filename)
