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
        AbsPos1 AbsPos2 Protein1 Protein2 Decoy1 Decoy2   Score
    0       779     779     Cas9     Cas9  FALSE  FALSE  119.83
    1       866     866     Cas9     Cas9  FALSE  FALSE  114.43
    2       677     677     Cas9     Cas9  FALSE  FALSE  200.98
    3       677      48     Cas9     Cas9  FALSE  FALSE   94.47
    4        34      34     Cas9     Cas9  FALSE  FALSE  110.48
    ..      ...     ...      ...      ...    ...    ...     ...
    248     396     396     Cas9     Cas9  FALSE  FALSE  305.63
    249     688     952     Cas9     Cas9  FALSE  FALSE  110.46
    250     793    1180     Cas9     Cas9  FALSE  FALSE  288.36
    251     575     688     Cas9     Cas9  FALSE  FALSE  376.15
    252    1180    1180     Cas9     Cas9  FALSE  FALSE  437.10
    [253 rows x 7 columns]

    >>> from pyXLMS.exporter import to_xiview
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> from pyXLMS.transform import filter_proteins
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> crosslinks = targets_only(pr)["crosslinks"]
    >>> cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    >>> df = to_xiview(cas9, filename=None)

    >>> from pyXLMS.exporter import to_xiview
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> from pyXLMS.transform import filter_proteins
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> crosslinks = targets_only(pr)["crosslinks"]
    >>> cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    >>> to_xiview(cas9, filename="crosslinks_xiVIEW.csv", minimal=False)
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
