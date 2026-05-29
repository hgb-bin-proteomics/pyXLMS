#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from ..data._crosslink import Crosslink
from ..data._util import check_input
from ._util import __get_filename
from ._to_msannika import get_msannika_crosslink_sequence

from typing import Optional
from typing import List


def __xls_to_xmas(
    xls: List[Crosslink],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports crosslinks to XMAS format.

    Parameters
    ----------
    xls : list of Crosslink
        A list of crosslinks.
    filename : str, or None
        If not None, the data will be written to a file with the specified filename.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame in XMAS format.

    Notes
    -----
    This function should not be called directly, it is called from ``to_xmas()``.
    """
    sequence_a = list()
    sequence_b = list()
    for xl in xls:
        sequence_a.append(
            get_msannika_crosslink_sequence(
                xl["alpha_peptide"], xl["alpha_peptide_crosslink_position"]
            )
        )
        sequence_b.append(
            get_msannika_crosslink_sequence(
                xl["beta_peptide"], xl["beta_peptide_crosslink_position"]
            )
        )
    xmas_df = pd.DataFrame(
        {
            "Sequence A": sequence_a,
            "Sequence B": sequence_b,
        }
    )
    if filename is not None:
        xmas_df.to_excel(
            __get_filename(filename, "xlsx"), engine="openpyxl", index=False
        )
    return xmas_df


def to_xmas(
    crosslinks: List[Crosslink],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports a list of crosslinks to XMAS format.

    Exports a list of crosslinks to XMAS format for visualization in ChimeraX. The tool XMAS
    is available from
    `github.com/ScheltemaLab/ChimeraX_XMAS_bundle <https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle>`_.

    Parameters
    ----------
    crosslinks : list of Crosslink
        A list of crosslinks.
    filename : str, or None
        If not None, the exported data will be written to a file with the specified filename.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing crosslinks in XMAS format.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If 'crosslinks' parameter contains elements of mixed data type.
    ValueError
        If the provided 'crosslinks' parameter contains no elements.

    Examples
    --------
    >>> from pyXLMS.exporter import to_xmas
    >>> from pyXLMS.data import create_crosslink_min
    >>> xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    >>> xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    >>> crosslinks = [xl1, xl2]
    >>> to_xmas(crosslinks, filename="crosslinks_xmas.xlsx")
       Sequence A  Sequence B
    0  [K]PEPTIDE  P[K]EPTIDE
    1  PE[K]PTIDE  PEP[K]TIDE

    >>> from pyXLMS.exporter import to_xmas
    >>> from pyXLMS.data import create_crosslink_min
    >>> xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    >>> xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    >>> crosslinks = [xl1, xl2]
    >>> to_xmas(crosslinks, filename=None)
       Sequence A  Sequence B
    0  [K]PEPTIDE  P[K]EPTIDE
    1  PE[K]PTIDE  PEP[K]TIDE
    """
    _ok = check_input(crosslinks, "crosslinks", list, Crosslink)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    if len(crosslinks) == 0:
        raise ValueError("Provided crosslinks contain no elements!")
    return __xls_to_xmas(crosslinks, filename)
