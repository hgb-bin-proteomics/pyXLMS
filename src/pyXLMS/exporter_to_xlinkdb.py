#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .transform_util import get_available_keys

from typing import Optional
from typing import Dict
from typing import Any
from typing import List


def __xls_to_xlinkdb(
    xls: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports crosslinks to XlinkDB format.

    Parameters
    ----------
    xls : list of dict of str, any
        A list of crosslinks.
    filename : str, or None
        If not None, the data will be written to a file with the specified filename.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame in XlinkDB format.

    Notes
    -----
    This function should not be called directly, it is called from ``to_xlinkdb()``.
    """
    peptide_a = list()
    protein_a = list()
    labeled_position_a = list()
    peptide_b = list()
    protein_b = list()
    labeled_position_b = list()
    probability = list()
    for xl in xls:
        peptide_a.append(xl["alpha_peptide"])
        protein_a.append(xl["alpha_proteins"][0])
        labeled_position_a.append(xl["alpha_peptide_crosslink_position"]-1)
        peptide_b.append(xl["beta_peptide"])
        protein_b.append(xl["beta_proteins"][0])
        labeled_position_a.append(xl["beta_peptide_crosslink_position"]-1)
        probability.append(1)
    xlinkdb_df = pd.DataFrame(
        {
            "Peptide A": peptide_a,
            "Protein A": protein_a,
            "Labeled Position A": labeled_position_a,
            "Peptide B": peptide_b,
            "Protein B": protein_b,
            "Labeled Position B": labeled_position_b,
            "Probability": probability,
        }
    )
    if filename is not None:
        xlinkdb_df.to_csv(filename, sep="\t", header=False, index=False)
    return xlinkdb_df


# wip
def to_xlinkdb(
    crosslinks: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports a list of crosslinks to XMAS format.

    Exports a list of crosslinks to XMAS format for visualization in ChimeraX. The tool XMAS
    is available from
    `here <https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle>`_.

    Parameters
    ----------
    crosslinks : list of dict of str, any
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
    _ok = check_input(crosslinks, "crosslinks", list, dict)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    if len(crosslinks) == 0:
        raise ValueError("Provided crosslinks contain no elements!")
    if "data_type" not in crosslinks[0] or crosslinks[0]["data_type"] != "crosslink":
        raise TypeError(
            "Unsupported data type for input crosslinks! Parameter crosslinks has to be a list of crosslinks!"
        )
    if not assert_data_type_same(crosslinks):
        raise TypeError("Not all elements in data have the same data type!")
    return __xls_to_xmas(crosslinks, filename)
