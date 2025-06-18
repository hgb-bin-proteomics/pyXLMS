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
        labeled_position_a.append(xl["alpha_peptide_crosslink_position"] - 1)
        peptide_b.append(xl["beta_peptide"])
        protein_b.append(xl["beta_proteins"][0])
        labeled_position_b.append(xl["beta_peptide_crosslink_position"] - 1)
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
        xlinkdb_df.to_csv(
            __get_filename(filename, "tsv"), sep="\t", header=False, index=False
        )
    return xlinkdb_df


# wip
def to_xlinkdb(
    crosslinks: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports a list of crosslinks to XlinkDB format.

    Exports a list of crosslinks to XlinkDB format. The tool XlinkDB is accessible
    via the link
    `here <https://xlinkdb.gs.washington.edu/xlinkdb/index.php>`_.
    Requires that "alpha_proteins" and "beta_proteins" fields are set for all crosslinks.

    Parameters
    ----------
    crosslinks : list of dict of str, any
        A list of crosslinks.
    filename : str, or None
        If not None, the exported data will be written to a file with the specified filename.
        The filename should not contain a file extension and consist only of alpha-numeric
        characters (a-Z, 0-9).

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing crosslinks in XlinkDB format.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If 'crosslinks' parameter contains elements of mixed data type.
    ValueError
        If the filename contains any non-alpha-numeric characters.
    ValueError
        If the provided 'crosslinks' parameter contains no elements.
    RuntimeError
        If not all of the required information is present in the input data.

    Notes
    -----
    XlinkDB input format requires a column with probabilities that the crosslinks are correct. Since that is not available
    from most crosslink search engines, this is simply set to a constant ``1``.

    Examples
    --------
    """
    _ok = check_input(crosslinks, "crosslinks", list, dict)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    if filename is not None and not filename.isalnum():
        raise ValueError(
            "Parameter filename must only contain alpha-numeric characters and no file extension!"
        )
    if len(crosslinks) == 0:
        raise ValueError("Provided crosslinks contain no elements!")
    if "data_type" not in crosslinks[0] or crosslinks[0]["data_type"] != "crosslink":
        raise TypeError(
            "Unsupported data type for input crosslinks! Parameter crosslinks has to be a list of crosslinks!"
        )
    available_keys = get_available_keys(crosslinks)
    if not available_keys["alpha_proteins"] or not available_keys["beta_proteins"]:
        raise RuntimeError(
            "Can't export to XlinkDB because not all necessary information is available!"
        )
    return __xls_to_xlinkdb(crosslinks, filename)
