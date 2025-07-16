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


def to_xlinkdb(
    crosslinks: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports a list of crosslinks to XlinkDB format.

    Exports a list of crosslinks to XlinkDB format. The tool XlinkDB is accessible
    via the link
    `xlinkdb.gs.washington.edu/xlinkdb <https://xlinkdb.gs.washington.edu/xlinkdb/index.php>`_.
    Requires that ``alpha_proteins`` and ``beta_proteins`` fields are set for all crosslinks.

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
    >>> from pyXLMS.exporter import to_xlinkdb
    >>> from pyXLMS.parser import read
    >>> pr = read("data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv", engine="xiSearch/xiFDR", crosslinker="DSS")
    >>> crosslinks = pr["crosslinks"]
    >>> to_xlinkdb(crosslinks, filename="crosslinksForXlinkDB")
                   Peptide A Protein A  Labeled Position A      Peptide B Protein B  Labeled Position B  Probability
    0            VVDELVKVMGR      Cas9                   6    VVDELVKVMGR      Cas9                   6            1
    1    MLASAGELQKGNELALPSK      Cas9                   9    VVDELVKVMGR      Cas9                   6            1
    2          MDGTEELLVKLNR      Cas9                   9  MDGTEELLVKLNR      Cas9                   9            1
    3           MTNFDKNLPNEK      Cas9                   5       SKLVSDFR      Cas9                   1            1
    4               DFQFYKVR      Cas9                   5    MIAKSEQEIGK      Cas9                   3            1
    ..                   ...       ...                 ...            ...       ...                 ...          ...
    222        LPKYSLFELENGR      Cas9                   2          SDKNR      Cas9                   2            1
    223               DKQSGK      Cas9                   1         DKQSGK      Cas9                   1            1
    224               AGFIKR      Cas9                   4   SDNVPSEEVVKK      Cas9                  10            1
    225                EKIEK      Cas9                   1          KVTVK      Cas9                   0            1
    226                LSKSR      Cas9                   2          LSKSR      Cas9                   2            1
    [227 rows x 7 columns]

    >>> from pyXLMS.exporter import to_xlinkdb
    >>> from pyXLMS.parser import read
    >>> pr = read("data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv", engine="xiSearch/xiFDR", crosslinker="DSS")
    >>> crosslinks = pr["crosslinks"]
    >>> df = to_xlinkdb(crosslinks, filename=None)
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
