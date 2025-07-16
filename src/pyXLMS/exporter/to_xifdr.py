#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from ..data import check_input
from ..transform.filter import filter_target_decoy
from ..transform.util import get_available_keys
from .util import __get_filename

from typing import Optional
from typing import Dict
from typing import Any
from typing import List


def __csms_to_xifdr(
    csms: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports crosslink-spectrum-matches to xiFDR format.

    Parameters
    ----------
    csms : list of dict of str, any
        A list of crosslink-spectrum-matches.
    filename : str, or None
        If not None, the data will be written to a file with the specified filename.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame in xiFDR format.

    Notes
    -----
    This function should not be called directly, it is called from ``to_xifdr()``.
    """
    run = list()
    scan = list()
    peptide1 = list()
    peptide2 = list()
    peptidelink1 = list()
    peptidelink2 = list()
    isdecoy1 = list()
    isdecoy2 = list()
    precursorcharge = list()
    accession1 = list()
    accession2 = list()
    peptideposition1 = list()
    peptideposition2 = list()
    score = list()
    for csm in csms:
        run.append(csm["spectrum_file"])
        scan.append(csm["scan_nr"])
        peptide1.append(csm["alpha_peptide"])
        peptide2.append(csm["beta_peptide"])
        peptidelink1.append(csm["alpha_peptide_crosslink_position"])
        peptidelink2.append(csm["beta_peptide_crosslink_position"])
        isdecoy1.append("true" if csm["alpha_decoy"] else "false")
        isdecoy2.append("true" if csm["beta_decoy"] else "false")
        precursorcharge.append(csm["charge"])
        accession1.append(";".join(csm["alpha_proteins"]))
        accession2.append(";".join(csm["beta_proteins"]))
        peptideposition1.append(
            ";".join([str(pos) for pos in csm["alpha_proteins_peptide_positions"]])
        )
        peptideposition2.append(
            ";".join([str(pos) for pos in csm["beta_proteins_peptide_positions"]])
        )
        score.append(csm["score"])
    xifdr_df = pd.DataFrame(
        {
            "run": run,
            "scan": scan,
            "peptide1": peptide1,
            "peptide2": peptide2,
            "peptide link 1": peptidelink1,
            "peptide link 2": peptidelink2,
            "is decoy 1": isdecoy1,
            "is decoy 2": isdecoy2,
            "precursor charge": precursorcharge,
            "accession1": accession1,
            "accession2": accession2,
            "peptide position 1": peptideposition1,
            "peptide position 2": peptideposition2,
            "score": score,
        }
    )
    if filename is not None:
        xifdr_df.to_csv(__get_filename(filename, "csv"), index=False)
    return xifdr_df


def to_xifdr(
    csms: List[Dict[str, Any]],
    filename: Optional[str],
) -> pd.DataFrame:
    r"""Exports a list of crosslink-spectrum-matches to xiFDR format.

    Exports a list of crosslinks to xiFDR format. The tool xiFDR is accessible
    via the link
    `rappsilberlab.org/software/xifdr <https://www.rappsilberlab.org/software/xifdr/>`_.
    Requires that ``alpha_proteins``, ``beta_proteins``, ``alpha_proteins_peptide_positions``,
    ``beta_proteins_peptide_positions``, ``alpha_decoy``, ``beta_decoy``, ``charge`` and ``score``
    fields are set for all crosslink-spectrum-matches.

    Parameters
    ----------
    csms : list of dict of str, any
        A list of crosslink-spectrum-matches.
    filename : str, or None
        If not None, the exported data will be written to a file with the specified filename.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing crosslink-spectrum-matches in xiFDR format.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If 'csms' parameter contains elements of mixed data type.
    ValueError
        If the provided 'csms' parameter contains no elements.
    RuntimeError
        If not all of the required information is present in the input data.

    Examples
    --------
    >>> from pyXLMS.exporter import to_xifdr
    >>> from pyXLMS.parser import read
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> to_xifdr(csms, filename="msannika_xiFDR.csv")
                                           run   scan          peptide1  ... peptide position 1  peptide position 2   score
    0    XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw   2257            GQKNSR  ...                777                 777  119.83
    1    XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw   2448            GQKNSR  ...                777                 693   13.91
    2    XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw   2561             SDKNR  ...                864                 864  114.43
    3    XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw   2719            DKQSGK  ...                676                 676  200.98
    4    XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw   2792            DKQSGK  ...                676                  45   94.47
    ..                                     ...    ...               ...  ...                ...                 ...     ...
    821  XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw  23297     MDGTEELLVKLNR  ...                387                 387  286.05
    822  XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw  23454  KIECFDSVEISGVEDR  ...                575                 682  376.15
    823  XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw  23581    SSFEKNPIDFLEAK  ...               1176                1176  412.44
    824  XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw  23683    SSFEKNPIDFLEAK  ...               1176                1176  437.10
    825  XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw  27087    MEDESKLHKFKDFK  ...                 99                1176   15.89
    [826 rows x 14 columns]

    >>> from pyXLMS.exporter import to_xifdr
    >>> from pyXLMS.parser import read
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> df = to_xifdr(csms, filename=None)
    """
    _ok = check_input(csms, "csms", list, dict)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    if len(csms) == 0:
        raise ValueError("Provided crosslink-spectrum-matches contain no elements!")
    if "data_type" not in csms[0] or csms[0]["data_type"] != "crosslink-spectrum-match":
        raise TypeError(
            "Unsupported data type for input csms! Parameter csms has to be a list of crosslink-spectrum-matches!"
        )
    available_keys = get_available_keys(csms)
    if (
        not available_keys["alpha_proteins"]
        or not available_keys["beta_proteins"]
        or not available_keys["alpha_proteins_crosslink_positions"]
        or not available_keys["beta_proteins_crosslink_positions"]
        or not available_keys["alpha_decoy"]
        or not available_keys["beta_decoy"]
        or not available_keys["charge"]
        or not available_keys["score"]
        or len(filter_target_decoy(csms)["Target-Decoy"]) == 0
    ):
        raise RuntimeError(
            "Can't export to xiFDR because not all necessary information is available!"
        )
    return __csms_to_xifdr(csms, filename)
