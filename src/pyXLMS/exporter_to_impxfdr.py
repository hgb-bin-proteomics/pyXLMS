#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .data import create_crosslink_from_csm
from .transform_util import get_available_keys
from .transform_filter import filter_target_decoy
from .exporter_to_msannika import to_msannika

from typing import Optional
from typing import Dict
from typing import Any
from typing import List


def to_impxfdr(
    data: List[Dict[str, Any]],
    filename: Optional[str],
    targets_only: bool = True,
) -> pd.DataFrame:
    r"""Exports a list of crosslinks or crosslink-spectrum-matches to IMP-X-FDR format.

    Exports a list of crosslinks or crosslink-spectrum-matches to IMP-X-FDR format for benchmarking purposes.
    The tool IMP-X-FDR is available from
    `here <https://github.com/vbc-proteomics-org/imp-x-fdr>`_.
    We recommend using version 1.1.0 and selecting "MS Annika" as input file format for the here exported file.
    A slightly modified version is available
    `here <https://github.com/hgb-bin-proteomics/MSAnnika_NC_Results/blob/master/Peplib_Beveridge/MS_Annika/Tools/IMP-X-FDR.v1.1.0.zip>_`.
    This version contains a few bug fixes and was used for the MS Annika 2.0 and MS Annika 3.0 publications.
    Requires that "alpha_proteins", "beta_proteins", "alpha_proteins_crosslink_positions" and "beta_proteins_crosslink_positions" fields
    are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslinks or crosslink-spectrum-matches.
    filename : str, or None, default = None
        If not None, the exported data will be written to a file with the specified filename.
        The filename should end in ".xlsx" as the file is exported to Microsoft Excel file format.
    targets_only : bool, default = True
        Whether or not only target crosslinks or crosslink-spectrum-matches should be exported. For
        benchmarking purposes this is usually the case. If the crosslinks or crosslink-spectrum-matches
        do not contain target-decoy labels this should be set to False.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing crosslinks or crosslink-spectrum-matches in IMP-X-FDR format.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If data contains elements of mixed data type.
    ValueError
        If the provided data contains no elements or if none of the data has target-decoy labels
        and parameter 'targets_only' is set to True.
    RuntimeError
        If not all of the required information is present in the input data.

    Examples
    --------
    >>> from pyXLMS.exporter import to_impxfdr
    >>> from pyXLMS.parser import read
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    _ok = check_input(targets_only, "targets_only", bool)
    if targets_only:
        data = filter_target_decoy(data)["Target-Target"]
    if "data_type" not in data[0] or data[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
        )
    available_keys = get_available_keys(data)
    if (
        not available_keys["alpha_proteins"]
        or not available_keys["alpha_proteins_crosslink_positions"]
        or not available_keys["beta_proteins"]
        or not available_keys["beta_proteins_crosslink_positions"]
    ):
        raise RuntimeError(
            "Can't export to IMP-X-FDR because not all necessary information is available!"
        )
    if data[0]["data_type"] == "crosslink":
        return to_msannika(data, filename, format="xlsx")
    return to_msannika(
        [create_crosslink_from_csm(csm) for csm in data], filename, format="xlsx"
    )
