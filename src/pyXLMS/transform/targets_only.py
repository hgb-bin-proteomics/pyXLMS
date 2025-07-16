#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data import check_input
from ..data import check_input_multi
from ..data import create_parser_result
from .filter import filter_target_decoy

from typing import Dict
from typing import List
from typing import Any


def targets_only(
    data: List[Dict[str, Any]] | Dict[str, Any],
) -> List[Dict[str, Any]] | Dict[str, Any]:
    r"""Get target crosslinks or crosslink-spectrum-matches.

    Get target crosslinks or crosslink-spectrum-matches from a list of target and decoy crosslinks or
    crosslink-spectrum-matches, or a parser_result. This effectively filters out any target-decoy and
    decoy-decoy matches and is essentially a convenience wrapper for
    ``transform.filter_target_decoy()["Target-Target"]``.

    Parameters
    ----------
    data : dict of str, any, or list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks, or a parser_result.

    Returns
    -------
    list of dict of str, any, or dict of str, any
        If a list of crosslink-spectrum-matches or crosslinks was provided, a list of target
        crosslink-spectrum-matches or crosslinks is returned. If a parser_result was provided,
        a parser_result with target crosslink-spectrum-matches and/or target crosslinks will
        be returned.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    RuntimeError
        If no target crosslinks or crosslink-spectrum-matches were found.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> result = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> targets = targets_only(result["crosslink-spectrum-matches"])
    >>> len(targets)
    786

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> result = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> targets = targets_only(result["crosslinks"])
    >>> len(targets)
    265

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> result = read(["data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx"], engine="MS Annika", crosslinker="DSS")
    >>> result_targets = targets_only(result)
    >>> len(result_targets["crosslink-spectrum-matches"])
    786
    >>> len(result_targets["crosslinks"])
    265
    """
    _ok = check_input_multi(data, "data", [dict, list])
    if isinstance(data, list):
        _ok = check_input(data, "data", list, dict)
        if len(data) == 0:
            return data
        if "data_type" not in data[0] or data[0]["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match, "
                "or a parser_result!"
            )
        targets = filter_target_decoy(data)["Target-Target"]
        if len(targets) == 0:
            raise RuntimeError(
                "No target matches found! Are you sure your data is labelled?"
            )
        return targets
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError(
            "Can't filter target matches for dict. Dict has to be a valid 'parser_result'!"
        )
    new_csms = (
        filter_target_decoy(data["crosslink-spectrum-matches"])["Target-Target"]
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        filter_target_decoy(data["crosslinks"])["Target-Target"]
        if data["crosslinks"] is not None
        else None
    )
    if new_csms is not None:
        if not isinstance(new_csms, list):
            raise RuntimeError(
                "Something went wrong while getting target crosslink-spectrum-matches.\n"
                f"Expected data type: list. Got: {type(new_csms)}."
            )
        if len(new_csms) == 0:
            raise RuntimeError(
                "No target crosslink-spectrum-matches found! Are you sure they are labelled?"
            )
    if new_xls is not None:
        if not isinstance(new_xls, list):
            raise RuntimeError(
                "Something went wrong while getting target crosslinks.\n"
                f"Expected data type: list. Got: {type(new_xls)}."
            )
        if len(new_xls) == 0:
            raise RuntimeError(
                "No target crosslinks found! Are you sure they are labelled?"
            )
    return create_parser_result(
        search_engine=data["search_engine"], csms=new_csms, crosslinks=new_xls
    )
