#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._parser_result import ParserResult
from ..data._util import check_input_multi
from ..data._parser_result import create_parser_result
from ._filter import filter_target_decoy
from ._util import assert_csms
from ._util import assert_xls
from ._util import assert_csms_or_xls

from typing import List


def targets_only(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult,
) -> List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult:
    r"""Get target crosslinks or crosslink-spectrum-matches.

    Get target crosslinks or crosslink-spectrum-matches from a list of target and decoy crosslinks or
    crosslink-spectrum-matches, or a parser_result. This effectively filters out any target-decoy and
    decoy-decoy matches and is essentially a convenience wrapper for
    ``transform.filter_target_decoy()["Target-Target"]``.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        A list of crosslink-spectrum-matches or crosslinks, or a parser_result.

    Returns
    -------
    list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
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
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> targets = targets_only(result["crosslink-spectrum-matches"])
    >>> len(targets)
    786

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> targets = targets_only(result["crosslinks"])
    >>> len(targets)
    265

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import targets_only
    >>> result = read(
    ...     [
    ...         "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...         "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
    ...     ],
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> result_targets = targets_only(result)
    >>> len(result_targets["crosslink-spectrum-matches"])
    786
    >>> len(result_targets["crosslinks"])
    265
    """
    _ok = check_input_multi(data, "data", [ParserResult, list])
    if isinstance(data, list):
        csms_or_xls = assert_csms_or_xls(data)
        if len(csms_or_xls) == 0:
            return csms_or_xls
        targets = filter_target_decoy(csms_or_xls)["Target-Target"]
        if len(targets) == 0:
            raise RuntimeError(
                "No target matches found! Are you sure your data is labelled?"
            )
        return targets
    new_csms = (
        assert_csms(
            filter_target_decoy(data["crosslink-spectrum-matches"])["Target-Target"]
        )
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        assert_xls(filter_target_decoy(data["crosslinks"])["Target-Target"])
        if data["crosslinks"] is not None
        else None
    )
    if new_csms is not None:
        if len(new_csms) == 0:
            raise RuntimeError(
                "No target crosslink-spectrum-matches found! Are you sure they are labelled?"
            )
    if new_xls is not None:
        if len(new_xls) == 0:
            raise RuntimeError(
                "No target crosslinks found! Are you sure they are labelled?"
            )
    return create_parser_result(
        search_engine=data["search_engine"],
        csms=new_csms,
        crosslinks=new_xls,
    )
