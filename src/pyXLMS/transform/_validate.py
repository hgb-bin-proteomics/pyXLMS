#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
import numpy as np
from tqdm import tqdm

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._util import check_input_multi
from ..data._parser_result import create_parser_result
from ._util import check_available_keys
from ._util import assert_csms
from ._util import assert_xls
from ._util import assert_csms_or_xls
from ._filter import filter_target_decoy
from ._filter import filter_crosslink_type

from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __verify_fdr_strict(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    fdr: float,
    cutoff: float,
    score: Literal["higher_better", "lower_better"],
) -> bool:
    r"""Verifies that a list of crosslinks or crosslink-spectrum-matches has the estimated false discovery rate using a strict approach.

    Verifies that a list of crosslinks or crosslink-spectrum-matches has the estimated false discovery rate (FDR) using the
    formula (TD+DD)/TT given a score cutoff.
    Requires that "score", "alpha_decoy" and "beta_decoy" fields are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to validate.
    fdr : float
        The target FDR, must be given as a real number between 0 and 1.
    cutoff : float
        Score cutoff that defines which crosslinks or crosslink-spectrum-matches fall within the FDR validated result.
    score : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    bool
        Returns True if the given score cutoff yields only crosslinks or crosslink-spectrum-matches
        that are within the target FDR. Returns False if the given score cutoff yields crosslinks or
        crosslink-spectrum-matches that produce a higher estimated FDR than the desired target FDR.

    Notes
    -----
    This function should not be called directly, it is called from ``__validate_strict()``.
    """
    D = 0
    T = 0
    for item in data:
        if score == "higher_better" and item["score"] >= cutoff:
            if not item["alpha_decoy"] and not item["beta_decoy"]:
                T += 1
            else:
                D += 1
        elif score == "lower_better" and item["score"] <= cutoff:
            if not item["alpha_decoy"] and not item["beta_decoy"]:
                T += 1
            else:
                D += 1
        else:
            # do nothing
            pass
    return D / T < fdr


def __validate_strict(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    fdr: float,
    score: Literal["higher_better", "lower_better"],
) -> List[CrosslinkSpectrumMatch] | List[Crosslink]:
    r"""Validate a list of crosslinks or crosslink-spectrum-matches by strict false discovery rate estimation.

    Validate a list of crosslinks or crosslink-spectrum-matches by strict false discovery rate (FDR) estimation using the
    formula (TD+DD)/TT.
    Requires that "score", "alpha_decoy" and "beta_decoy" fields are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to validate.
    fdr : float
        The target FDR, must be given as a real number between 0 and 1. The default of 0.01 corresponds to 1% FDR.
    score : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    list of dict of str, any
        A list of validated crosslink-spectrum-matches or crosslinks.

    Warns
    -----
    RuntimeWarning
        If none of the data passes the desired FDR threshold.

    Notes
    -----
    This function should not be called directly, it is called from ``validate()``.
    """
    scores = list()
    td = list()
    for item in data:
        scores.append(item["score"])
        if not item["alpha_decoy"] and not item["beta_decoy"]:
            td.append(0)
        else:
            td.append(1)
    scores = np.array(scores)
    cutoff = 0
    if score == "higher_better":
        td = np.array(td)[np.argsort(scores, kind="stable")]
        scores = scores[np.argsort(scores, kind="stable")]
        cutoff = scores[0]  # scores.max()
    else:
        td = np.array(td)[np.argsort(scores, kind="stable")[::-1]]
        scores = scores[np.argsort(scores, kind="stable")[::-1]]
        cutoff = scores[0]  # scores.min()
    nr_items = len(td)
    for i in tqdm(
        range(nr_items),
        total=nr_items,
        desc="Iterating over scores for FDR calculation...",
    ):
        if (nr_items - i - td[i:].sum()) <= 0.0:
            warnings.warn(
                RuntimeWarning(
                    "None of the data passes the desired FDR threshold! This is usually due to decoys with very good scores."
                )
            )
            return []
        if td[i:].sum() / (nr_items - i - td[i:].sum()) < fdr:
            # we need to verify in this case because there might be multiple
            # items with the same score
            if __verify_fdr_strict(data, fdr, scores[i], score):
                cutoff = scores[i]
                break
    validated_items = list()
    for item in data:
        if score == "higher_better" and item["score"] >= cutoff:
            validated_items.append(item)
        elif score == "lower_better" and item["score"] <= cutoff:
            validated_items.append(item)
        else:
            # do nothing
            pass
    return validated_items


def __verify_fdr_relaxed(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    fdr: float,
    cutoff: float,
    score: Literal["higher_better", "lower_better"],
) -> bool:
    r"""Verifies that a list of crosslinks or crosslink-spectrum-matches has the estimated false discovery rate using a relaxed approach.

    Verifies that a list of crosslinks or crosslink-spectrum-matches has the estimated false discovery rate (FDR) using the
    formula (TD-DD)/TT given a score cutoff.
    Requires that "score", "alpha_decoy" and "beta_decoy" fields are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to validate.
    fdr : float
        The target FDR, must be given as a real number between 0 and 1.
    cutoff : float
        Score cutoff that defines which crosslinks or crosslink-spectrum-matches fall within the FDR validated result.
    score : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    bool
        Returns True if the given score cutoff yields only crosslinks or crosslink-spectrum-matches
        that are within the target FDR. Returns False if the given score cutoff yields crosslinks or
        crosslink-spectrum-matches that produce a higher estimated FDR than the desired target FDR.

    Raises
    ------
    RuntimeError
        If the number of DD matches exceeds the number of TD matches.
        FDR can not be estimated with the formula '(TD-DD)/TT' in these cases.

    Notes
    -----
    This function should not be called directly, it is called from ``__validate_relaxed()``.
    """
    D = 0
    DT = 0
    T = 0
    for item in data:
        if score == "higher_better" and item["score"] >= cutoff:
            if not item["alpha_decoy"] and not item["beta_decoy"]:
                T += 1
            elif item["alpha_decoy"] and item["beta_decoy"]:
                D += 1
            else:
                DT += 1
        elif score == "lower_better" and item["score"] <= cutoff:
            if not item["alpha_decoy"] and not item["beta_decoy"]:
                T += 1
            elif item["alpha_decoy"] and item["beta_decoy"]:
                D += 1
            else:
                DT += 1
        else:
            # do nothing
            pass
    if (DT - D) < 0.0:
        raise RuntimeError(
            f"Number of DD matches is greater than the number of TD matches for score {cutoff}! "
            "Invalid FDR estimation! Please use formula 'D/T' instead!"
        )
    return (DT - D) / T < fdr


def __validate_relaxed(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    fdr: float,
    score: Literal["higher_better", "lower_better"],
) -> List[CrosslinkSpectrumMatch] | List[Crosslink]:
    r"""Validate a list of crosslinks or crosslink-spectrum-matches by relaxed false discovery rate estimation.

    Validate a list of crosslinks or crosslink-spectrum-matches by relaxed false discovery rate (FDR) estimation using the
    formula (TD-DD)/TT.
    Requires that "score", "alpha_decoy" and "beta_decoy" fields are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to validate.
    fdr : float
        The target FDR, must be given as a real number between 0 and 1.
    score : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    list of dict of str, any
        A list of validated crosslink-spectrum-matches or crosslinks.

    Raises
    ------
    RuntimeError
        If the number of DD matches exceeds the number of TD matches.
        FDR can not be estimated with the formula '(TD-DD)/TT' in these cases.

    Warns
    -----
    RuntimeWarning
        If none of the data passes the desired FDR threshold.

    Notes
    -----
    This function should not be called directly, it is called from ``validate()``.
    """
    scores = list()
    td = list()
    tdd = list()
    for item in data:
        scores.append(item["score"])
        if not item["alpha_decoy"] and not item["beta_decoy"]:
            td.append(0)
            tdd.append(0)
        elif item["alpha_decoy"] and item["beta_decoy"]:
            td.append(1)
            tdd.append(-1)
        else:
            td.append(1)
            tdd.append(1)
    scores = np.array(scores)
    cutoff = 0
    if score == "higher_better":
        td = np.array(td)[np.argsort(scores, kind="stable")]
        tdd = np.array(tdd)[np.argsort(scores, kind="stable")]
        scores = scores[np.argsort(scores, kind="stable")]
        cutoff = scores[0]  # scores.max()
    else:
        td = np.array(td)[np.argsort(scores, kind="stable")[::-1]]
        tdd = np.array(tdd)[np.argsort(scores, kind="stable")[::-1]]
        scores = scores[np.argsort(scores, kind="stable")[::-1]]
        cutoff = scores[0]  # scores.min()
    nr_items = len(td)
    for i in tqdm(
        range(nr_items),
        total=nr_items,
        desc="Iterating over scores for FDR calculation...",
    ):
        if tdd[i:].sum() < 0.0:
            raise RuntimeError(
                f"Number of DD matches is greater than the number of TD matches for score {scores[i]}! "
                "Invalid FDR estimation! Please use formula 'D/T' instead!"
            )
        if (nr_items - i - td[i:].sum()) <= 0.0:
            warnings.warn(
                RuntimeWarning(
                    "None of the data passes the desired FDR threshold! This is usually due to decoys with very good scores."
                )
            )
            return []
        if tdd[i:].sum() / (nr_items - i - td[i:].sum()) < fdr:
            # we need to verify in this case because there might be multiple
            # items with the same score
            if __verify_fdr_relaxed(data, fdr, scores[i], score):
                cutoff = scores[i]
                break
    validated_items = list()
    for item in data:
        if score == "higher_better" and item["score"] >= cutoff:
            validated_items.append(item)
        elif score == "lower_better" and item["score"] <= cutoff:
            validated_items.append(item)
        else:
            # do nothing
            pass
    return validated_items


def validate(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult,
    fdr: float = 0.01,
    formula: Literal["D/T", "(TD+DD)/TT", "(TD-DD)/TT"] = "D/T",
    score: Literal["higher_better", "lower_better"] = "higher_better",
    separate_intra_inter: bool = False,
    ignore_missing_labels: bool = False,
) -> List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult:
    r"""Validate a list of crosslinks or crosslink-spectrum-matches, or a parser_result by estimating false discovery rate.

    Validate a list of crosslinks or crosslink-spectrum-matches, or a parser_result by estimating false discovery rate (FDR) using the defined formula.
    Requires that "score", "alpha_decoy" and "beta_decoy" fields are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        A list of crosslink-spectrum-matches or crosslinks to validate, or a parser_result.
    fdr : float, default = 0.01
        The target FDR, must be given as a real number between 0 and 1. The default of 0.01 corresponds to 1% FDR.
    formula : str, one of "D/T", "(TD+DD)/TT", or "(TD-DD)/TT", default = "D/T"
        Which formula to use to estimate FDR. 'D' denotes any decoy matches including decoy-decoy, decoy-target and target-decoy matches.
        'DD' denotes decoy-decoy matches, 'T' and 'TT' denote target-target matches, and 'TD' denotes target-decoy
        and decoy-target matches.
    score : str, one of "higher_better" or "lower_better", default = "higher_better"
        If a higher score is considered better, or a lower score is considered better.
    separate_intra_inter : bool, default = False
        If FDR should be estimated separately for intra and inter matches.
    ignore_missing_labels : bool, default = False
        If crosslinks and crosslink-spectrum-matches should be ignored if they don't have target and decoy labels. By default and error is
        thrown if any unlabelled data is encountered.

    Returns
    -------
    list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        If a list of crosslink-spectrum-matches or crosslinks was provided, a list of validated
        crosslink-spectrum-matches or crosslinks is returned. If a parser_result was provided,
        an parser_result with validated crosslink-spectrum-matches and/or validated crosslinks will
        be returned.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If parameter formula is not one of 'D/T', '(TD+DD)/TT', or '(TD-DD)/TT'.
    TypeError
        If parameter score is not one of 'higher_better' or 'lower_better'.
    ValueError
        If parameter fdr is outside of the supported range.
    ValueError
        If the number of DD matches exceeds the number of TD matches for formula '(TD-DD)/TT'.
        FDR cannot be estimated with the formula '(TD-DD)/TT' in these cases.

    Warnings
    --------
    Please note that the by default selected FDR formula of 'D/T' actually denotes '(TD+DD)/TT' where 'TD' includes both 'TD' and 'DT' matches
    since we consider any peptide pair with at least one decoy hit a decoy. Estimating FDR via 'DD/TT' is not valid for crosslinking mass spectrometry
    as it severly underestimated the actual error and is therefore not supported by pyXLMS!

    Notes
    -----
    Please note that progress bars will usually not complete when running this function. This is by design as it is not
    necessary to iterate over all scores to estimate FDR.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import validate
    >>> pr = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> len(csms)
    826
    >>> validated = validate(csms)
    >>> len(validated)
    705

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import validate
    >>> pr = read(
    ...     [
    ...         "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...         "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
    ...     ],
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> len(pr["crosslink-spectrum-matches"])
    826
    >>> len(pr["crosslinks"])
    300
    >>> validated = validate(pr)
    >>> len(validated["crosslink-spectrum-matches"])
    705
    >>> len(validated["crosslinks"])
    226

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import validate
    >>> pr = read(
    ...     [
    ...         "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...         "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
    ...     ],
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> len(pr["crosslink-spectrum-matches"])
    826
    >>> len(pr["crosslinks"])
    300
    >>> validated = validate(pr, fdr=0.05)
    >>> len(validated["crosslink-spectrum-matches"])
    825
    >>> len(validated["crosslinks"])
    260
    """
    _ok = check_input_multi(data, "data", [ParserResult, list])
    _ok = check_input(fdr, "fdr", float)
    _ok = check_input(formula, "formula", str)
    _ok = check_input(score, "score", str)
    _ok = check_input(separate_intra_inter, "separate_intra_inter", bool)
    _ok = check_input(ignore_missing_labels, "ignore_missing_labels", bool)
    if fdr >= 1.0 or fdr <= 0.0:
        raise ValueError(
            "FDR must be given as a real number between 0 and 1, e.g. 0.01 corresponds to 1% FDR!"
        )
    if formula not in ["D/T", "(TD+DD)/TT", "(TD-DD)/TT"]:
        raise TypeError(
            "Parameter 'formula' has to be one of 'D/T', '(TD+DD)/TT' or '(TD-DD)/TT'! Where D and DD is the number of decoys, T and TT the number of targets, "
            "and TD the number of target-decoys!"
        )
    if score not in ["higher_better", "lower_better"]:
        raise TypeError(
            "Parameter 'score' has to be one of 'higher_better' or 'lower_better'!"
        )
    if isinstance(data, list):
        csms_or_xls = assert_csms_or_xls(data)
        if len(csms_or_xls) == 0:
            return csms_or_xls
        if ignore_missing_labels:
            csms_or_xls = [
                item
                for item in csms_or_xls
                if item["alpha_decoy"] is not None and item["beta_decoy"] is not None
            ]
        csms_or_xls = assert_csms_or_xls(csms_or_xls)
        _ok = check_available_keys(["score", "alpha_decoy", "beta_decoy"], csms_or_xls)
        if formula == "(TD-DD)/TT":
            if len(filter_target_decoy(csms_or_xls)["Target-Decoy"]) == 0:
                raise ValueError(
                    "Can't estimate FDR with formula '(TD-DD)/TT' when there are no TD matches! Please select the default formula instead!"
                )
            if separate_intra_inter:
                separate = filter_crosslink_type(csms_or_xls)
                return assert_csms_or_xls(
                    __validate_relaxed(separate["Intra"], fdr, score)
                    + __validate_relaxed(separate["Inter"], fdr, score)
                )
            return __validate_relaxed(csms_or_xls, fdr, score)
        if separate_intra_inter:
            separate = filter_crosslink_type(csms_or_xls)
            return assert_csms_or_xls(
                __validate_strict(separate["Intra"], fdr, score)
                + __validate_strict(separate["Inter"], fdr, score)
            )
        return __validate_strict(csms_or_xls, fdr, score)
    new_csms = (
        assert_csms(
            validate(
                data["crosslink-spectrum-matches"],
                fdr,
                formula,
                score,
                separate_intra_inter,
                ignore_missing_labels,
            )
        )
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        assert_xls(
            validate(
                data["crosslinks"],
                fdr,
                formula,
                score,
                separate_intra_inter,
                ignore_missing_labels,
            )
        )
        if data["crosslinks"] is not None
        else None
    )
    return create_parser_result(
        search_engine=data["search_engine"],
        csms=new_csms,
        crosslinks=new_xls,
    )
