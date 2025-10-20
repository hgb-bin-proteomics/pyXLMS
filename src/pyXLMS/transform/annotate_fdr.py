#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
from tqdm import tqdm

from ..data import check_input
from ..data import check_input_multi
from ..data import create_parser_result
from .util import get_available_keys
from .filter import filter_target_decoy
from .filter import filter_crosslink_type

from typing import Dict
from typing import List
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __annotate_fdr_strict(
    data: List[Dict[str, Any]],
    score: Literal["higher_better", "lower_better"],
) -> List[Dict[str, Any]]:
    r"""Annotates a list of crosslinks or crosslink-spectrum-matches with their false discover rate by strict false discovery rate estimation.

    Annotates a list of crosslinks or crosslink-spectrum-matches with their false discover rate (FDR) by strict false discovery rate estimation
    using the formula (TD+DD)/TT. The annotated FDR is written to "additional_information" for every crosslink or crosslink-spectrum-match and
    is available via the accessor "pyXLMS_annotated_FDR".
    Requires that "score", "alpha_decoy" and "beta_decoy" fields are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    score: str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    list of dict of str, any
        A list of FDR annotated crosslink-spectrum-matches or crosslinks.

    Warns
    -----
    RuntimeWarning
        If one of the false discovery rates can't be computed because of not enough target-target matches.

    Notes
    -----
    This function should not be called directly, it is called from ``annotate_fdr()``.
    """
    data_type = (
        "crosslinks"
        if data[0]["data_type"] == "crosslink"
        else "crosslink-spectrum-matches"
    )
    for current_item in tqdm(
        data,
        total=len(data),
        desc=f"Annotating FDR for {data_type}...",
    ):
        tt = 0
        td = 0
        dd = 0
        for item in data:
            if score == "higher_better":
                if item["score"] >= current_item["score"]:
                    if item["alpha_decoy"] and item["beta_decoy"]:
                        dd += 1
                    elif not item["alpha_decoy"] and not item["beta_decoy"]:
                        tt += 1
                    else:
                        td += 1
            else:
                if item["score"] <= current_item["score"]:
                    if item["alpha_decoy"] and item["beta_decoy"]:
                        dd += 1
                    elif not item["alpha_decoy"] and not item["beta_decoy"]:
                        tt += 1
                    else:
                        td += 1
        fdr = float("nan")
        if tt == 0:
            warnings.warn(
                RuntimeWarning(
                    f"Could not calculate FDR {data[0]['data_type']} because the number of target-target matches "
                    + f"at score {current_item['score']} is 0!"
                )
            )
        else:
            fdr = (td + dd) / tt
        if current_item["additional_information"] is None:
            current_item["additional_information"] = {"pyXLMS_annotated_FDR": fdr}
        else:
            current_item["additional_information"]["pyXLMS_annotated_FDR"] = fdr
    return data


def __annotate_fdr_relaxed(
    data: List[Dict[str, Any]],
    score: Literal["higher_better", "lower_better"],
) -> List[Dict[str, Any]]:
    r"""Annotates a list of crosslinks or crosslink-spectrum-matches with their false discover rate by relaxed false discovery rate estimation.

    Annotates a list of crosslinks or crosslink-spectrum-matches with their false discover rate (FDR) by relaxed false discovery rate estimation
    using the formula (TD-DD)/TT. The annotated FDR is written to "additional_information" for every crosslink or crosslink-spectrum-match and
    is available via the accessor "pyXLMS_annotated_FDR".
    Requires that "score", "alpha_decoy" and "beta_decoy" fields are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    score: str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    list of dict of str, any
        A list of FDR annotated crosslink-spectrum-matches or crosslinks.

    Warns
    -----
    RuntimeWarning
        If one of the false discovery rates can't be computed because of not enough target-target matches.
    RuntimeWarning
        If one of the false discovery rates can't be computed because the number of decoy-decoy matches exceeds
        the number of target-decoy matches.

    Notes
    -----
    This function should not be called directly, it is called from ``annotate_fdr()``.
    """
    data_type = (
        "crosslinks"
        if data[0]["data_type"] == "crosslink"
        else "crosslink-spectrum-matches"
    )
    for current_item in tqdm(
        data,
        total=len(data),
        desc=f"Annotating FDR for {data_type}...",
    ):
        tt = 0
        td = 0
        dd = 0
        for item in data:
            if score == "higher_better":
                if item["score"] >= current_item["score"]:
                    if item["alpha_decoy"] and item["beta_decoy"]:
                        dd += 1
                    elif not item["alpha_decoy"] and not item["beta_decoy"]:
                        tt += 1
                    else:
                        td += 1
            else:
                if item["score"] <= current_item["score"]:
                    if item["alpha_decoy"] and item["beta_decoy"]:
                        dd += 1
                    elif not item["alpha_decoy"] and not item["beta_decoy"]:
                        tt += 1
                    else:
                        td += 1
        fdr = float("nan")
        if tt == 0:
            warnings.warn(
                RuntimeWarning(
                    f"Could not calculate FDR {data[0]['data_type']} because the number of target-target matches "
                    + f"at score {current_item['score']} is 0!"
                )
            )
        elif td - dd < 0:
            warnings.warn(
                RuntimeWarning(
                    f"Could not calculate FDR {data[0]['data_type']} because number the of decoy-decoy matches "
                    + f"exceeds the number of target-decoy matches at score {current_item['score']}!"
                )
            )
        else:
            fdr = (td + dd) / tt
        if current_item["additional_information"] is None:
            current_item["additional_information"] = {"pyXLMS_annotated_FDR": fdr}
        else:
            current_item["additional_information"]["pyXLMS_annotated_FDR"] = fdr
    return data


def annotate_fdr(
    data: List[Dict[str, Any]] | Dict[str, Any],
    formula: Literal["D/T", "(TD+DD)/TT", "(TD-DD)/TT"] = "D/T",
    score: Literal["higher_better", "lower_better"] = "higher_better",
    separate_intra_inter: bool = False,
    ignore_missing_labels: bool = False,
) -> List[Dict[str, Any]] | Dict[str, Any]:
    r"""Annotates a list of crosslinks or crosslink-spectrum-matches, or a parser_result with their false dicovery rate by estimating the false discovery rate.

    Annotates a list of crosslinks or crosslink-spectrum-matches, or a parser_result with their false discovery rate (FDR) by estimating false discovery rate
    using the defined formula. The annotated FDR is written to "additional_information" for every crosslink or crosslink-spectrum-match and is available via
    the accessor "pyXLMS_annotated_FDR".
    Requires that "score", "alpha_decoy" and "beta_decoy" fields are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of dict of str, any, or dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate, or a parser_result.
    formula : str, one of "D/T", "(TD+DD)/TT", or "(TD-DD)/TT", default = "D/T"
        Which formula to use to estimate FDR. D and DD denote decoy matches, T and TT denote target matches, and TD denotes target-decoy
        and decoy-target matches.
    score: str, one of "higher_better" or "lower_better", default = "higher_better"
        If a higher score is considered better, or a lower score is considered better.
    separate_intra_inter : bool, default = False
        If FDR should be estimated separately for intra and inter matches.
    ignore_missing_labels : bool, default = False
        If crosslinks and crosslink-spectrum-matches should be ignored if they don't have target and decoy labels. By default and error is
        thrown if any unlabelled data is encountered.

    Returns
    -------
    list of dict of str, any, or dict of str, any
        If a list of crosslink-spectrum-matches or crosslinks was provided, a list of annotated
        crosslink-spectrum-matches or crosslinks is returned. If a parser_result was provided,
        an parser_result with annotated crosslink-spectrum-matches and/or annotated crosslinks will
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
        If attribute 'score' is not available for any of the data.
    ValueError
        If attribute 'alpha_decoy' or 'beta_decoy' is not available for any of the data and parameter ignore_missing_labels
        is set to False.
    ValueError
        If the number of DD matches exceeds the number of TD matches for formula '(TD-DD)/TT'.
        FDR cannot be annotated with the formula '(TD-DD)/TT' in these cases.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import annotate_fdr
    >>> pr = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> csms = annotate_fdr(csms)
    >>> validated_csms = [
    ...     csm
    ...     for csm in csms
    ...     if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ... ]
    >>> len(validated_csms)
    705

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import annotate_fdr
    >>> pr = read(
    ...     [
    ...         "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...         "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
    ...     ],
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> pr = annotate_fdr(pr)
    >>> validated_csms = [
    ...     csm
    ...     for csm in pr["crosslink-spectrum-matches"]
    ...     if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ... ]
    >>> len(validated_csms)
    705
    >>> validated_xls = [
    ...     xl
    ...     for xl in pr["crosslinks"]
    ...     if xl["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ... ]
    >>> len(validated_xls)
    226
    """
    _ok = check_input_multi(data, "data", [dict, list])
    _ok = check_input(formula, "formula", str)
    _ok = check_input(score, "score", str)
    _ok = check_input(separate_intra_inter, "separate_intra_inter", bool)
    _ok = check_input(ignore_missing_labels, "ignore_missing_labels", bool)
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
        if ignore_missing_labels:
            data = [
                item
                for item in data
                if item["alpha_decoy"] is not None and item["beta_decoy"] is not None
            ]
        available_keys = get_available_keys(data)
        if (
            not available_keys["score"]
            or not available_keys["alpha_decoy"]
            or not available_keys["beta_decoy"]
        ):
            raise ValueError(
                "Can't annotate data if 'score' or target/decoy labels are missing! Selecting 'ignore_missing_labels = True' will ignore crosslinks and crosslink-spectrum-matches "
                "that don't have a valid target/decoy label and filter them out!"
            )
        if formula == "(TD-DD)/TT":
            if len(filter_target_decoy(data)["Target-Decoy"]) == 0:
                raise ValueError(
                    "Can't annotate FDR with formula '(TD-DD)/TT' when there are no TD matches! Please select the default formula instead!"
                )
            if separate_intra_inter:
                separate = filter_crosslink_type(data)
                return __annotate_fdr_relaxed(
                    separate["Intra"], score
                ) + __annotate_fdr_relaxed(separate["Inter"], score)
            return __annotate_fdr_relaxed(data, score)
        if separate_intra_inter:
            separate = filter_crosslink_type(data)
            return __annotate_fdr_strict(
                separate["Intra"], score
            ) + __annotate_fdr_strict(separate["Inter"], score)
        return __annotate_fdr_strict(data, score)
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError("Can't validate dict. Dict has to be a valid 'parser_result'!")
    new_csms = (
        annotate_fdr(
            data["crosslink-spectrum-matches"],
            formula,
            score,
            separate_intra_inter,
            ignore_missing_labels,
        )
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        annotate_fdr(
            data["crosslinks"],
            formula,
            score,
            separate_intra_inter,
            ignore_missing_labels,
        )
        if data["crosslinks"] is not None
        else None
    )
    if new_csms is not None:
        if not isinstance(new_csms, list):
            raise RuntimeError(
                "Something went wrong while annotating crosslink-spectrum-matches.\n"
                f"Expected data type: list. Got: {type(new_csms)}."
            )
    if new_xls is not None:
        if not isinstance(new_xls, list):
            raise RuntimeError(
                "Something went wrong while annotating crosslinks.\n"
                f"Expected data type: list. Got: {type(new_xls)}."
            )
    return create_parser_result(
        search_engine=data["search_engine"], csms=new_csms, crosslinks=new_xls
    )
