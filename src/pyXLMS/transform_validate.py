#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import numpy as np
from tqdm import tqdm

from .data import check_input
from .data import check_input_multi
from .data import create_parser_result
from .transform_util import get_available_keys

from typing import Optional
from typing import Dict
from typing import List
from typing import Tuple
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __verify_fdr_strict(
    data: List[Dict[str, Any]], 
    fdr: float,
    cutoff: float,
    score: Literal["higher_better", "lower_better"]
) -> bool:
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
    data: List[Dict[str, Any]],
    fdr: float,
    score: Literal["higher_better", "lower_better"]
) -> List[Dict[str, Any]]:
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
        td = np.array(td)[np.argsort(scores, stable=True)]
        scores = scores[np.argsort(scores, stable=True)]
        cutoff = scores[0] #scores.max()
    else:
        td = np.array(td)[np.argsort(scores, stable=True)[::-1]]
        scores = scores[np.argsort(scores, stable=True)[::-1]]
        cutoff = scores[0] #scores.min()
    nr_items = len(td)
    for i in tqdm(range(nr_items), total=nr_items, desc="Iterating over scores for FDR calculation..."):
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
    data: List[Dict[str, Any]], 
    fdr: float,
    cutoff: float,
    score: Literal["higher_better", "lower_better"]
) -> bool:
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
    if (DT-D) < 0.0:
        raise RuntimeError()
    return (DT-D) / T < fdr


def __validate_relaxed(
    data: List[Dict[str, Any]],
    fdr: float,
    score: Literal["higher_better", "lower_better"]
) -> List[Dict[str, Any]]:
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
        td = np.array(td)[np.argsort(scores, stable=True)]
        tdd = np.array(tdd)[np.argsort(scores, stable=True)]
        scores = scores[np.argsort(scores, stable=True)]
        cutoff = scores[0] #scores.max()
    else:
        td = np.array(td)[np.argsort(scores, stable=True)[::-1]]
        tdd = np.array(tdd)[np.argsort(scores, stable=True)[::-1]]
        scores = scores[np.argsort(scores, stable=True)[::-1]]
        cutoff = scores[0] #scores.min()
    nr_items = len(td)
    for i in tqdm(range(nr_items), total=nr_items, desc="Iterating over scores for FDR calculation..."):
        if tdd[i:].sum() < 0.0:
            raise RuntimeError()
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
    data: List[Dict[str, Any]] | Dict[str, Any],
    fdr: float = 0.01,
    formula: Literal["D/T", "(TD+DD)/TT", "(TD-DD)/TT"] = "D/T",
    score: Literal["higher_better", "lower_better"] = "higher_better":
    separate_intra_inter: bool = False,
    ignore_missing_labels: bool = False,
) -> List[Dict[str, Any]] | Dict[str, Any]:
    _ok = check_input_multi(data, "data", [dict, list])
    _ok = check_input(fdr, "fdr", float)
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
            "Parameter 'score' has to be one of 'higher_better' or 'lower_better'! If two identical crosslinks or crosslink-spectrum"
            "-matches are found, the one with the higher score is kept if 'higher_better' is selected, and vice versa."
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
            data = [item for item in data if item["alpha_decoy"] is not None and item["beta_decoy"] is not None]
        available_keys = get_available_keys(data)
        if not available_keys["score"] or not available_keys["alpha_decoy"] or not available_keys["beta_decoy"]:
            raise ValueError(
                "Can't validate data if 'score' or target/decoy labels are missing! Selecting 'ignore_missing_labels = True' will ignore crosslinks and crosslink-spectrum-matches "
                "that don't have a valid target/decoy label and filter them out!"
            )
        if formula == "(TD-DD)/TT":
            if separate_intra_inter:
                intra = list()
                inter = list()
                for item in data:
                    if item["crosslink_type"] == "intra":
                        intra.append(item)
                    else:
                        inter.append(item)
                return __validate_relaxed(intra, fdr, score) + __validate_relaxed(inter, fdr, score)
            return __validate_relaxed(data, fdr, score)
        if separate_intra_inter:
            intra = list()
            inter = list()
            for item in data:
                if item["crosslink_type"] == "intra":
                    intra.append(item)
                else:
                    inter.append(item)
            return __validate_strict(intra, fdr, score) + __validate_strict(inter, fdr, score)
        return __validate_strict(data, fdr, score)
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError(
            "Can't validate dict. Dict has to be a valid 'parser_result'!"
        )
    new_csms = (
        validate(data["crosslink-spectrum-matches"], fdr, formula, score, separate_intra_inter, ignore_missing_labels)
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        validate(data["crosslinks"], fdr, formula, score, separate_intra_inter, ignore_missing_labels)
        if data["crosslinks"] is not None
        else None
    )
    if new_csms is not None:
        if not isinstance(new_csms, list):
            raise RuntimeError(
                "Something went wrong while validating crosslink-spectrum-matches.\n"
                f"Expected data type: list. Got: {type(new_csms)}."
            )
    if new_xls is not None:
        if not isinstance(new_xls, list):
            raise RuntimeError(
                "Something went wrong while validating crosslinks.\n"
                f"Expected data type: list. Got: {type(new_xls)}."
            )
    return create_parser_result(
        search_engine=data["search_engine"], csms=new_csms, crosslinks=new_xls
    )
