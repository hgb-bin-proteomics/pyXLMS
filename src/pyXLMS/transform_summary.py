#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from .data import check_input
from .data import check_input_multi

from typing import Dict
from typing import List
from typing import Any


def __summary_csm(data: List[Dict[str, Any]], **kwargs) -> Dict[str, Any]:
    # todo
    # nr
    # nr unique
    # nr intra
    # nr inter
    # nr tt
    # nr td
    # nr dd
    # min score
    # max score
    return {}


def __summary_xl(data: List[Dict[str, Any]], **kwargs) -> Dict[str, Any]:
    # todo
    # nr
    # nr unique by seq
    # nr unique by pos
    # nr intra
    # nr inter
    # nr tt
    # nr td
    # nr dd
    # min score
    # max score
    return {}


def summary(data: List[Dict[str, Any]] | Dict[str, Any], **kwargs) -> Dict[str, Any]:
    r""" """
    _ok = check_input_multi(data, "data", [dict, list])
    if isinstance(data, list):
        _ok = check_input(data, "data", list, dict)
        if "data_type" not in data[0] or data[0]["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match, "
                "or a parser_result!"
            )
        if data[0]["data_type"] == "crosslink-spectrum-match":
            return __summary_csm(data, **kwargs)
        return __summary_xl(data, **kwargs)
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError(
            "Can't annotate positions for dict. Dict has to be a valid 'parser_result'!"
        )
    csm_summary = (
        __summary_csm(data["crosslink-spectrum-matches"], **kwargs)
        if data["crosslink-spectrum-matches"] is not None
        else {}
    )
    xl_summary = (
        __summary_xl(data["crosslinks"], **kwargs)
        if data["crosslinks"] is not None
        else {}
    )
    return {**csm_summary, **xl_summary}
