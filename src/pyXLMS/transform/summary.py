#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data import check_input
from ..data import check_input_multi
from .aggregate import unique
from .filter import filter_crosslink_type
from .filter import filter_target_decoy

from typing import Dict
from typing import List
from typing import Any


def __summary_csm(data: List[Dict[str, Any]]) -> Dict[str, float]:
    r"""Extracts summary stats from a list of crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches.

    Returns
    -------
    dict of str, float
        A dictionary with summary stats of the list of crosslink-spectrum-matches.

    Notes
    -----
    This function should not be called directly, it is called from ``summary()``.
    """
    # number of CSMs
    nr = float(len(data))
    # number of unique CSMs
    nr_unique = float(len(unique(data)))
    csm_types = filter_crosslink_type(data)
    # number of intra CSMs
    nr_intra = float(len(csm_types["Intra"]))
    # number of inter CSMs
    nr_inter = float(len(csm_types["Inter"]))
    target_decoys = filter_target_decoy(data)
    # number of TT CSMs
    nr_tt = float(len(target_decoys["Target-Target"]))
    # nr of TD CSMs
    nr_td = float(len(target_decoys["Target-Decoy"]))
    # nr of DD CSMs
    nr_dd = float(len(target_decoys["Decoy-Decoy"]))
    scores = [csm["score"] for csm in data if csm["score"] is not None]
    # min CSM score
    min_score = float(min(scores)) if len(scores) > 0 else float("nan")
    # max CSM score
    max_score = float(max(scores)) if len(scores) > 0 else float("nan")
    return {
        "Number of CSMs": nr,
        "Number of unique CSMs": nr_unique,
        "Number of intra CSMs": nr_intra,
        "Number of inter CSMs": nr_inter,
        "Number of target-target CSMs": nr_tt,
        "Number of target-decoy CSMs": nr_td,
        "Number of decoy-decoy CSMs": nr_dd,
        "Minimum CSM score": min_score,
        "Maximum CSM score": max_score,
    }


def __summary_xl(data: List[Dict[str, Any]]) -> Dict[str, float]:
    r"""Extracts summary stats from a list of crosslinks.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslinks.

    Returns
    -------
    dict of str, float
        A dictionary with summary stats of the list of crosslinks.

    Notes
    -----
    This function should not be called directly, it is called from ``summary()``.
    """
    # number of crosslinks
    nr = float(len(data))
    # number of unique crosslinks by peptide
    nr_unique_peptide = float(len(unique(data, by="peptide")))
    # number of unique crosslinks by protein
    nr_unique_protein = float("nan")
    try:
        nr_unique_protein = float(len(unique(data, by="protein")))
    except Exception as _e:
        pass
    xl_types = filter_crosslink_type(data)
    # number of intra crosslinks
    nr_intra = float(len(xl_types["Intra"]))
    # number of inter crosslinks
    nr_inter = float(len(xl_types["Inter"]))
    target_decoys = filter_target_decoy(data)
    # number of TT crosslinks
    nr_tt = float(len(target_decoys["Target-Target"]))
    # nr of TD crosslinks
    nr_td = float(len(target_decoys["Target-Decoy"]))
    # nr of DD crosslinks
    nr_dd = float(len(target_decoys["Decoy-Decoy"]))
    scores = [xl["score"] for xl in data if xl["score"] is not None]
    # min crosslink score
    min_score = float(min(scores)) if len(scores) > 0 else float("nan")
    # max crosslink score
    max_score = float(max(scores)) if len(scores) > 0 else float("nan")
    return {
        "Number of crosslinks": nr,
        "Number of unique crosslinks by peptide": nr_unique_peptide,
        "Number of unique crosslinks by protein": nr_unique_protein,
        "Number of intra crosslinks": nr_intra,
        "Number of inter crosslinks": nr_inter,
        "Number of target-target crosslinks": nr_tt,
        "Number of target-decoy crosslinks": nr_td,
        "Number of decoy-decoy crosslinks": nr_dd,
        "Minimum crosslink score": min_score,
        "Maximum crosslink score": max_score,
    }


def summary(data: List[Dict[str, Any]] | Dict[str, Any]) -> Dict[str, float]:
    r"""Extracts summary stats from a list of crosslinks or crosslink-spectrum-matches, or a
    parser_result.

    Extracts summary statistics from a list of crosslinks or crosslink-spectrum-matches, or a parser_result.
    The statistic depend on the supplied data type, if a list of crosslinks is supplied a dictionary with the
    following statistics and keys is returned:

    - Number of crosslinks
    - Number of unique crosslinks by peptide
    - Number of unique crosslinks by protein
    - Number of intra crosslinks
    - Number of inter crosslinks
    - Number of target-target crosslinks
    - Number of target-decoy crosslinks
    - Number of decoy-decoy crosslinks
    - Minimum crosslink score
    - Maximum crosslink score

    If a list of crosslink-spectrum-matches is supplied dictionary with the following statistics and keys is
    returned:

    - Number of CSMs
    - Number of unique CSMs
    - Number of intra CSMs
    - Number of inter CSMs
    - Number of target-target CSMs
    - Number of target-decoy CSMs
    - Number of decoy-decoy CSMs
    - Minimum CSM score
    - Maximum CSM score

    If a parser_result is supplied, a dictionary with both containing all of these is returned - if they are available.
    A parser_result that only contains crosslinks will only yield a dictionary with crosslink statistics, and vice versa
    a parser_result that only contains crosslink-spectrum-matches will only yield a dictionary with crosslink-spectrum-
    match statistics. If the parser_result result contains both, then both dictionaries will be merged and returned.
    Please note that in this case a single dictionary is returned, that contains both the keys for crosslinks and
    crosslink-spectrum-matches.

    Statistics are also printed to ``stdout``.

    Parameters
    ----------
    data : list of dict of str, any, or dict of str, any
        A list of crosslinks or crosslink-spectrum-matches, or a parser_result.

    Returns
    -------
    dict of str, float
        A dictionary with summary statistics.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import summary
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> stats = summary(csms)
    Number of CSMs: 826.0
    Number of unique CSMs: 826.0
    Number of intra CSMs: 803.0
    Number of inter CSMs: 23.0
    Number of target-target CSMs: 786.0
    Number of target-decoy CSMs: 39.0
    Number of decoy-decoy CSMs: 1.0
    Minimum CSM score: 1.11
    Maximum CSM score: 452.99

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import summary
    >>> pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx", engine="MS Annika", crosslinker="DSS")
    >>> stats = summary(pr)
    Number of crosslinks: 300.0
    Number of unique crosslinks by peptide: 300.0
    Number of unique crosslinks by protein: 298.0
    Number of intra crosslinks: 279.0
    Number of inter crosslinks: 21.0
    Number of target-target crosslinks: 265.0
    Number of target-decoy crosslinks: 0.0
    Number of decoy-decoy crosslinks: 35.0
    Minimum crosslink score: 1.11
    Maximum crosslink score: 452.99
    """
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
            csm_summary = __summary_csm(data)
            for k, v in csm_summary.items():
                print(f"{k}: {v}")
            return csm_summary
        xl_summary = __summary_xl(data)
        for k, v in xl_summary.items():
            print(f"{k}: {v}")
        return xl_summary
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError(
            "Can't annotate positions for dict. Dict has to be a valid 'parser_result'!"
        )
    csm_summary = (
        __summary_csm(data["crosslink-spectrum-matches"])
        if data["crosslink-spectrum-matches"] is not None
        else {}
    )
    for k, v in csm_summary.items():
        print(f"{k}: {v}")
    xl_summary = (
        __summary_xl(data["crosslinks"]) if data["crosslinks"] is not None else {}
    )
    for k, v in xl_summary.items():
        print(f"{k}: {v}")
    return {**csm_summary, **xl_summary}
