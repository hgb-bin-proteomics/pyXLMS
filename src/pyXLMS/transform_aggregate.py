#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from .data import check_input
from .data import check_input_multi
from .data import create_csm
from .data import create_crosslink
from .data import create_parser_result
from .transform_util import get_available_keys

from typing import Optional
from typing import BinaryIO
from typing import Callable
from typing import Dict
from typing import Tuple
from typing import List
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __score_better(score: float, reference: float, function: Literal["higher_better", "lower_better"]) -> bool:
    r"""Checks if the score is better than the provided reference score.
    
    Checks if the score is better than the provided reference score using the given scoring scheme.
    
    Parameters
    ----------
    score : float
        The score that should be compared.
    reference : float
        The reference score to compare to.
    function : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.
        
    Returns
    -------
    bool
        If the given score is better than the reference score.
    """
    if function == "higher_better":
        return score > reference
    return score < reference


def __get_csm_key(csm: Dict[str, Any]) -> str:
    r"""Get the unique key for a crosslink-spectrum-match.
    
    Parameters
    ----------
    csm : dict of str, any
        A pyXLMS crosslink-spectrum-match object.
        
    Returns
    -------
    str
        The unique key for the crosslink-spectrum-match.
    """
    return f"{csm['spectrum_file']}_{csm['scan_nr']}"


def __get_xl_key(xl: Dict[str, Any], by: Literal["peptide", "protein"]) -> str:
    r"""Get the unique key for a crosslink.
    
    Parameters
    ----------
    xl : dict of str, any
        A pyXLMS crosslink object.
    by : str, one of "peptide" or "protein"
        If peptide or protein crosslink position should be used for determining if a crosslink is unique.
        
    Returns
    -------
    str
        The unique key for the crosslink.
        
    Notes
    -----
    This function should not be called directly, it is called from ``__unique_xls()``.
    """
    if by == "peptide":
        return f"{xl['alpha_peptide']}_{xl['alpha_peptide_crosslink_position']}-{xl['beta_peptide']}_{xl['beta_peptide_crosslink_position']}"
    prot_pos_a = "-".join(sorted([f"{xl['alpha_proteins'][i]}_{xl['alpha_proteins_crosslink_positions'][i]}" for i in range(len(xl["alpha_proteins"]))]))
    prot_pos_b = "-".join(sorted([f"{xl['beta_proteins'][i]}_{xl['beta_proteins_crosslink_positions'][i]}" for i in range(len(xl["beta_proteins"]))]))
    return ":".join(sorted([prot_pos_a, prot_pos_b]))


def __unique_csms(
    csms: List[Dict[str, Any]],
    has_scores: bool,
    score: Literal["higher_better", "lower_better"]
) -> List[Dict[str, Any]]:
    unique_csms = dict()
    for csm in csms:
        key = __get_csm_key(csm)
        if key not in unique_csms:
            unique_csms[key] = csm
        elif has_scores and __score_better(csm["score"], unique_csms[key]["score"], score):
            unique_csms[key] = csm
        else:
            # do nothing
            pass
    return list(unique_csms.values())


def __unique_xls(
    xls: List[Dict[str, Any]],
    by: Literal["peptide", "protein"],
    has_scores: bool,
    score: Literal["higher_better", "lower_better"]
) -> List[Dict[str, Any]]:
    unique_xls = dict()
    for xl in xls:
        key = __get_xl_key(xl, by)
        if key not in unique_xls:
            unique_xls[key] = xl
        elif has_scores and __score_better(xl["score"], unique_xls[key]["score"], score):
            unique_xls[key] = xl
        else:
            # do nothing
            pass
    return list(unique_xls.values())


def unique(
    data: Dict[str, Any] | List[Dict[str, Any]],
    by: Literal["peptide", "protein"] = "protein",
    score: Literal["higher_better", "lower_better"] = "higher_better"
) -> Dict[str, Any] | List[Dict[str, Any]]:
    _ok = check_input_multi(data, "data", [dict, list])
    _ok = check_input(by, "by", str)
    _ok = check_input(score, "score", str)
    if by not in ["peptide", "protein"]:
        raise TypeError(
            "Parameter 'by' has to be one of 'peptide' or 'protein'! Option 'peptide' will group by peptide sequence and "
            "peptide crosslink position while option 'protein' will group by protein identifier and protein crosslink position."
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
        available_keys = get_available_keys(data)
        unique_items = list()
        if data[0]["data_type"] == "crosslink" and by == "protein":
            if available_keys["alpha_proteins"] and available_keys["alpha_proteins_crosslink_positions"] and available_keys["beta_proteins"] and available_keys["beta_proteins_crosslink_positions"]:
                unique_items += __unique_xls(data, by, available_keys["score"], score)
            else:
                raise ValueError(
                    "Grouping by protein crosslink position is only available if all crosslinks have defined protein crosslink positions!\n"
                    "This error might be fixable with 'transform.reannotate_positions()'"!
                )
        elif data[0]["data_type"] == "crosslink":
            unique_items += __unique_xls(data, by, available_keys["score"], score)
        else:
            unique_items += __unique_csms(data, available_keys["score"], score)
        return unique_items
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError(
            "Can't annotate positions for dict. Dict has to be a valid 'parser_result'!"
        )
    new_csms = (
        unique(
            data["crosslink-spectrum-matches"], by, score
        )
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        unique(data["crosslinks"], by, score)
        if data["crosslinks"] is not None
        else None
    )
    if new_csms is not None:
        if not isinstance(new_csms, list):
            raise RuntimeError(
                "Something went wrong while getting unique crosslink-spectrum-matches.\n"
                f"Expected data type: list. Got: {type(new_csms)}."
            )
    if new_xls is not None:
        if not isinstance(new_xls, list):
            raise RuntimeError(
                "Something went wrong while getting unique crosslinks.\n"
                f"Expected data type: list. Got: {type(new_xls)}."
            )
    return create_parser_result(
        search_engine=data["search_engine"], csms=new_csms, crosslinks=new_xls
    )
