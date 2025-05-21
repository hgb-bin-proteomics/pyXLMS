#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from .data import check_input

from typing import Optional
from typing import Dict
from typing import Tuple
from typing import List
from typing import Any


def modifications_to_str(
    modifications: Optional[Dict[int, Tuple[str, float]]],
) -> str | None:
    r"""Returns the string representation of a modifications dictionary.

    Parameters
    ----------
    modifications : dict of [str, tuple], or None
        The modifications of a peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.
        ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications should be denoted with position ``len(peptide) + 1``.

    Returns
    -------
    str, or None
        The string representation of the modifications (or ``None`` if no modification was provided).

    Examples
    --------
    >>> from pyXLMS.transform import modifications_to_str
    >>> modifications_to_str({1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)})
    '(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])'
    """
    ## check input
    _ok = (
        check_input(modifications, "modifcations", dict, tuple)
        if modifications is not None
        else True
    )

    modifications_str = ""
    if modifications is None:
        return None
    for modification_pos in sorted(modifications.keys()):
        modifications_str += f"({modification_pos}:[{modifications[modification_pos][0]}|{modifications[modification_pos][1]}]);"
    return modifications_str.rstrip(";")


def assert_data_type_same(data_list: List[Dict[str, Any]]) -> bool:
    _ok = check_input(data_list, "data_list", list, dict)
    data_type = l[0]["data_type"]
    for item in data_list:
        if item["data_type"] != data_type:
            return False
    return True


def get_available_keys(data_list: List[Dict[str, Any]]) -> Dict[str, bool]:
    if not assert_data_type_same(data_list):
        raise TypeError()
    data_type = data_list[0]["data_type"]
    # available keys
    modifications_a = True
    proteins_a = True,
    xl_position_proteins_a = True,
    pep_position_proteins_a = True,
    score_a = True,
    decoy_a = True,
    modifications_b = True,
    proteins_b = True,
    xl_position_proteins_b = True,
    pep_position_proteins_b = True,
    score_b = True,
    decoy_b = True,
    score = True,
    charge = True,
    rt = True,
    im_cv = True
    additional_information = True
    # parse available keys
    if data_type == "crosslink":
        for data in data_list:
            if data["completeness"] != "full":
                if data["alpha_proteins"] is None:
                    proteins_a = False
                if data["alpha_proteins_crosslink_positions"] is None:
                    xl_position_proteins_a = False
                if data["alpha_decoy"] is None:
                    decoy_a = False
                if data["beta_proteins"] is None:
                    proteins_b = False
                if data["beta_proteins_crosslink_positions"] is None:
                    xl_position_proteins_b = False
                if data["beta_decoy"] is None:
                    decoy_b = False
                if data["score"] is None:
                    score = False
                if data["additional_information"] is None:
                    additional_information = False
        return {
            "data_type": True,
            "completeness": True,
            "alpha_peptide": True,
            "alpha_peptide_crosslink_position": True,
            "alpha_proteins": proteins_a,
            "alpha_proteins_crosslink_positions": xl_position_proteins_a,
            "alpha_decoy": decoy_a,
            "beta_peptide": True,
            "beta_peptide_crosslink_position": True,
            "beta_proteins": proteins_b,
            "beta_proteins_crosslink_positions": xl_position_proteins_b,
            "beta_decoy": decoy_b,
            "crosslink_type": True,
            "score": score,
            "additional_information": additional_information,
        }
    if data_type == "crosslink-spectrum-match":
        pass
    raise TypeError()
    return {"err": True}
