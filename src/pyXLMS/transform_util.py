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
    peptide_a = True,
    modifications_a = True
    xl_position_peptide_a = True,
    proteins_a = True,
    xl_position_proteins_a = True,
    pep_position_proteins_a = True,
    score_a = True,
    decoy_a = True,
    peptide_b = True,
    modifications_b = True,
    xl_position_peptide_b = True,
    proteins_b = True,
    xl_position_proteins_b = True,
    pep_position_proteins_b = True,
    score_b = True,
    decoy_b = True,
    score = True,
    spectrum_file = True,
    scan_nr = True,
    charge = True,
    rt = True,
    im_cv = True
    # parse available keys
    if data_type == "crosslink":
        pass
    if data_type == "crosslink-spectrum-match":
        pass
    raise TypeError()
    return {"err": True}
