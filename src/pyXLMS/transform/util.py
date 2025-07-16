#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data import check_input

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
    r"""Checks that all data is of the same data type.

    Verifies that all elements in the provided list are of the same data type.

    Parameters
    ----------
    data_list : list of dict of str, any
        A list of dictionaries with the ``data_type`` key.

    Returns
    -------
    bool
        If all elements are of the same data type.

    Examples
    --------
    >>> from pyXLMS.transform import assert_data_type_same
    >>> from pyXLMS import data
    >>> data_list = [data.create_crosslink_min("PEPK", 4, "PKEP", 2), data.create_crosslink_min("KPEP", 1, "PEKP", 3)]
    >>> assert_data_type_same(data_list)
    True

    >>> from pyXLMS.transform import assert_data_type_same
    >>> from pyXLMS import data
    >>> data_list = [data.create_crosslink_min("PEPK", 4, "PKEP", 2), data.create_csm_min("KPEP", 1, "PEKP", 3, "RUN_1", 1)]
    >>> assert_data_type_same(data_list)
    False
    """
    _ok = check_input(data_list, "data_list", list, dict)
    data_type = data_list[0]["data_type"]
    for item in data_list[1:]:
        if item["data_type"] != data_type:
            return False
    return True


def get_available_keys(data_list: List[Dict[str, Any]]) -> Dict[str, bool]:
    r"""Checks which data is available from a list of crosslinks or crosslink-spectrum-matches.

    Verifies which data fields have been set for all crosslinks or crosslink-spectrum-matches in the
    given list. Will return a dictionary structured the same as a crosslink or crosslink-spectrum-match,
    but instead of the data it will return either True or False, depending if the field was set or not.

    Parameters
    ----------
    data_list : list of dict of str, any
        A list of crosslinks or crosslink-spectrum-matches.

    Returns
    -------
    dict of str, bool
        - If a list of crosslinks was provided, a dictionary with the following keys will be returned, where the value
          of each key denotes if the data field is available for all crosslinks in ``data_list``.
          Keys: ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_peptide_crosslink_position``,
          ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_decoy``, ``beta_peptide``, ``beta_peptide_crosslink_position``,
          ``beta_proteins``, ``beta_proteins_crosslink_positions``, ``beta_decoy``, ``crosslink_type``, ``score``, and ``additional_information``.
        - If a list of crosslink-spectrum-matches was provided, a dictionary with the following keys will be returned, where the value
          of each key denotes if the data field is available for all crosslink-spectrum-matches in ``data_list``.
          Keys: ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_modifications``,
          ``alpha_peptide_crosslink_position``, ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_proteins_peptide_positions``,
          ``alpha_score``, ``alpha_decoy``, ``beta_peptide``, ``beta_modifications``, ``beta_peptide_crosslink_position``, ``beta_proteins``,
          ``beta_proteins_crosslink_positions``, ``beta_proteins_peptide_positions``, ``beta_score``, ``beta_decoy``, ``crosslink_type``, ``score``,
          ``spectrum_file``, ``scan_nr``, ``retention_time``, ``ion_mobility``, and ``additional_information``.

    Raises
    ------
    TypeError
        If not all elements in ``data_list`` are of the same data type.
    TypeError
        If one or more elements in the list are of an unsupported data type.

    Examples
    --------
    >>> from pyXLMS.transform import get_available_keys
    >>> from pyXLMS import data
    >>> data_list = [data.create_crosslink_min("PEPK", 4, "PKEP", 2), data.create_crosslink_min("KPEP", 1, "PEKP", 3)]
    >>> available_keys = get_available_keys(data_list)
    >>> available_keys["alpha_peptide"]
    True
    >>> available_keys["score"]
    False
    """
    if not assert_data_type_same(data_list):
        raise TypeError("Not all elements of the list have the same data type!")
    data_type = data_list[0]["data_type"]
    # available keys
    modifications_a = True
    proteins_a = True
    xl_position_proteins_a = True
    pep_position_proteins_a = True
    score_a = True
    decoy_a = True
    modifications_b = True
    proteins_b = True
    xl_position_proteins_b = True
    pep_position_proteins_b = True
    score_b = True
    decoy_b = True
    score = True
    charge = True
    rt = True
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
        for data in data_list:
            if data["completeness"] != "full":
                if data["alpha_modifications"] is None:
                    modifications_a = False
                if data["alpha_proteins"] is None:
                    proteins_a = False
                if data["alpha_proteins_crosslink_positions"] is None:
                    xl_position_proteins_a = False
                if data["alpha_proteins_peptide_positions"] is None:
                    pep_position_proteins_a = False
                if data["alpha_score"] is None:
                    score_a = False
                if data["alpha_decoy"] is None:
                    decoy_a = False
                if data["beta_modifications"] is None:
                    modifications_b = False
                if data["beta_proteins"] is None:
                    proteins_b = False
                if data["beta_proteins_crosslink_positions"] is None:
                    xl_position_proteins_b = False
                if data["beta_proteins_peptide_positions"] is None:
                    pep_position_proteins_b = False
                if data["beta_score"] is None:
                    score_b = False
                if data["beta_decoy"] is None:
                    decoy_b = False
                if data["score"] is None:
                    score = False
                if data["charge"] is None:
                    charge = False
                if data["retention_time"] is None:
                    rt = False
                if data["ion_mobility"] is None:
                    im_cv = False
                if data["additional_information"] is None:
                    additional_information = False
        return {
            "data_type": True,
            "completeness": True,
            "alpha_peptide": True,
            "alpha_modifications": modifications_a,
            "alpha_peptide_crosslink_position": True,
            "alpha_proteins": proteins_a,
            "alpha_proteins_crosslink_positions": xl_position_proteins_a,
            "alpha_proteins_peptide_positions": pep_position_proteins_a,
            "alpha_score": score_a,
            "alpha_decoy": decoy_a,
            "beta_peptide": True,
            "beta_modifications": modifications_b,
            "beta_peptide_crosslink_position": True,
            "beta_proteins": proteins_b,
            "beta_proteins_crosslink_positions": xl_position_proteins_b,
            "beta_proteins_peptide_positions": pep_position_proteins_b,
            "beta_score": score_b,
            "beta_decoy": decoy_b,
            "crosslink_type": True,
            "score": score,
            "spectrum_file": True,
            "scan_nr": True,
            "charge": charge,
            "retention_time": rt,
            "ion_mobility": im_cv,
            "additional_information": additional_information,
        }
    raise TypeError(
        f"Unknown data type {data_type}. Data type must be 'crosslink' or 'crosslink-spectrum-match'!"
    )
    return {"err": True}
