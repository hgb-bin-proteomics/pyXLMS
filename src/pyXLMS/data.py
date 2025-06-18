#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any


def check_input(
    parameter: Any,
    parameter_name: str,
    supported_class: Any,
    supported_subclass: Optional[Any] = None,
) -> bool:
    r"""Checks if the given parameter is of the specified type.

    Function that checks if a given parameter is of the specified type and if iterable, all elements are of the specified element type.
    This is mostly an input check function to catch any errors arising from not supported inputs early.

    Parameters
    ----------
    parameter : any
        Parameter to check class of.
    parameter_name : str
        Name of the parameter.
    supported_class : any
        Class the parameter has to be of.
    supported_subclass : any, or None, default = None
        Class of the values in case the parameter is a list or dict.

    Returns
    -------
    bool
        If the given input is okay.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.

    Examples
    --------
    >>> from pyXLMS.data import check_input
    >>> check_input("PEPTIDE", "peptide_a", str)
    True

    >>> from pyXLMS.data import check_input
    >>> check_input([1, 2], "xl_position_proteins_a", list, int)
    True
    """
    if not isinstance(parameter, supported_class):
        raise TypeError(f"{parameter_name} must be {supported_class}!")
    if isinstance(parameter, list) and supported_subclass is not None:
        for value in parameter:
            if not isinstance(value, supported_subclass):
                raise TypeError(
                    f"List values of {parameter_name} must be {supported_subclass}!"
                )
    if isinstance(parameter, dict) and supported_subclass is not None:
        for key in parameter:
            if not isinstance(parameter[key], supported_subclass):
                raise TypeError(
                    f"Dict values of {parameter_name} must be {supported_subclass}!"
                )
    return True


def check_input_multi(
    parameter: Any,
    parameter_name: str,
    supported_classes: List[Any],
    supported_subclass: Optional[Any] = None,
) -> bool:
    r"""Checks if the given parameter is of one of the specified types.

    Function that checks if a given parameter is of one of the specified types and if iterable, all elements are of the specified element type.
    This is mostly an input check function to catch any errors arising from not supported inputs early.

    Parameters
    ----------
    parameter : any
        Parameter to check class of.
    parameter_name : str
        Name of the parameter.
    supported_class : list of any
        Classes the parameter has to be of.
    supported_subclass : any, or None, default = None
        Class of the values in case the parameter is a list or dict.

    Returns
    -------
    bool
        If the given input is okay.

    Raises
    ------
    TypeError
        If the parameter is not of one of the given classes.

    Examples
    --------
    >>> from pyXLMS.data import check_input_multi
    >>> check_input_multi("PEPTIDE", "peptide_a", [str, list])
    True
    """
    if not isinstance(parameter, tuple(supported_classes)):
        raise TypeError(
            f"{parameter_name} must be one of {','.join([str(c) for c in supported_classes])}!"
        )
    if isinstance(parameter, list) and supported_subclass is not None:
        for value in parameter:
            if not isinstance(value, supported_subclass):
                raise TypeError(
                    f"List values of {parameter_name} must be {supported_subclass}!"
                )
    if isinstance(parameter, dict) and supported_subclass is not None:
        for key in parameter:
            if not isinstance(parameter[key], supported_subclass):
                raise TypeError(
                    f"Dict values of {parameter_name} must be {supported_subclass}!"
                )
    return True


def check_indexing(value: int | List[int]) -> bool:
    r"""Checks that the given value is not 0-based.

    Parameters
    ----------
    value : int, or list of int
        The value(s) to check.

    Returns
    -------
    bool
        If the given value(s) is/are okay.

    Raises
    ------
    ValueError
        If any of the values are smaller than one.

    Examples
    --------
    >>> from pyXLMS.data import check_indexing
    >>> check_indexing([1, 2, 3])
    True
    """
    check_input_multi(value, "value", [int, list], int)
    if isinstance(value, int):
        if value < 1:
            raise ValueError(
                "0-based value found! All positions must use 1-based indexing!"
            )
    else:
        for val in value:
            if val < 1:
                raise ValueError(
                    "0-based value found! All positions must use 1-based indexing!"
                )
    return True


def create_crosslink(
    peptide_a: str,
    xl_position_peptide_a: int,
    proteins_a: Optional[List[str]],
    xl_position_proteins_a: Optional[List[int]],
    decoy_a: Optional[bool],
    peptide_b: str,
    xl_position_peptide_b: int,
    proteins_b: Optional[List[str]],
    xl_position_proteins_b: Optional[List[int]],
    decoy_b: Optional[bool],
    score: Optional[float],
    additional_information: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    r"""Creates a crosslink data structure.

    Contains minimal data necessary for representing a single crosslink. The returned crosslink data structure is a dictionary with keys
    as detailed in the return section.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    proteins_a : list of str, or None
        The accessions of proteins that the first peptide is associated with.
    xl_position_proteins_a : list of int, or None
        Positions of the crosslink in the proteins of the first peptide (1-based).
    decoy_a : bool, or None
        Whether the alpha peptide is from the decoy database or not.
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    proteins_b : list of str, or None
        The accessions of proteins that the second peptide is associated with.
    xl_position_proteins_b : list of int, or None
        Positions of the crosslink in the proteins of the second peptide (1-based).
    decoy_b : bool, or None
        Whether the beta peptide is from the decoy database or not.
    score: float, or None
        Score of the crosslink.
    additional_information: dict with str keys, or None, default = None
        A dictionary with additional information associated with the crosslink.

    Returns
    -------
    dict
        The dictionary representing the crosslink with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_peptide_crosslink_position``,
        ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_decoy``, ``beta_peptide``, ``beta_peptide_crosslink_position``,
        ``beta_proteins``, ``beta_proteins_crosslink_positions``, ``beta_decoy``, ``crosslink_type``, ``score``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.
    ValueError
        If the length of crosslink positions is not equal to the length of proteins.

    Notes
    -----
    The minimum required data for creating a crosslink is:

    - ``peptide_a``: The unmodified amino acid sequence of the first peptide.
    - ``peptide_b``: The unmodified amino acid sequence of the second peptide.
    - ``xl_position_peptide_a``: The position of the crosslinker in the sequence of the first peptide (1-based).
    - ``xl_position_peptide_b``: The position of the crosslinker in the sequence of the second peptide (1-based).

    Examples
    --------
    >>> from pyXLMS.data import create_crosslink
    >>> minimal_crosslink = create_crosslink("PEPTIDEA", 1, None, None, None, "PEPTIDEB", 5, None, None, None, None)
    >>> crosslink = create_crosslink("PEPTIDEA", 1, ["PROTEINA"], [1], False, "PEPTIDEB", 5, ["PROTEINB"], [3], False, 34.5)
    """
    ## input checks
    full = check_input(peptide_a, "peptide_a", str)
    full = check_input(peptide_b, "peptide_b", str)
    full = check_input(xl_position_peptide_a, "xl_position_peptide_a", int)
    full = check_input(xl_position_peptide_b, "xl_position_peptide_b", int)
    full = (
        full and check_input(proteins_a, "proteins_a", list, str)
        if proteins_a is not None
        else False
    )
    full = (
        full and check_input(proteins_b, "proteins_b", list, str)
        if proteins_b is not None
        else False
    )
    full = (
        full
        and check_input(xl_position_proteins_a, "xl_position_proteins_a", list, int)
        if xl_position_proteins_a is not None
        else False
    )
    full = (
        full
        and check_input(xl_position_proteins_b, "xl_position_proteins_b", list, int)
        if xl_position_proteins_b is not None
        else False
    )
    full = (
        full and check_input(decoy_a, "decoy_a", bool) if decoy_a is not None else False
    )
    full = (
        full and check_input(decoy_b, "decoy_b", bool) if decoy_b is not None else False
    )
    full = full and check_input(score, "score", float) if score is not None else False
    if proteins_a is not None and xl_position_proteins_a is not None:
        if len(proteins_a) != len(xl_position_proteins_a):
            raise ValueError(
                "Crosslink position has to be given for every protein! Length of proteins_a and xl_position_proteins_a has to match!"
            )
    if proteins_b is not None and xl_position_proteins_b is not None:
        if len(proteins_b) != len(xl_position_proteins_b):
            raise ValueError(
                "Crosslink position has to be given for every protein! Length of proteins_b and xl_position_proteins_b has to match!"
            )
    _ok = check_indexing(xl_position_peptide_a)
    _ok = check_indexing(xl_position_peptide_b)
    _ok = (
        check_indexing(xl_position_proteins_a)
        if xl_position_proteins_a is not None
        else True
    )
    _ok = (
        check_indexing(xl_position_proteins_b)
        if xl_position_proteins_b is not None
        else True
    )
    ## processing
    key_a = f"{peptide_a.strip()}{xl_position_peptide_a}"
    key_b = f"{peptide_b.strip()}{xl_position_peptide_b}"
    # if homomeric crosslink
    if key_a == key_b:
        key_a += "_0"
        key_b += "_1"
    crosslink = {
        key_a: {
            "peptide": peptide_a,
            "xl_position_peptide": xl_position_peptide_a,
            "proteins": proteins_a,
            "xl_position_proteins": xl_position_proteins_a,
            "decoy": decoy_a,
        },
        key_b: {
            "peptide": peptide_b,
            "xl_position_peptide": xl_position_peptide_b,
            "proteins": proteins_b,
            "xl_position_proteins": xl_position_proteins_b,
            "decoy": decoy_b,
        },
    }
    keys = sorted(list(crosslink.keys()))
    alpha_proteins = (
        [protein.strip() for protein in crosslink[keys[0]]["proteins"]]
        if crosslink[keys[0]]["proteins"] is not None
        else []
    )
    beta_proteins = (
        [protein.strip() for protein in crosslink[keys[1]]["proteins"]]
        if crosslink[keys[1]]["proteins"] is not None
        else []
    )
    return {
        "data_type": "crosslink",
        "completeness": "full" if full else "partial",
        "alpha_peptide": crosslink[keys[0]]["peptide"].strip(),
        "alpha_peptide_crosslink_position": crosslink[keys[0]]["xl_position_peptide"],
        "alpha_proteins": alpha_proteins if len(alpha_proteins) > 0 else None,
        "alpha_proteins_crosslink_positions": crosslink[keys[0]][
            "xl_position_proteins"
        ],
        "alpha_decoy": crosslink[keys[0]]["decoy"],
        "beta_peptide": crosslink[keys[1]]["peptide"].strip(),
        "beta_peptide_crosslink_position": crosslink[keys[1]]["xl_position_peptide"],
        "beta_proteins": beta_proteins if len(beta_proteins) > 0 else None,
        "beta_proteins_crosslink_positions": crosslink[keys[1]]["xl_position_proteins"],
        "beta_decoy": crosslink[keys[1]]["decoy"],
        "crosslink_type": "intra"
        if len(set(alpha_proteins).intersection(set(beta_proteins))) > 0
        else "inter",
        "score": score if not pd.isna(score) else None,  # pyright: ignore[reportGeneralTypeIssues]
        "additional_information": additional_information,
    }


def create_crosslink_min(
    peptide_a: str,
    xl_position_peptide_a: int,
    peptide_b: str,
    xl_position_peptide_b: int,
    **kwargs,
) -> Dict[str, Any]:
    r"""Creates a crosslink data structure from minimal input.

    Contains minimal data necessary for representing a single crosslink. This is an alias for
    ``data.create_crosslink()``that sets all optional parameters to ``None`` for convenience.
    The returned crosslink data structure is a dictionary with keys as detailed in the return
    section.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    **kwargs
        Any additional parameters will be passed to ``data.create_crosslink()``.

    Returns
    -------
    dict
        The dictionary representing the crosslink with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_peptide_crosslink_position``,
        ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_decoy``, ``beta_peptide``, ``beta_peptide_crosslink_position``,
        ``beta_proteins``, ``beta_proteins_crosslink_positions``, ``beta_decoy``, ``crosslink_type``, ``score``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Notes
    -----
    See also ``data.create_crosslink()``.

    Examples
    --------
    >>> from pyXLMS.data import create_crosslink_min
    >>> minimal_crosslink = create_crosslink_min("PEPTIDEA", 1, "PEPTIDEB", 5)
    """
    return create_crosslink(
        peptide_a=peptide_a,
        xl_position_peptide_a=xl_position_peptide_a,
        proteins_a=kwargs["proteins_a"] if "proteins_a" in kwargs else None,
        xl_position_proteins_a=kwargs["xl_position_proteins_a"]
        if "xl_position_proteins_a" in kwargs
        else None,
        decoy_a=kwargs["decoy_a"] if "decoy_a" in kwargs else None,
        peptide_b=peptide_b,
        xl_position_peptide_b=xl_position_peptide_b,
        proteins_b=kwargs["proteins_b"] if "proteins_b" in kwargs else None,
        xl_position_proteins_b=kwargs["xl_position_proteins_b"]
        if "xl_position_proteins_b" in kwargs
        else None,
        decoy_b=kwargs["decoy_b"] if "decoy_b" in kwargs else None,
        score=kwargs["score"] if "score" in kwargs else None,
        additional_information=kwargs["additional_information"]
        if "additional_information" in kwargs
        else None,
    )


def create_crosslink_from_csm(csm: Dict[str, Any]) -> Dict[str, Any]:
    r"""Creates a crosslink data structure from a crosslink-spectrum-match.

    Creates a crosslink data structure from a crosslink-spectrum-match. The returned crosslink data structure is a dictionary with keys
    as detailed in the return section.

    Parameters
    ----------
    csm : dict of str
        The crosslink-spectrum-match item to be converted to a crosslink item.

    Returns
    -------
    dict
        The dictionary representing the crosslink with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_peptide_crosslink_position``,
        ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_decoy``, ``beta_peptide``, ``beta_peptide_crosslink_position``,
        ``beta_proteins``, ``beta_proteins_crosslink_positions``, ``beta_decoy``, ``crosslink_type``, ``score``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Raises
    ------
    TypeError
        If parameter ``csm`` is not a valid crosslink-spectrum-match.

    Notes
    -----
    See also ``data.create_crosslink()``.

    Examples
    --------
    >>> from pyXLMS.data import create_csm_min, create_crosslink_from_csm
    >>> csm = create_csm_min("PEPTIDEA", 1, "PEPTIDEB", 5, "RUN_1", 1)
    >>> crosslink = create_crosslink_from_csm(csm)
    """
    _ok = check_input(csm, "csm", dict)
    if "data_type" not in csm or csm["data_type"] != "crosslink-spectrum-match":
        raise TypeError("Parameter csm is not a valid crosslink-spectrum-match!")
    return create_crosslink(
        peptide_a=csm["alpha_peptide"],
        xl_position_peptide_a=csm["alpha_peptide_crosslink_position"],
        proteins_a=csm["alpha_proteins"],
        xl_position_proteins_a=csm["alpha_proteins_crosslink_positions"],
        decoy_a=csm["alpha_decoy"],
        peptide_b=csm["beta_peptide"],
        xl_position_peptide_b=csm["beta_peptide_crosslink_position"],
        proteins_b=csm["beta_proteins"],
        xl_position_proteins_b=csm["beta_proteins_crosslink_positions"],
        decoy_b=csm["beta_decoy"],
        score=csm["score"],
        additional_information=csm["additional_information"],
    )


def create_csm(
    peptide_a: str,
    modifications_a: Optional[Dict[int, Tuple[str, float]]],
    xl_position_peptide_a: int,
    proteins_a: Optional[List[str]],
    xl_position_proteins_a: Optional[List[int]],
    pep_position_proteins_a: Optional[List[int]],
    score_a: Optional[float],
    decoy_a: Optional[bool],
    peptide_b: str,
    modifications_b: Optional[Dict[int, Tuple[str, float]]],
    xl_position_peptide_b: int,
    proteins_b: Optional[List[str]],
    xl_position_proteins_b: Optional[List[int]],
    pep_position_proteins_b: Optional[List[int]],
    score_b: Optional[float],
    decoy_b: Optional[bool],
    score: Optional[float],
    spectrum_file: str,
    scan_nr: int,
    charge: Optional[int],
    rt: Optional[float],
    im_cv: Optional[float],
    additional_information: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    r"""Creates a crosslink-spectrum-match data structure.

    Contains minimal data necessary for representing a single crosslink-spectrum-match. The returned crosslink-spectrum-match data structure
    is a dictionary with keys as detailed in the return section.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    modifications_a : dict of [int, tuple], or None
        The modifications of the first peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.
        ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications should be denoted with position ``len(peptide) + 1``.
        If the peptide is not modified an empty dictionary should be given.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    proteins_a : list of str, or None
        The accessions of proteins that the first peptide is associated with.
    xl_position_proteins_a : list of int, or None
        Positions of the crosslink in the proteins of the first peptide (1-based).
    pep_position_proteins_a : list of int, or None
        Positions of the first peptide in the corresponding proteins (1-based).
    score_a : float, or None
        Identification score of the first peptide.
    decoy_a : bool, or None
        Whether the alpha peptide is from the decoy database or not.
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    modifications_b : dict of [int, tuple], or None
        The modifications of the second peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.
        ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications should be denoted with position ``len(peptide) + 1``.
        If the peptide is not modified an empty dictionary should be given.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    proteins_b : list of str, or None
        The accessions of proteins that the second peptide is associated with.
    xl_position_proteins_b : list of int, or None
        Positions of the crosslink in the proteins of the second peptide (1-based).
    pep_position_proteins_b : list of int, or None
        Positions of the second peptide in the corresponding proteins (1-based).
    score_b : float, or None
        Identification score of the second peptide.
    decoy_b : bool, or None
        Whether the beta peptide is from the decoy database or not.
    score: float, or None
        Score of the crosslink-spectrum-match.
    spectrum_file : str
        Name of the spectrum file the crosslink-spectrum-match was identified in.
    scan_nr : int
        The corresponding scan number of the crosslink-spectrum-match.
    charge : int, or None
        The precursor charge of the corresponding mass spectrum of the crosslink-spectrum-match.
    rt : float, or None
        The retention time of the corresponding mass spectrum of the crosslink-spectrum-match in seconds.
    im_cv : float, or None
        The ion mobility or compensation voltage of the corresponding mass spectrum of the crosslink-spectrum-match.
    additional_information: dict with str keys, or None, default = None
        A dictionary with additional information associated with the crosslink-spectrum-match.

    Returns
    -------
    dict
        The dictionary representing the crosslink-spectrum-match with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_modifications``,
        ``alpha_peptide_crosslink_position``, ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_proteins_peptide_positions``,
        ``alpha_score``, ``alpha_decoy``, ``beta_peptide``, ``beta_modifications``, ``beta_peptide_crosslink_position``, ``beta_proteins``,
        ``beta_proteins_crosslink_positions``, ``beta_proteins_peptide_positions``, ``beta_score``, ``beta_decoy``, ``crosslink_type``, ``score``,
        ``spectrum_file``, ``scan_nr``, ``retention_time``, ``ion_mobility``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.
    ValueError
        If the length of crosslink positions or peptide positions is not equal to the length of proteins.

    Notes
    -----
    The minimum required data for creating a crosslink-spectrum-match is:

    - ``peptide_a``: The unmodified amino acid sequence of the first peptide.
    - ``peptide_b``: The unmodified amino acid sequence of the second peptide.
    - ``xl_position_peptide_a``: The position of the crosslinker in the sequence of the first peptide (1-based).
    - ``xl_position_peptide_b``: The position of the crosslinker in the sequence of the second peptide (1-based).
    - ``spectrum_file``: Name of the spectrum file the crosslink-spectrum-match was identified in.
    - ``scan_nr``: The corresponding scan number of the crosslink-spectrum-match.

    Examples
    --------
    >>> from pyXLMS.data import create_csm
    >>> minimal_csm = create_csm("PEPTIDEA", {}, 1, None, None, None, None, None, "PEPTIDEB", {}, 5, None, None, None, None, None, None, "MS_EXP1", 1, None, None, None)
    >>> csm = create_csm("PEPTIDEA", {1: ("Oxidation", 15.994915)}, 1, ["PROTEINA"], [1], [1], 20.1, False, "PEPTIDEB", {}, 5, ["PROTEINB"], [3], [1], 33.7, False, 20.1, "MS_EXP1", 1, 3, 13.5, -50)
    """
    ## input checks
    full = check_input(peptide_a, "peptide_a", str)
    full = check_input(peptide_b, "peptide_b", str)
    full = check_input(xl_position_peptide_a, "xl_position_peptide_a", int)
    full = check_input(xl_position_peptide_b, "xl_position_peptide_b", int)
    full = (
        full and check_input(modifications_a, "modifications_a", dict, tuple)
        if modifications_a is not None
        else False
    )
    full = (
        full and check_input(modifications_b, "modifications_b", dict, tuple)
        if modifications_b is not None
        else False
    )
    full = (
        full and check_input(proteins_a, "proteins_a", list, str)
        if proteins_a is not None
        else False
    )
    full = (
        full and check_input(proteins_b, "proteins_b", list, str)
        if proteins_b is not None
        else False
    )
    full = (
        full
        and check_input(xl_position_proteins_a, "xl_position_proteins_a", list, int)
        if xl_position_proteins_a is not None
        else False
    )
    full = (
        full
        and check_input(xl_position_proteins_b, "xl_position_proteins_b", list, int)
        if xl_position_proteins_b is not None
        else False
    )
    full = (
        full
        and check_input(pep_position_proteins_a, "pep_position_proteins_a", list, int)
        if pep_position_proteins_a is not None
        else False
    )
    full = (
        full
        and check_input(pep_position_proteins_b, "pep_position_proteins_b", list, int)
        if pep_position_proteins_b is not None
        else False
    )
    full = (
        full and check_input(score_a, "score_a", float)
        if score_a is not None
        else False
    )
    full = (
        full and check_input(score_b, "score_b", float)
        if score_b is not None
        else False
    )
    full = (
        full and check_input(decoy_a, "decoy_a", bool) if decoy_a is not None else False
    )
    full = (
        full and check_input(decoy_b, "decoy_b", bool) if decoy_b is not None else False
    )
    full = full and check_input(score, "score", float) if score is not None else False
    full = full and check_input(spectrum_file, "spectrum_file", str)
    full = full and check_input(scan_nr, "scan_nr", int)
    full = full and check_input(charge, "charge", int) if charge is not None else False
    full = full and check_input(rt, "rt", float) if rt is not None else False
    full = full and check_input(im_cv, "im_cv", float) if im_cv is not None else False
    if proteins_a is not None and xl_position_proteins_a is not None:
        if len(proteins_a) != len(xl_position_proteins_a):
            raise ValueError(
                "Crosslink position has to be given for every protein! Length of proteins_a and xl_position_proteins_a has to match!"
            )
    if proteins_b is not None and xl_position_proteins_b is not None:
        if len(proteins_b) != len(xl_position_proteins_b):
            raise ValueError(
                "Crosslink position has to be given for every protein! Length of proteins_b and xl_position_proteins_b has to match!"
            )
    if proteins_a is not None and pep_position_proteins_a is not None:
        if len(proteins_a) != len(pep_position_proteins_a):
            raise ValueError(
                "Peptide position has to be given for every protein! Length of proteins_a and pep_position_proteins_a has to match!"
            )
    if proteins_b is not None and pep_position_proteins_b is not None:
        if len(proteins_b) != len(pep_position_proteins_b):
            raise ValueError(
                "Peptide position has to be given for every protein! Length of proteins_b and pep_position_proteins_b has to match!"
            )
    _ok = check_indexing(xl_position_peptide_a)
    _ok = check_indexing(xl_position_peptide_b)
    _ok = (
        check_indexing(xl_position_proteins_a)
        if xl_position_proteins_a is not None
        else True
    )
    _ok = (
        check_indexing(xl_position_proteins_b)
        if xl_position_proteins_b is not None
        else True
    )
    _ok = (
        check_indexing(pep_position_proteins_a)
        if pep_position_proteins_a is not None
        else True
    )
    _ok = (
        check_indexing(pep_position_proteins_b)
        if pep_position_proteins_b is not None
        else True
    )
    ## validity
    if xl_position_proteins_a is not None and pep_position_proteins_a is not None:
        for i in range(len(xl_position_proteins_a)):
            if (
                xl_position_proteins_a[i] - pep_position_proteins_a[i] + 1
                != xl_position_peptide_a
            ):
                _ok = check_indexing(0)
    if xl_position_proteins_b is not None and pep_position_proteins_b is not None:
        for i in range(len(xl_position_proteins_b)):
            if (
                xl_position_proteins_b[i] - pep_position_proteins_b[i] + 1
                != xl_position_peptide_b
            ):
                _ok = check_indexing(0)
    ## processing
    key_a = f"{peptide_a.strip()}{xl_position_peptide_a}"
    key_b = f"{peptide_b.strip()}{xl_position_peptide_b}"
    # if homomeric crosslink
    if key_a == key_b:
        key_a += "_0"
        key_b += "_1"
    crosslink = {
        key_a: {
            "peptide": peptide_a,
            "modifications": {
                int(key): (
                    modifications_a[key][0].strip(),
                    float(modifications_a[key][1]),
                )
                for key in modifications_a.keys()
            }
            if modifications_a is not None
            else None,
            "xl_position_peptide": xl_position_peptide_a,
            "proteins": proteins_a,
            "xl_position_proteins": xl_position_proteins_a,
            "pep_position_proteins": pep_position_proteins_a,
            "score": score_a,
            "decoy": decoy_a,
        },
        key_b: {
            "peptide": peptide_b,
            "modifications": {
                int(key): (
                    modifications_b[key][0].strip(),
                    float(modifications_b[key][1]),
                )
                for key in modifications_b.keys()
            }
            if modifications_b is not None
            else None,
            "xl_position_peptide": xl_position_peptide_b,
            "proteins": proteins_b,
            "xl_position_proteins": xl_position_proteins_b,
            "pep_position_proteins": pep_position_proteins_b,
            "score": score_b,
            "decoy": decoy_b,
        },
    }
    keys = sorted(list(crosslink.keys()))
    alpha_proteins = (
        [protein.strip() for protein in crosslink[keys[0]]["proteins"]]
        if crosslink[keys[0]]["proteins"] is not None
        else []
    )
    beta_proteins = (
        [protein.strip() for protein in crosslink[keys[1]]["proteins"]]
        if crosslink[keys[1]]["proteins"] is not None
        else []
    )
    return {
        "data_type": "crosslink-spectrum-match",
        "completeness": "full" if full else "partial",
        "alpha_peptide": crosslink[keys[0]]["peptide"].strip(),
        "alpha_modifications": crosslink[keys[0]]["modifications"],
        "alpha_peptide_crosslink_position": crosslink[keys[0]]["xl_position_peptide"],
        "alpha_proteins": alpha_proteins if len(alpha_proteins) > 0 else None,
        "alpha_proteins_crosslink_positions": crosslink[keys[0]][
            "xl_position_proteins"
        ],
        "alpha_proteins_peptide_positions": crosslink[keys[0]]["pep_position_proteins"],
        "alpha_score": crosslink[keys[0]]["score"]
        if not pd.isna(crosslink[keys[0]]["score"])
        else None,  # pyright: ignore[reportGeneralTypeIssues]
        "alpha_decoy": crosslink[keys[0]]["decoy"],
        "beta_peptide": crosslink[keys[1]]["peptide"].strip(),
        "beta_modifications": crosslink[keys[1]]["modifications"],
        "beta_peptide_crosslink_position": crosslink[keys[1]]["xl_position_peptide"],
        "beta_proteins": beta_proteins if len(beta_proteins) > 0 else None,
        "beta_proteins_crosslink_positions": crosslink[keys[1]]["xl_position_proteins"],
        "beta_proteins_peptide_positions": crosslink[keys[1]]["pep_position_proteins"],
        "beta_score": crosslink[keys[1]]["score"]
        if not pd.isna(crosslink[keys[1]]["score"])
        else None,  # pyright: ignore[reportGeneralTypeIssues]
        "beta_decoy": crosslink[keys[1]]["decoy"],
        "crosslink_type": "intra"
        if len(set(alpha_proteins).intersection(set(beta_proteins))) > 0
        else "inter",
        "score": score if not pd.isna(score) else None,  # pyright: ignore[reportGeneralTypeIssues]
        "spectrum_file": spectrum_file.strip(),
        "scan_nr": scan_nr,
        "charge": charge,
        "retention_time": rt if not pd.isna(rt) else None,  # pyright: ignore[reportGeneralTypeIssues]
        "ion_mobility": im_cv if not pd.isna(im_cv) else None,  # pyright: ignore[reportGeneralTypeIssues]
        "additional_information": additional_information,
    }


def create_csm_min(
    peptide_a: str,
    xl_position_peptide_a: int,
    peptide_b: str,
    xl_position_peptide_b: int,
    spectrum_file: str,
    scan_nr: int,
    **kwargs,
) -> Dict[str, Any]:
    r"""Creates a crosslink-spectrum-match data structure from minimal input.

    Contains minimal data necessary for representing a single crosslink-spectrum-match. This
    is an alias for ``data.create_csm()``that sets all optional parameters to ``None`` for convenience.
    The returned crosslink-spectrum-match data structure is a dictionary with keys as detailed in the
    return section.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    spectrum_file : str
        Name of the spectrum file the crosslink-spectrum-match was identified in.
    scan_nr : int
        The corresponding scan number of the crosslink-spectrum-match.
    **kwargs
        Any additional parameters will be passed to ``data.create_csm()``.

    Returns
    -------
    dict
        The dictionary representing the crosslink-spectrum-match with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_modifications``,
        ``alpha_peptide_crosslink_position``, ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_proteins_peptide_positions``,
        ``alpha_score``, ``alpha_decoy``, ``beta_peptide``, ``beta_modifications``, ``beta_peptide_crosslink_position``, ``beta_proteins``,
        ``beta_proteins_crosslink_positions``, ``beta_proteins_peptide_positions``, ``beta_score``, ``beta_decoy``, ``crosslink_type``, ``score``,
        ``spectrum_file``, ``scan_nr``, ``retention_time``, ``ion_mobility``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Notes
    -----
    See also ``data.create_csm()``.

    Examples
    --------
    >>> from pyXLMS.data import create_csm_min
    >>> minimal_csm = create_csm("PEPTIDEA", 1, "PEPTIDEB", 5, "MS_EXP1", 1)
    """
    return create_csm(
        peptide_a=peptide_a,
        modifications_a=kwargs["modifications_a"]
        if "modifications_a" in kwargs
        else None,
        xl_position_peptide_a=xl_position_peptide_a,
        proteins_a=kwargs["proteins_a"] if "proteins_a" in kwargs else None,
        xl_position_proteins_a=kwargs["xl_position_proteins_a"]
        if "xl_position_proteins_a" in kwargs
        else None,
        pep_position_proteins_a=kwargs["pep_position_proteins_a"]
        if "pep_position_proteins_a" in kwargs
        else None,
        score_a=kwargs["score_a"] if "score_a" in kwargs else None,
        decoy_a=kwargs["decoy_a"] if "decoy_a" in kwargs else None,
        peptide_b=peptide_b,
        modifications_b=kwargs["modifications_b"]
        if "modifications_b" in kwargs
        else None,
        xl_position_peptide_b=xl_position_peptide_b,
        proteins_b=kwargs["proteins_b"] if "proteins_b" in kwargs else None,
        xl_position_proteins_b=kwargs["xl_position_proteins_b"]
        if "xl_position_proteins_b" in kwargs
        else None,
        pep_position_proteins_b=kwargs["pep_position_proteins_b"]
        if "pep_position_proteins_b" in kwargs
        else None,
        score_b=kwargs["score_b"] if "score_b" in kwargs else None,
        decoy_b=kwargs["decoy_b"] if "decoy_b" in kwargs else None,
        score=kwargs["score"] if "score" in kwargs else None,
        spectrum_file=spectrum_file,
        scan_nr=scan_nr,
        charge=kwargs["charge"] if "charge" in kwargs else None,
        rt=kwargs["rt"] if "rt" in kwargs else None,
        im_cv=kwargs["im_cv"] if "im_cv" in kwargs else None,
        additional_information=kwargs["additional_information"]
        if "additional_information" in kwargs
        else None,
    )


def create_parser_result(
    search_engine: str,
    csms: Optional[List[Dict[str, Any]]],
    crosslinks: Optional[List[Dict[str, Any]]],
) -> Dict[str, Any]:
    r"""Creates a parser result data structure.

    Contains all necessary data elements that should be contained in a result returned by a crosslink search engine result parser.

    Parameters
    ----------
    search_engine : str
        Name of the identifying crosslink search engine.
    csms : list of dict, or None
        List of crosslink-spectrum-matches as created by ``data.create_csm()``.
    crosslinks : list of dict, or None
        List of crosslinks as created by ``data.create_crosslink()``.

    Returns
    -------
    dict
        The parser result data structure which is a dictionary with keys ``data_type``, ``completeness``, ``search_engine``, ``crosslink-spectrum-matches`` and
        ``crosslinks``.

    Examples
    --------
    >>> from pyXLMS.data import create_parser_result
    >>> result = create_parser_result("MS Annika", None, None)
    >>> result["data_type"]
    'parser_result'
    >>> result["completeness"]
    'empty'
    >>> result["search_engine"]
    'MS Annika'
    """
    completeness = "partial"
    if csms is None and crosslinks is None:
        completeness = "empty"
    if csms is not None and crosslinks is not None:
        completeness = "full"
    return {
        "data_type": "parser_result",
        "completeness": completeness,
        "search_engine": search_engine,
        "crosslink-spectrum-matches": csms,
        "crosslinks": crosslinks,
    }
