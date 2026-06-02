#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from typing import Optional
from typing import Dict
from typing import List
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
    supported_classes : list of any
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


def __get_modified_peptide(
    sequence: str,
    modifications: Optional[Dict[int, Tuple[str, float]]],
    crosslink_position: int,
    crosslinker: Optional[str | float],
) -> str:
    r"""Returns the Proforma string for a single peptide.

    Parameters
    ----------
    sequence : str
        The unmodified peptide sequence.
    modifications : dict of int, tuple of str and float
        The pyXLMS specific modifications object. See ``data.create_csm()`` for reference.
    crosslink_position : int
        Crosslink position in the peptide sequence (1-based).
    crosslinker : str, or float, or None
        Optional name or mass of the crosslink reagent. If the name is given, it should be a valid
        name from XLMOD.

    Returns
    -------
    str
        The Proforma string of the peptidoform.

    Notes
    -----
    - This function should not be called directly, it is called from ``__to_proforma_csm()`` and ``__to_proforma_xl``.
    - Modifications with unknown mass are skipped.
    - If no modifications are given, only the crosslink modification will be encoded in the Proforma.
    - If no modifications are given and no crosslinker is given, the unmodified peptide Proforma will be returned.
    """
    if isinstance(crosslinker, float):
        crosslinker = f"+{crosslinker}" if crosslinker > 0.0 else f"{crosslinker}"
    pep_len = len(sequence)
    if modifications is not None:
        new_modifications = dict()
        for pos, mod in modifications.items():
            if not pd.isna(mod[1]):
                new_modifications[pos] = (
                    mod[0],
                    f"+{mod[1]}" if mod[1] > 0.0 else f"{mod[1]}",
                )
        if crosslink_position not in new_modifications and crosslinker is not None:
            new_modifications[crosslink_position] = ("", f"{crosslinker}")
        for pos in sorted(new_modifications.keys(), reverse=True):
            if pos == 0:
                sequence = f"[{new_modifications[pos][1]}]-" + sequence
            elif pos == pep_len + 1:
                sequence = sequence + f"-[{new_modifications[pos][1]}]"
            else:
                sequence = (
                    sequence[:pos] + f"[{new_modifications[pos][1]}]" + sequence[pos:]
                )
        return sequence
    if crosslinker is not None:
        sequence = (
            sequence[:crosslink_position]
            + f"[{crosslinker}]"
            + sequence[crosslink_position:]
        )
        return sequence
    return sequence
