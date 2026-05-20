#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations


from typing import Optional
from typing import List
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
