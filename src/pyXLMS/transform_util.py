#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from .data import check_input

from typing import Optional
from typing import Dict
from typing import Tuple


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
