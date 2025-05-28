#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from .data import check_input
from .data import check_input_multi

from typing import Optional
from typing import Dict
from typing import List
from typing import Tuple
from typing import Any


def __get_modified_peptide(
    sequence: str,
    modifications: Optional[Dict[int, Tuple[str, float]]],
    crosslink_position: int,
    crosslinker: Optional[str | float]
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
    - Modifications with unknown mass are skipped.
    - If no modifications are given, only the crosslink modification will be encoded in the Proforma.
    - If no modifications are given and no crosslinker is given, the unmodified peptide Proforma will be returned.
    """
    if isinstance(crosslinker, float):
        crosslinker = f"+{crosslinker}"
    pep_len = len(sequence)
    if modifications is not None:
        modifications = dict(modifications)
        if crosslink_position not in modifications and crosslinker is not None:
            modifications[crosslink_position] = ("", crosslinker)
        for pos in sorted(modifications.keys(), reverse=True):
            if not pd.isna(modifications[pos][1]):
                if pos == 0:
                    sequence = f"[+{modifications[pos][1]}]-" + sequence
                elif pos == pep_len + 1:
                    sequence = sequence + f"-[+{modifications[pos][1]}]"
                else:
                    sequence = sequence[:pos] + f"[+{modifications[pos][1]}]" + sequence[pos:]
        return sequence
    if crosslinker is not None:
        sequence = sequence[:crosslink_position] + f"[{crosslinker}]" + sequence[crosslink_position:]
        return sequence
    return sequence


def __to_proforma_csm(csm: Dict[str, Any], crosslinker: Optional[str | float]) -> str:
    peptide_a = __get_modified_peptide(csm["alpha_peptide"], csm["alpha_modifications"], csm["alpha_peptide_crosslink_position"], crosslinker)
    peptide_b = __get_modified_peptide(csm["beta_peptide"], csm["beta_modifications"], csm["beta_peptide_crosslink_position"], crosslinker)
    if csm["charge"] is not None:
        return f"{peptide_a}//{peptide_b}/{csm['charge']}"
    return f"{peptide_a}//{peptide_b}"


def __to_proforma_xl(xl: Dict[str, Any], crosslinker: Optional[str | float]) -> str:
    peptide_a = __get_modified_peptide(xl["alpha_peptide"], None, xl["alpha_peptide_crosslink_position"], crosslinker)
    peptide_b = __get_modified_peptide(xl["beta_peptide"], None, xl["beta_peptide_crosslink_position"], crosslinker)
    return f"{peptide_a}//{peptide_b}"


def to_proforma(
    data: Dict[str, Any] | List[Dict[str, Any]],
    crosslinker: Optional[str | float]
) -> str | List[str]:
    r"""
    """
    _ok = check_input_multi(crosslinker, "crosslinker", [str, float]) if crosslinker is not None else True
    if isinstance(data, list):
        _ok = check_input(data, "data", list, dict)
        proforma = list()
        for item in data:
            if "data_type" not in item or item["data_type"] not in ["crosslink", "crosslink-spectrum-match"]:
                raise TypeError()
            if item["data_type"] == "crosslink":
                proforma.append(__to_proforma_xl(item, crosslinker))
            else:
                proforma.append(__to_proforma_csm(item, crosslinker))
        return proforma
    _ok = check_input(data, "data", dict)
    if "data_type" not in data or data["data_type"] not in ["crosslink", "crosslink-spectrum-match"]:
        raise TypeError()
    if data["data_type"] == "crosslink":
        return __to_proforma_xl(data, crosslinker)
    return __to_proforma_csm(data, crosslinker)
