#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from typing import Dict
from typing import Any

def create_crosslink(peptide_a: str, xl_position_a: int,
                     peptide_b: str, xl_position_b: int) -> Dict[str, Any]:
    """Returns a crosslink dictionary.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    xl_position_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    xl_position_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).

    Returns
    -------
    crosslink : dict
        The dictionary representing the crosslink with keys alpha_peptide, alpha_crosslink_position, beta_peptide, beta_crosslink_position.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.
    """
    crosslink = sorted([f"{peptide_a};{xl_position_a}", f"{peptide_b};{xl_position_b}"])
    return {"alpha_peptide": crosslink[0].split(";")[0].strip(),
            "alpha_crosslink_position": int(crosslink[0].split(";")[1].strip()),
            "beta_peptide": crosslink[1].split(";")[0].strip(),
            "beta_crosslink_position": int(crosslink[1].split(";")[1].strip())}
