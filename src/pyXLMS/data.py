#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from typing import List
from typing import Dict
from typing import Any

def check_input(parameter: Any,
                parameter_name: str,
                supported_class: Any,
                supported_subclass: Any = None) -> bool:
    """Checks if the given parameter is of the specified type.

    Parameters
    ----------
    parameter : any
        Parameter to check class of.
    parameter_name : str
        Name of the parameter.
    supported_class : any
        Class the parameter has to be of.
    supported_subclass : any
        Class of the values in case the parameter is a list.

    Returns
    -------
    bool
        If the given input is okay.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.
    """
    if type(parameter) != supported_class:
        raise TypeError(f"{parameter_name} must be {supported_class}!")
    if type(parameter) == list and supported_subclass is not None:
        for value in parameter:
            if type(value) != supported_subclass:
                raise TypeError(f"List values of {parameter_name} must be {supported_subclass}")
    return True

def create_crosslink(peptide_a: str,
                     xl_position_peptide_a: int,
                     proteins_a: List[str],
                     xl_position_proteins_a: List[int],
                     peptide_b: str,
                     xl_position_peptide_b: int,
                     proteins_b: List[str],
                     xl_position_proteins_b: List[int],
                     score: float) -> Dict[str, Any]:
    """Returns a crosslink dictionary.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    proteins_a: list of str
        The accessions of proteins that the first peptide is associated with.
    xl_position_proteins_a: list of int
        Positions of the crosslink in the proteins of the first peptide (1-based).
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    proteins_b: list of str
        The accessions of proteins that the second peptide is associated with.
    xl_position_proteins_b: list of int
        Positions of the crosslink in the proteins of the second peptide (1-based).
    score: float
        Score of the crosslink.

    Returns
    -------
    dict
        The dictionary representing the crosslink with keys data_type, alpha_peptide, alpha_peptide_crosslink_position,
        alpha_proteins, alpha_proteins_crosslink_positions, beta_peptide, beta_peptide_crosslink_position, beta_proteins,
        beta_proteins_crosslink_positions, and score.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.
    """
    ## input checks
    check_input(peptide_a, "peptide_a", str)
    check_input(peptide_b, "peptide_b", str)
    check_input(xl_position_peptide_a, "xl_position_peptide_a", int)
    check_input(xl_position_peptide_b, "xl_position_peptide_b", int)
    check_input(proteins_a, "proteins_a", list, str)
    check_input(proteins_b, "proteins_b", list, str)
    check_input(xl_position_proteins_a, "xl_position_proteins_a", list, int)
    check_input(xl_position_proteins_b, "xl_position_proteins_b", list, int)
    check_input(score, "score", float)
    if len(proteins_a) != len(xl_position_proteins_a):
        raise ValueError("Crosslink position has to be given for every protein! Length of proteins_a and xl_position_proteins_a has to match!")
    if len(proteins_b) != len(xl_position_proteins_b):
        raise ValueError("Crosslink position has to be given for every protein! Length of proteins_b and xl_position_proteins_b has to match!")
    ## processing
    crosslink = {f"{peptide_a.strip()}{xl_position_peptide_a}":
                    {
                        "peptide": peptide_a,
                        "xl_position_peptide": xl_position_peptide_a,
                        "proteins": proteins_a,
                        "xl_position_proteins": xl_position_proteins_a
                    },
                 f"{peptide_b.strip()}{xl_position_peptide_b}":
                    {
                        "peptide": peptide_b,
                        "xl_position_peptide": xl_position_peptide_b,
                        "proteins": proteins_b,
                        "xl_position_proteins": xl_position_proteins_b
                    }
                }
    keys = sorted(list(crosslink.keys()))
    return {"data_type": "crosslink",
            "alpha_peptide": crosslink[keys[0]]["peptide"].strip(),
            "alpha_peptide_crosslink_position": crosslink[keys[0]]["xl_position_peptide"],
            "alpha_proteins": [protein.strip() for protein in crosslink[keys[0]]["proteins"]],
            "alpha_proteins_crosslink_positions": crosslink[keys[0]]["xl_position_proteins"],
            "beta_peptide": crosslink[keys[1]]["peptide"].strip(),
            "beta_peptide_crosslink_position": crosslink[keys[1]]["xl_position_peptide"],
            "beta_proteins": [protein.strip() for protein in crosslink[keys[1]]["proteins"]],
            "beta_proteins_crosslink_positions": crosslink[keys[1]]["xl_position_proteins"],
            "score": score}
