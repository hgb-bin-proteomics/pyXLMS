#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from typing import List
from typing import Dict
from typing import Tuple
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
    if type(parameter) == dict and supported_subclass is not None:
        for key in parameter:
            if type(parameter[key]) != supported_subclass:
                raise TypeError(f"Dict values of {parameter_name} must be {supported_subclass}")
    return True

def create_crosslink(peptide_a: str,
                     xl_position_peptide_a: int,
                     proteins_a: List[str],
                     xl_position_proteins_a: List[int],
                     decoy_a: bool,
                     peptide_b: str,
                     xl_position_peptide_b: int,
                     proteins_b: List[str],
                     xl_position_proteins_b: List[int],
                     decoy_b: bool,
                     score: float) -> Dict[str, Any]:
    """Returns a crosslink dictionary.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    proteins_a : list of str
        The accessions of proteins that the first peptide is associated with.
    xl_position_proteins_a : list of int
        Positions of the crosslink in the proteins of the first peptide (1-based).
    decoy_a : bool
        Whether the alpha peptide is from the decoy database or not.
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    proteins_b : list of str
        The accessions of proteins that the second peptide is associated with.
    xl_position_proteins_b : list of int
        Positions of the crosslink in the proteins of the second peptide (1-based).
    decoy_b : bool
        Whether the beta peptide is from the decoy database or not.
    score: float
        Score of the crosslink.

    Returns
    -------
    dict
        The dictionary representing the crosslink with keys data_type, alpha_peptide, alpha_peptide_crosslink_position,
        alpha_proteins, alpha_proteins_crosslink_positions, alpha_decoy, beta_peptide, beta_peptide_crosslink_position,
        beta_proteins, beta_proteins_crosslink_positions, beta_decoy, and score.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.
    ValueError
        If the length of crosslink positions is not equal to the length of proteins.
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
    check_input(decoy_a, "decoy_a", bool)
    check_input(decoy_b, "decoy_b", bool)
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
                        "xl_position_proteins": xl_position_proteins_a,
                        "decoy": decoy_a
                    },
                 f"{peptide_b.strip()}{xl_position_peptide_b}":
                    {
                        "peptide": peptide_b,
                        "xl_position_peptide": xl_position_peptide_b,
                        "proteins": proteins_b,
                        "xl_position_proteins": xl_position_proteins_b,
                        "decoy": decoy_b
                    }
                }
    keys = sorted(list(crosslink.keys()))
    return {"data_type": "crosslink",
            "alpha_peptide": crosslink[keys[0]]["peptide"].strip(),
            "alpha_peptide_crosslink_position": crosslink[keys[0]]["xl_position_peptide"],
            "alpha_proteins": [protein.strip() for protein in crosslink[keys[0]]["proteins"]],
            "alpha_proteins_crosslink_positions": crosslink[keys[0]]["xl_position_proteins"],
            "alpha_decoy": crosslink[keys[0]]["decoy"],
            "beta_peptide": crosslink[keys[1]]["peptide"].strip(),
            "beta_peptide_crosslink_position": crosslink[keys[1]]["xl_position_peptide"],
            "beta_proteins": [protein.strip() for protein in crosslink[keys[1]]["proteins"]],
            "beta_proteins_crosslink_positions": crosslink[keys[1]]["xl_position_proteins"],
            "beta_decoy": crosslink[keys[1]]["decoy"],
            "score": score}

def create_csm(peptide_a: str,
               modifications_a: Dict[int, Tuple[str, float]],
               xl_position_peptide_a: int,
               proteins_a: List[str],
               xl_position_proteins_a: List[int],
               pep_position_proteins_a: List[int],
               score_a: float,
               decoy_a: bool,
               peptide_b: str,
               modifications_b: Dict[int, Tuple[str, float]],
               xl_position_peptide_b: int,
               proteins_b: List[str],
               xl_position_proteins_b: List[int],
               pep_position_proteins_b: List[int],
               score_b: float,
               decoy_b: bool,
               score: float,
               spectrum_file: str,
               scan_nr: int,
               charge: int,
               rt: float,
               im_cv: float) -> Dict[str, Any]:
    """Returns a crosslink-spectrum-match dictionary.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    modifications_a : dict of str, tuple
        The modifications of the first peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    proteins_a : list of str
        The accessions of proteins that the first peptide is associated with.
    xl_position_proteins_a : list of int
        Positions of the crosslink in the proteins of the first peptide (1-based).
    pep_position_proteins_a : list of int
        Positions of the first peptide in the corresponding proteins (1-based).
    score_a : float
        Identification score of the first peptide.
    decoy_a : bool
        Whether the alpha peptide is from the decoy database or not.
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    modifications_b : dict of str, tuple
        The modifications of the second peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    proteins_b : list of str
        The accessions of proteins that the second peptide is associated with.
    xl_position_proteins_b : list of int
        Positions of the crosslink in the proteins of the second peptide (1-based).
    pep_position_proteins_b : list of int
        Positions of the second peptide in the corresponding proteins (1-based).
    score_b : float
        Identification score of the second peptide.
    decoy_b : bool
        Whether the beta peptide is from the decoy database or not.
    score: float
        Score of the crosslink-spectrum-match.
    spectrum_file : str
        Name of the spectrum file the crosslink-spectrum-match was identified in.
    scan_nr : int
        The corresponding scan number of the crosslink-spectrum-match.
    charge : int
        The precursor charge of the corresponding mass spectrum of the crosslink-spectrum-match.
    rt : float
        The retention time of the corresponding mass spectrum of the crosslink-spectrum-match.
    im_cv : float
        The ion mobility or compensation voltage of the corresponding mass spectrum of the crosslink-spectrum-match.

    Returns
    -------
    dict
        The dictionary representing the crosslink-spectrum-match with keys data_type, alpha_peptide, alpha_modifications,
        alpha_peptide_crosslink_position, alpha_proteins, alpha_proteins_crosslink_positions, alpha_proteins_peptide_positions,
        alpha_score, alpha_decoy, beta_peptide, beta_modifications, beta_peptide_crosslink_position, beta_proteins,
        beta_proteins_crosslink_positions, beta_proteins_peptide_positions, beta_score, beta_decoy, score, spectrum_file, scan_nr,
        retention_time, and ion_mobility.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.
    ValueError
        If the length of crosslink positions or peptide positions is not equal to the length of proteins.
    """
    ## input checks
    check_input(peptide_a, "peptide_a", str)
    check_input(peptide_b, "peptide_b", str)
    check_input(modifications_a, "modifications_a", dict, tuple)
    check_input(modifications_b, "modifications_b", dict, tuple)
    check_input(xl_position_peptide_a, "xl_position_peptide_a", int)
    check_input(xl_position_peptide_b, "xl_position_peptide_b", int)
    check_input(proteins_a, "proteins_a", list, str)
    check_input(proteins_b, "proteins_b", list, str)
    check_input(xl_position_proteins_a, "xl_position_proteins_a", list, int)
    check_input(xl_position_proteins_b, "xl_position_proteins_b", list, int)
    check_input(pep_position_proteins_a, "pep_position_proteins_a", list, int)
    check_input(pep_position_proteins_b, "pep_position_proteins_b", list, int)
    check_input(score_a, "score_a", float)
    check_input(score_b, "score_b", float)
    check_input(decoy_a, "decoy_a", bool)
    check_input(decoy_b, "decoy_b", bool)
    check_input(score, "score", float)
    check_input(spectrum_file, "spectrum_file", str)
    check_input(scan_nr, "scan_nr", int)
    check_input(charge, "charge", int)
    check_input(rt, "rt", float)
    check_input(im_cv, "im_cv", float)
    if len(proteins_a) != len(xl_position_proteins_a):
        raise ValueError("Crosslink position has to be given for every protein! Length of proteins_a and xl_position_proteins_a has to match!")
    if len(proteins_b) != len(xl_position_proteins_b):
        raise ValueError("Crosslink position has to be given for every protein! Length of proteins_b and xl_position_proteins_b has to match!")
    if len(proteins_a) != len(pep_position_proteins_a):
        raise ValueError("Peptide position has to be given for every protein! Length of proteins_a and pep_position_proteins_a has to match!")
    if len(proteins_b) != len(pep_position_proteins_b):
        raise ValueError("Peptide position has to be given for every protein! Length of proteins_b and pep_position_proteins_b has to match!")
    ## processing
    crosslink = {f"{peptide_a.strip()}{xl_position_peptide_a}":
                    {
                        "peptide": peptide_a,
                        "modifications": {key : (modifications_a[key][0].strip(), float(modifications_a[key][1])) for key in modifications_a.keys()},
                        "xl_position_peptide": xl_position_peptide_a,
                        "proteins": proteins_a,
                        "xl_position_proteins": xl_position_proteins_a,
                        "pep_position_proteins": pep_position_proteins_a,
                        "score": score_a,
                        "decoy": decoy_a
                    },
                 f"{peptide_b.strip()}{xl_position_peptide_b}":
                    {
                        "peptide": peptide_b,
                        "modifications": {key : (modifications_b[key][0].strip(), float(modifications_b[key][1])) for key in modifications_b.keys()},
                        "xl_position_peptide": xl_position_peptide_b,
                        "proteins": proteins_b,
                        "xl_position_proteins": xl_position_proteins_b,
                        "pep_position_proteins": pep_position_proteins_b,
                        "score": score_b,
                        "decoy": decoy_b
                    }
                }
    keys = sorted(list(crosslink.keys()))
    return {"data_type": "crosslink-spectrum-match",
            "alpha_peptide": crosslink[keys[0]]["peptide"].strip(),
            "alpha_modifications": crosslink[keys[0]]["modifications"],
            "alpha_peptide_crosslink_position": crosslink[keys[0]]["xl_position_peptide"],
            "alpha_proteins": [protein.strip() for protein in crosslink[keys[0]]["proteins"]],
            "alpha_proteins_crosslink_positions": crosslink[keys[0]]["xl_position_proteins"],
            "alpha_proteins_peptide_positions": crosslink[keys[0]]["pep_position_proteins"],
            "alpha_score": crosslink[keys[0]]["score"],
            "alpha_decoy": crosslink[keys[0]]["decoy"],
            "beta_peptide": crosslink[keys[1]]["peptide"].strip(),
            "beta_modifications": crosslink[keys[1]]["modifications"],
            "beta_peptide_crosslink_position": crosslink[keys[1]]["xl_position_peptide"],
            "beta_proteins": [protein.strip() for protein in crosslink[keys[1]]["proteins"]],
            "beta_proteins_crosslink_positions": crosslink[keys[1]]["xl_position_proteins"],
            "beta_proteins_peptide_positions": crosslink[keys[1]]["pep_position_proteins"],
            "beta_score": crosslink[keys[1]]["score"],
            "beta_decoy": crosslink[keys[1]]["decoy"],
            "score": score,
            "spectrum_file": spectrum_file.strip(),
            "scan_nr": scan_nr,
            "charge": charge,
            "retention_time": rt,
            "ion_mobility": im_cv}
