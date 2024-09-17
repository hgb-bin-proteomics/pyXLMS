#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pandas as pd
from .data import check_input

from typing import List
from typing import Dict
from typing import Any

def modifications_to_str(modifications: Dict[int, Tuple[str, float]]) -> str:
    """Returns the string representation of a modifications dictionary.

    Parameters
    ----------
    modifications : dict of str, tuple
        The modifications of a peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.

    Returns
    -------
    str
        The string representation of the modifications.
    """
    modifications_str = ""
    for modification_pos in modifications.keys():
        modifications_str += f"({modification_pos}:[{modifications[modification_pos][0]}|{modifications[modification_pos][1]}]);"
    return modifications_str.rstrip(";")

def __cc(l: List[Any], sep: str = ";") -> str:
    """Concatenates list elements to a string using the defined seperator.

    Parameters
    ----------
    l : list
        The list to concatenate.
    sep : str, default = ";"
        The seperator to use for concatentation.

    Returns
    -------
    str
        The concatenated string of the list.
    """
    s = ""
    for i in l:
        s += str(i).strip() + sep
    return s.rstrip(";")

def __crosslinks_to_dataframe(data: List[Dict[str, Any]]) -> pd.DataFrame:
    """Returns a pandas DataFrame of the given crosslinks.

    Parameters
    ----------
    data : list
        A list of crosslinks.

    Returns
    -------
    pandas.DataFrame
        The pandas DataFrame created from the list of input crosslinks.

    Raises
    ------
    TypeError
        If the list does not contain crosslinks.
    ValueError
        If the list does not contain any objects.

    Notes
    -----
    This function should not be called directly, it is called from 'to_dataframe'.
    """
    ## columns
    alpha_peptide = list()
    alpha_peptide_crosslink_position = list()
    alpha_proteins = list()
    alpha_proteins_crosslink_positions = list()
    alpha_decoy = list()
    beta_peptide = list()
    beta_peptide_crosslink_position = list()
    beta_proteins = list()
    beta_proteins_crosslink_positions = list()
    beta_decoy = list()
    score = list()
    ## assign values
    for crosslink in data:
        alpha_peptide.append(crosslink["alpha_peptide"])
        alpha_peptide_crosslink_position.append(crosslink["alpha_peptide_crosslink_position"])
        alpha_proteins.append(__cc(crosslink["alpha_proteins"]))
        alpha_proteins_crosslink_positions.append(__cc(crosslink["alpha_proteins_crosslink_positions"]))
        alpha_decoy.append(crosslink["alpha_decoy"])
        beta_peptide.append(crosslink["beta_peptide"])
        beta_peptide_crosslink_position.append(crosslink["beta_peptide_crosslink_position"])
        beta_proteins.append(__cc(crosslink["beta_proteins"]))
        beta_proteins_crosslink_positions.append(__cc(crosslink["beta_proteins_crosslink_positions"]))
        beta_decoy.append(crosslink["beta_decoy"])
        score.append(crosslink["score"])
    return pd.DataFrame({"Alpha Peptide": alpha_peptide,
                         "Alpha Peptide Crosslink Position": alpha_peptide_crosslink_position,
                         "Alpha Proteins": alpha_proteins,
                         "Alpha Proteins Crosslink Positions": alpha_proteins_crosslink_positions,
                         "Alpha Decoy": alpha_decoy,
                         "Beta Peptide": beta_peptide,
                         "Beta Peptide Crosslink Position": beta_peptide_crosslink_position,
                         "Beta Proteins": beta_proteins,
                         "Beta Proteins Crosslink Positions": beta_proteins_crosslink_positions,
                         "Beta Decoy": beta_decoy,
                         "Crosslink Score": score})

def __csms_to_dataframe(data: List[Dict[str, Any]]) -> pd.DataFrame:
    """Returns a pandas DataFrame of the given crosslink-spectrum-matches.

    Parameters
    ----------
    data : list
        A list of crosslink-spectrum-matches.

    Returns
    -------
    pandas.DataFrame
        The pandas DataFrame created from the list of input crosslink-spectrum-matches.

    Raises
    ------
    TypeError
        If the list does not contain crosslink-spectrum-matches.
    ValueError
        If the list does not contain any objects.

    Notes
    -----
    This function should not be called directly, it is called from 'to_dataframe'.
    """
    ## columns
    alpha_peptide = list()
    alpha_modifications = list()
    alpha_peptide_crosslink_position = list()
    alpha_proteins = list()
    alpha_proteins_crosslink_positions = list()
    alpha_proteins_peptide_positions = list()
    alpha_score = list()
    alpha_decoy = list()
    beta_peptide = list()
    beta_modifications = list()
    beta_peptide_crosslink_position = list()
    beta_proteins = list()
    beta_proteins_crosslink_positions = list()
    beta_proteins_peptide_positions = list()
    beta_score = list()
    beta_decoy = list()
    score = list()
    spectrum_file = list()
    scan_nr = list()
    charge = list()
    retention_time = list()
    ion_mobility = list()
    ## assign values
    for csm in data:
        alpha_peptide.append(csm["alpha_peptide"])
        alpha_modifications.append(modifications_to_str(csm["alpha_modifications"]))
        alpha_peptide_crosslink_position.append(csm["alpha_peptide_crosslink_position"])
        alpha_proteins.append(__cc(csm["alpha_proteins"]))
        alpha_proteins_crosslink_positions.append(__cc(csm["alpha_proteins_crosslink_positions"]))
        alpha_proteins_peptide_positions.append(__cc(csm["alpha_proteins_peptide_positions"]))
        alpha_score.append(csm["alpha_score"])
        alpha_decoy.append(csm["alpha_decoy"])
        beta_peptide.append(csm["beta_peptide"])
        beta_modifications.append(modifications_to_str(csm["beta_modifications"]))
        beta_peptide_crosslink_position.append(csm["beta_peptide_crosslink_position"])
        beta_proteins.append(__cc(csm["beta_proteins"]))
        beta_proteins_crosslink_positions.append(__cc(csm["beta_proteins_crosslink_positions"]))
        beta_proteins_peptide_positions.append(__cc(csm["beta_proteins_peptide_positions"]))
        beta_score.append(csm["beta_score"])
        beta_decoy.append(csm["beta_decoy"])
        score.append(csm["score"])
        spectrum_file.append(csm["spectrum_file"])
        scan_nr.append(csm["scan_nr"])
        charge.append(csm["charge"])
        retention_time.append(csm["retention_time"])
        ion_mobility.append(csm["ion_mobility"])
    return pd.DataFrame({"Alpha Peptide": alpha_peptide,
                         "Alpha Peptide Modifications": alpha_modifications,
                         "Alpha Peptide Crosslink Position": alpha_peptide_crosslink_position,
                         "Alpha Proteins": alpha_proteins,
                         "Alpha Proteins Crosslink Positions": alpha_proteins_crosslink_positions,
                         "Alpha Proteins Peptide Positions": alpha_proteins_peptide_positions,
                         "Alpha Score": alpha_score,
                         "Alpha Decoy": alpha_decoy,
                         "Beta Peptide": beta_peptide,
                         "Beta Peptide Modifications": beta_modifications,
                         "Beta Peptide Crosslink Position": beta_peptide_crosslink_position,
                         "Beta Proteins": beta_proteins,
                         "Beta Proteins Crosslink Positions": beta_proteins_crosslink_positions,
                         "Beta Proteins Peptide Positions": beta_proteins_peptide_positions,
                         "Beta Score": beta_score,
                         "Beta Decoy": beta_decoy,
                         "CSM Score": score,
                         "Spectrum File": spectrum_file,
                         "Scan Nr": scan_nr,
                         "Precursor Charge": charge,
                         "Retention Time": retention_time,
                         "Ion Mobility": ion_mobility})

def to_dataframe(data: List[Dict[str, Any]]) -> pd.DataFrame:
    """Returns a pandas DataFrame of the given crosslinks or crosslink-spectrum-matches.

    Parameters
    ----------
    data : list
        A list of crosslinks or crosslink-spectrum-matches.

    Returns
    -------
    pandas.DataFrame
        The pandas DataFrame created from the list of input crosslinks or crosslink-spectrum-matches.

    Raises
    ------
    TypeError
        If the list does not contain crosslinks or crosslink-spectrum-matches.
    ValueError
        If the list does not contain any objects.
    """
    ## input checks
    check_input(data, "data", list, dict)
    ## function calls
    if len(data) > 0:
        if "data_type" in data[0] and data[0]["data_type"] == "crosslink":
            return __crosslinks_to_dataframe(data)
        elif "data_type" in data[0] and data[0]["data_type"] == "crosslink-spectrum-match":
            return __csms_to_dataframe(data)
        else:
            raise TypeError("The given data object is not supported!")
    else:
        raise ValueError("Parameter data has to be at least of length one!")
