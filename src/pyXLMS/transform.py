#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pandas as pd
from .data import check_input

from typing import List
from typing import Dict
from typing import Any

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
        alpha_proteins.append(crosslink["alpha_proteins"])
        alpha_proteins_crosslink_positions.append(crosslink["alpha_proteins_crosslink_positions"])
        alpha_decoy.append(crosslink["alpha_decoy"])
        beta_peptide.append(crosslink["beta_peptide"])
        beta_peptide_crosslink_position.append(crosslink["beta_peptide_crosslink_position"])
        beta_proteins.append(crosslink["beta_proteins"])
        beta_proteins_crosslink_positions.append(crosslink["beta_proteins_crosslink_positions"])
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

    return

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
