#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data import check_input
from ..data import check_input_multi

from typing import Dict
from typing import List
from typing import Set
from typing import Any


def filter_target_decoy(data: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    r"""Seperate crosslinks or crosslink-spectrum-matches based on target and decoy matches.

    Seperates crosslinks or crosslink-spectrum-matches based on if both peptides match to the
    target database, or if both match to the decoy database, or if one of them matches to the
    target database and the other to the decoy database. The first we denote as "Target-Target"
    or "TT" matches, the second as "Decoy-Decoy" or "DD" matches, and the third as "Target-Decoy"
    or "TD" matches.

    Parameters
    ----------
    data : list of dict of str, any
        A list of pyXLMS crosslinks or crosslink-spectrum-matches.

    Returns
    -------
    dict
        Returns a dictionary with key ``Target-Target`` which contains all TT matches, key ``Target-Decoy``
        which contains all TD matches, and key ``Decoy-Decoy`` which contains all DD matches.

    Raises
    ------
    TypeError
        If an unsupported data type is provided.

    Notes
    -----
    Any crosslinks or crosslink-spectrum-matches with missing 'alpha_decoy' or 'beta_decoy' attributes will be
    filtered out and not returned.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_target_decoy
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> target_and_decoys = filter_target_decoy(result["crosslink-spectrum-matches"])
    >>> len(target_and_decoys["Target-Target"])
    786
    >>> len(target_and_decoys["Target-Decoy"])
    39
    >>> len(target_and_decoys["Decoy-Decoy"])
    1

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_target_decoy
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> target_and_decoys = filter_target_decoy(result["crosslinks"])
    >>> len(target_and_decoys["Target-Target"])
    265
    >>> len(target_and_decoys["Target-Decoy"])
    0
    >>> len(target_and_decoys["Decoy-Decoy"])
    35
    """
    _ok = check_input(data, "data", list, dict)
    tt = list()
    td = list()
    dd = list()
    for item in data:
        if "data_type" not in item or item["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
            )
        if item["alpha_decoy"] is not None and item["beta_decoy"] is not None:
            if item["alpha_decoy"] and item["beta_decoy"]:
                dd.append(item)
            elif not item["alpha_decoy"] and not item["beta_decoy"]:
                tt.append(item)
            else:
                td.append(item)
    return {"Target-Target": tt, "Target-Decoy": td, "Decoy-Decoy": dd}


def filter_proteins(
    data: List[Dict[str, Any]], proteins: Set[str] | List[str]
) -> Dict[str, List[Any]]:
    r"""Get all crosslinks or crosslink-spectrum-matches originating from proteins of interest.

    Gets all crosslinks or crosslink-spectrum-matches originating from a list of proteins of interest and
    returns a list of crosslinks or crosslink-spectrum-matches where both peptides come from a protein of
    interest and a list of crosslinks or crosslink-spectrum-matches where one of the peptides comes from a
    protein of interest.

    Parameters
    ----------
    data : list of dict of str, any
        A list of pyXLMS crosslinks or crosslink-spectrum-matches.
    proteins : set of str, or list of str
        A set of protein accessions of interest.

    Returns
    -------
    dict
        Returns a dictionary with key ``Proteins`` which contains the list of proteins of interest,
        key ``Both`` which contains all crosslinks or crosslink-spectrum-matches where both peptides
        are originating from a protein of interest, and key ``One`` which contains all crosslinks or
        crosslink-spectrum-matches where one of the two peptides is originating from a protein of
        interest.

    Raises
    ------
    TypeError
        If an unsupported data type is provided.

    Notes
    -----
    Any crosslinks or crosslink-spectrum-matches with missing 'alpha_proteins' or 'beta_proteins' attributes will be
    filtered out and not returned.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_proteins
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> proteins_csms = filter_proteins(result["crosslink-spectrum-matches"], ["Cas9"])
    >>> proteins_csms["Proteins"]
    ['Cas9']
    >>> len(proteins_csms["Both"])
    798
    >>> len(proteins_csms["One"])
    23

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_proteins
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> proteins_xls = filter_proteins(result["crosslinks"], ["Cas9"])
    >>> proteins_xls["Proteins"]
    ['Cas9']
    >>> len(proteins_xls["Both"])
    274
    >>> len(proteins_xls["One"])
    21
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input_multi(proteins, "proteins", [set, list], str)
    proteins = set(proteins)
    intra = list()
    inter = list()
    for item in data:
        if "data_type" not in item or item["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
            )
        if item["alpha_proteins"] is not None and item["beta_proteins"] is not None:
            a = set(item["alpha_proteins"])
            b = set(item["beta_proteins"])
            if len(proteins.intersection(a)) > 0 and len(proteins.intersection(b)) > 0:
                intra.append(item)
            elif (
                len(proteins.intersection(a)) == 0
                and len(proteins.intersection(b)) == 0
            ):
                continue
            else:
                inter.append(item)
    return {"Proteins": list(proteins), "Both": intra, "One": inter}


def filter_protein_distribution(
    data: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    r"""Get all crosslinks or crosslink-spectrum-matches sorted by their associated proteins.

    Sorts all crosslinks or crosslink-spectrum-matches into a dictionary that maps protein
    accessions to all crosslinks or crosslink-spectrum-matches that are associated with that
    protein.

    Parameters
    ----------
    data : list of dict of str, any
        A list of pyXLMS crosslinks or crosslink-spectrum-matches.

    Returns
    -------
    dict of str, list of dict of str, any
        Returns a dictionary that maps proteins accessions (keys) to a list of crosslinks or
        crosslink-spectrum-matches (values) that are associated with that protein.

    Raises
    ------
    TypeError
        If an unsupported data type is provided.

    Notes
    -----
    Any crosslinks or crosslink-spectrum-matches with missing 'alpha_proteins' or 'beta_proteins' attributes will be
    filtered out and not returned. Please also note that the total number of crosslinks or crosslink-spectrum-matches
    returned will be greater than the size of the input because they might match to more than one protein.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_protein_distribution
    >>> result = read(
    ...     "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    ... )
    >>> proteins_csms = filter_protein_distribution(
    ...     result["crosslink-spectrum-matches"]
    ... )
    >>> list(proteins_csms.keys())  # proteins found
    ['Cas9', 'sp|MYG_HUMAN|', 'sp|CAH1_HUMAN|', 'sp|RETBP_HUMAN|', 'sp|K1C15_SHEEP|']
    >>> len(proteins_csms["Cas9"])  # number of CSMs for protein Cas9
    728
    """
    _ok = check_input(data, "data", list, dict)
    proteins = dict()
    for item in data:
        if "data_type" not in item or item["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
            )
        if item["alpha_proteins"] is not None and item["beta_proteins"] is not None:
            current_proteins = set(item["alpha_proteins"]).union(
                set(item["beta_proteins"])
            )
            for protein in current_proteins:
                if protein in proteins:
                    proteins[protein].append(item)
                else:
                    proteins[protein] = [item]
    return proteins


def filter_crosslink_type(data: List[Dict[str, Any]]) -> Dict[str, List[Any]]:
    r"""Separate crosslinks and crosslink-spectrum-matches by their crosslink type.

    Gets all crosslinks or crosslink-spectrum-matches depending on crosslink type. Will separate based
    on if a crosslink or crosslink-spectrum-match is of type "intra" or "inter" crosslink.

    Parameters
    ----------
    data : list of dict of str, any
        A list of pyXLMS crosslinks or crosslink-spectrum-matches.

    Returns
    -------
    dict of str, list of dict
        Returns a dictionary with key ``Intra`` which contains all crosslinks or crosslink-spectrum-
        matches with crosslink type = "intra", and key ``Inter`` which contains all crosslinks or
        crosslink-spectrum-matches with crosslink type = "inter".

    Raises
    ------
    TypeError
        If an unsupported data type is provided.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_crosslink_type
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> crosslink_type_filtered_csms = filter_crosslink_type(
    ...     result["crosslink-spectrum-matches"]
    ... )
    >>> len(crosslink_type_filtered_csms["Intra"])
    803
    >>> len(crosslink_type_filtered_csms["Inter"])
    23

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_crosslink_type
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> crosslink_type_filtered_crosslinks = filter_crosslink_type(result["crosslinks"])
    >>> len(crosslink_type_filtered_crosslinks["Intra"])
    279
    >>> len(crosslink_type_filtered_crosslinks["Inter"])
    21
    """
    _ok = check_input(data, "data", list, dict)
    intra = list()
    inter = list()
    for item in data:
        if "data_type" not in item or item["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
            )
        if item["crosslink_type"] == "intra":
            intra.append(item)
        else:
            inter.append(item)
    return {"Intra": intra, "Inter": inter}


def filter_peptide_pair_distribution(
    data: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    r"""Get all crosslink-spectrum-matches sorted by their peptide pair.

    Sorts all crosslink-spectrum-matches into a dictionary that maps peptide pairs denoted as their
    amino acid sequences delimited by a hyphen (e.g. "MTNFDKNLPNEK-SKLVSDFR") to their associated
    crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of dict of str, any
        A list of pyXLMS crosslink-spectrum-matches.

    Returns
    -------
    dict of str, list of dict
        Returns a dictionary that maps peptide pairs denoted as their amino acid sequences delimited by a
        hyphen to their associated crosslink-spectrum-matches.

    Raises
    ------
    TypeError
        If an unsupported data type is provided.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_peptide_pair_distribution
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> peptide_pairs = filter_peptide_pair_distribution(
    ...     result["crosslink-spectrum-matches"]
    ... )
    >>> list(peptide_pairs.keys())[:5]  # first 5 found peptide pairs
    ['GQKNSR-GQKNSR', 'GQKNSR-GSQKDR', 'SDKNR-SDKNR', 'DKQSGK-DKQSGK', 'DKQSGK-HSIKK']
    >>> len(
    ...     peptide_pairs["MTNFDKNLPNEK-SKLVSDFR"]
    ... )  # number of CSMs for peptide pair MTNFDKNLPNEK-SKLVSDFR
    21
    """
    _ok = check_input(data, "data", list, dict)
    peptide_pairs = dict()
    for item in data:
        if "data_type" not in item or item["data_type"] != "crosslink-spectrum-match":
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink-spectrum-match!"
            )
        peptide_pair = f"{item['alpha_peptide']}-{item['beta_peptide']}"
        if peptide_pair in peptide_pairs:
            peptide_pairs[peptide_pair].append(item)
        else:
            peptide_pairs[peptide_pair] = [item]
    return peptide_pairs
