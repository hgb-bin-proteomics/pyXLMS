#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._util import check_input
from ..data._util import check_input_multi
from ._util import check_available_keys
from ._util import assert_csms
from ._util import assert_csms_or_xls

from typing import Dict
from typing import List
from typing import Set


def filter_target_decoy(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
) -> Dict[str, List[CrosslinkSpectrumMatch]] | Dict[str, List[Crosslink]]:
    r"""Seperate crosslinks or crosslink-spectrum-matches based on target and decoy matches.

    Seperates crosslinks or crosslink-spectrum-matches based on if both peptides match to the
    target database, or if both match to the decoy database, or if one of them matches to the
    target database and the other to the decoy database. The first we denote as "Target-Target"
    or "TT" matches, the second as "Decoy-Decoy" or "DD" matches, and the third as "Target-Decoy"
    or "TD" matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of pyXLMS crosslink-spectrum-matches or crosslinks.

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
    _ok = check_input(data, "data", list)
    tt = list()
    td = list()
    dd = list()
    data = assert_csms_or_xls(data)
    for item in data:
        if item["alpha_decoy"] is not None and item["beta_decoy"] is not None:
            if item["alpha_decoy"] and item["beta_decoy"]:
                dd.append(item)
            elif not item["alpha_decoy"] and not item["beta_decoy"]:
                tt.append(item)
            else:
                td.append(item)
    return {"Target-Target": tt, "Target-Decoy": td, "Decoy-Decoy": dd}


def filter_proteins(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink], proteins: Set[str] | List[str]
) -> (
    Dict[str, List[str] | List[CrosslinkSpectrumMatch]]
    | Dict[str, List[str] | List[Crosslink]]
):
    r"""Get all crosslinks or crosslink-spectrum-matches originating from proteins of interest.

    Gets all crosslinks or crosslink-spectrum-matches originating from a list of proteins of interest and
    returns a list of crosslinks or crosslink-spectrum-matches where both peptides come from a protein of
    interest and a list of crosslinks or crosslink-spectrum-matches where one of the peptides comes from a
    protein of interest.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of pyXLMS crosslink-spectrum-matches or crosslinks.
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
    _ok = check_input(data, "data", list)
    _ok = check_input_multi(proteins, "proteins", [set, list], str)
    data = assert_csms_or_xls(data)
    proteins = set(proteins)
    intra = list()
    inter = list()
    for item in data:
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
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
) -> Dict[str, List[CrosslinkSpectrumMatch]] | Dict[str, List[Crosslink]]:
    r"""Get all crosslinks or crosslink-spectrum-matches sorted by their associated proteins.

    Sorts all crosslinks or crosslink-spectrum-matches into a dictionary that maps protein
    accessions to all crosslinks or crosslink-spectrum-matches that are associated with that
    protein.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of pyXLMS crosslink-spectrum-matches or crosslinks.

    Returns
    -------
    dict
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
    _ok = check_input(data, "data", list)
    data = assert_csms_or_xls(data)
    proteins = dict()
    for item in data:
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


def filter_crosslink_type(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
) -> Dict[str, List[CrosslinkSpectrumMatch]] | Dict[str, List[Crosslink]]:
    r"""Separate crosslinks and crosslink-spectrum-matches by their crosslink type.

    Gets all crosslinks or crosslink-spectrum-matches depending on crosslink type. Will separate based
    on if a crosslink or crosslink-spectrum-match is of type "intra" or "inter" crosslink.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of pyXLMS crosslink-spectrum-matches or crosslinks.

    Returns
    -------
    dict
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
    _ok = check_input(data, "data", list)
    data = assert_csms_or_xls(data)
    intra = list()
    inter = list()
    for item in data:
        if item["crosslink_type"] == "intra":
            intra.append(item)
        else:
            inter.append(item)
    return {"Intra": intra, "Inter": inter}


def filter_peptide_pair_distribution(
    data: List[CrosslinkSpectrumMatch],
    prefix_decoys: bool = True,
) -> Dict[str, List[CrosslinkSpectrumMatch]]:
    r"""Get all crosslink-spectrum-matches sorted by their peptide pair.

    Sorts all crosslink-spectrum-matches into a dictionary that maps peptide pairs denoted as their
    amino acid sequences plus their crosslink positions delimited by a hyphen (e.g. "MTNFDKNLPNEK:6-SKLVSDFR:2")
    to their associated crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch
        A list of pyXLMS crosslink-spectrum-matches.
    prefix_decoys : bool, default = True
        Whether decoy peptides should be prefixed with a "DECOY\_" string.

    Returns
    -------
    dict of str, list of CrosslinkSpectrumMatch
        Returns a dictionary that maps peptide pairs denoted as their amino acid sequences plus their
        crosslink positions delimited by a hyphen to their associated crosslink-spectrum-matches.

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
    ['GQKNSR:3-GQKNSR:3', 'GQKNSR:3-DECOY_GSQKDR:4', 'SDKNR:3-SDKNR:3', 'DKQSGK:2-DKQSGK:2', 'DKQSGK:2-HSIKK:4']
    >>> len(
    ...     peptide_pairs["MTNFDKNLPNEK:6-SKLVSDFR:2"]
    ... )  # number of CSMs for peptide pair MTNFDKNLPNEK:6-SKLVSDFR:2
    21
    """
    _ok = check_input(data, "data", list)
    data = assert_csms(data)
    peptide_pairs = dict()
    for item in data:
        peptide_pair = (
            f"{'DECOY_' if prefix_decoys and item['alpha_decoy'] else ''}{item['alpha_peptide']}:{item['alpha_peptide_crosslink_position']}-"
            f"{'DECOY_' if prefix_decoys and item['beta_decoy'] else ''}{item['beta_peptide']}:{item['beta_peptide_crosslink_position']}"
        )
        if peptide_pair in peptide_pairs:
            peptide_pairs[peptide_pair].append(item)
        else:
            peptide_pairs[peptide_pair] = [item]
    return peptide_pairs


def filter_residue_pair_distribution(
    data: List[CrosslinkSpectrumMatch],
    prefix_decoys: bool = True,
) -> Dict[str, List[CrosslinkSpectrumMatch]]:
    r"""Get all crosslink-spectrum-matches sorted by their protein residue pair.

    Sorts all crosslink-spectrum-matches into a dictionary that maps protein residue pairs denoted as their
    protein accessions plus their protein crosslink positions delimited by a hyphen (e.g. "Cas9:48-Cas9:677")
    to their associated crosslink-spectrum-matches. If a peptide matches to more than one protein, the residues
    are delimited by commas (e.g. "Cas9:48,ALBU:36-Cas9:677").
    Requires that ``alpha_proteins``, ``beta_proteins``, ``alpha_proteins_crosslink_positions``, and
    ``beta_proteins_crosslink_positions`` fields are set for all crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch
        A list of pyXLMS crosslink-spectrum-matches.
    prefix_decoys : bool, default = True
        Whether decoy residues/proteins should be prefixed with a "DECOY\_" string.

    Returns
    -------
    dict of str, list of CrosslinkSpectrumMatch
        Returns a dictionary that maps protein residue pairs denoted as their protein accessions plus their protein
        crosslink positions delimited by a hyphen to their associated crosslink-spectrum-matches. If a peptide matches
        to more than one protein, the residues are delimited by commas.

    Raises
    ------
    TypeError
        If an unsupported data type is provided.
    RuntimeError
        If any of the crosslink-spectrum-matches do not have associated proteins or protein crosslink positions.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import filter_residue_pair_distribution
    >>> result = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> residue_pairs = filter_residue_pair_distribution(
    ...     result["crosslink-spectrum-matches"]
    ... )
    >>> list(residue_pairs.keys())[:5]  # first 5 found residue pairs
    ['Cas9:779-Cas9:779', 'Cas9:779-DECOY_Cas9:696', 'Cas9:866-Cas9:866', 'Cas9:677-Cas9:677', 'Cas9:48-Cas9:677']
    >>> len(
    ...     residue_pairs["Cas9:1122-Cas9:884"]
    ... )  # number of CSMs for residue pair Cas9:1122-Cas9:884
    22
    """
    _ok = check_input(data, "data", list)
    data = assert_csms(data)
    _ok = check_available_keys(
        [
            "alpha_proteins",
            "beta_proteins",
            "alpha_proteins_crosslink_positions",
            "beta_proteins_crosslink_positions",
        ],
        data,
    )
    residue_pairs = dict()
    for item in data:
        alpha_residue = "DECOY_" if prefix_decoys and item["alpha_decoy"] else ""
        alpha_residue += ",".join(
            sorted(
                [
                    f"{item['alpha_proteins'][i]}:{item['alpha_proteins_crosslink_positions'][i]}"
                    for i in range(len(item["alpha_proteins"]))
                ]
            )
        )
        beta_residue = "DECOY_" if prefix_decoys and item["beta_decoy"] else ""
        beta_residue += ",".join(
            sorted(
                [
                    f"{item['beta_proteins'][i]}:{item['beta_proteins_crosslink_positions'][i]}"
                    for i in range(len(item["beta_proteins"]))
                ]
            )
        )
        residue_pair = "-".join(sorted([alpha_residue, beta_residue]))
        if residue_pair in residue_pairs:
            residue_pairs[residue_pair].append(item)
        else:
            residue_pairs[residue_pair] = [item]
    return residue_pairs
