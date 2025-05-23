#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .data import check_input
from .data import create_csm
from .data import create_crosslink
from .data import create_parser_result
from .transform_util import assert_data_type_same

from typing import Optional
from typing import BinaryIO
from typing import Callable
from typing import Dict
from typing import Tuple
from typing import List
from typing import Any


def __get_proteins_and_positions(peptide: str, protein_db: Dict[str, str]) -> Tuple[List[str], List[int]]:
    r"""Retrieve matching proteins and peptide positions for a specific peptide.

    Matches the specified peptide against the given protein database and returns all proteins
    that contain the peptides, as well as the corresponding peptide positions in those proteins.
    Uses 0-based indexing.

    Parameters
    ----------
    peptide : str
        Unmodified peptide sequence.
    protein_db : dict of str, str
        A dictionary that maps protein accessions to their sequences.

    Returns
    -------
    tuple of list of str, list of int
        List of protein accessions, and list of peptide positions.

    Raises
    ------
    RuntimeError
        If the peptide could not be matched to any protein.

    Notes
    -----
    This function should not be called directly, it is called from ``reannotate_positions()``.

    Warnings
    --------
    Contrary to most functions in pyXLMS, this function uses 0-based indexing.
    """
    proteins = list()
    positions = list()
    for id, seq in protein_db.items():
        if peptide in seq:
            proteins.append(id)
            positions.append(seq.index(peptide))
    if len(proteins) == 0:
        raise RuntimeError(f"No match found for peptide {peptide}!")
    return (proteins, positions)


def fasta_title_to_accession(title: str) -> str:
    if "|" in title:
        return title.split("|")[1].strip()
    return title.strip()


def reannotate_positions(
    data: List[Dict[str, Any]] | Dict[str, Any],
    fasta: str | BinaryIO,
    title_to_accession: Optional[Callable[[str], str]] = None
) -> List[Dict[str, Any]] | Dict[str, Any]:
    r"""
    """
    if title_to_accession is not None:
        _ok = check_input(title_to_accession, "title_to_accession", Callable)
    else:
        title_to_accession = fasta_title_to_accession
    if isinstance(data, list):
        _ok = assert_data_type_same(data)
        protein_db = dict()
        reannoted = list()
        # read fasta file
        if isinstance(fasta, str):
            with open(fasta, "r", encoding = "utf-8") as f:
                for i, item in enumerate(SimpleFastaParser(f)):
                    protein_db[title_to_accession(item[0])] = item[1]
            if len(protein_db) != i + 1:
                warnings.warn(
                    f"Possible duplicates found in fasta file! Read {i + 1} sequences but only stored {len(protein_db)}."
                )
        else:
            for i, item in enumerate(SimpleFastaParser(fasta)):
                protein_db[title_to_accession(item[0])] = item[1]
            if len(protein_db) != i + 1:
                warnings.warn(
                    f"Possible duplicates found in fasta file! Read {i + 1} sequences but only stored {len(protein_db)}."
                )
        # annotate crosslinks
        if data[0]["data_type"] == "crosslink":
            for xl in tqdm(data, total=len(data), desc="Annotating crosslinks..."):
                proteins_a, pep_position0_proteins_a = __get_proteins_and_positions(xl["alpha_peptide"], protein_db)
                proteins_b, pep_position0_proteins_b = __get_proteins_and_positions(xl["beta_peptide"], protein_db)
                reannoted.append(create_crosslink(
                    peptide_a=xl["alpha_peptide"],
                    xl_position_peptide_a=xl["alpha_peptide_crosslink_position"],
                    proteins_a=proteins_a,
                    xl_position_proteins_a=[pos + xl["alpha_peptide_crosslink_position"] for pos in pep_position0_proteins_a],
                    decoy_a=xl["alpha_decoy"],
                    peptide_b=xl["beta_peptide"],
                    xl_position_peptide_b=xl["beta_peptide_crosslink_position"],
                    proteins_b=proteins_b,
                    xl_position_proteins_b=[pos + xl["beta_peptide_crosslink_position"] for pos in pep_position0_proteins_b],
                    decoy_b=xl["beta_decoy"],
                    score=xl["score"],
                    additional_information=xl["additional_information"],
                ))
        # annotate csms
        elif data[0]["data_type"] == "crosslink-spectrum-match":
            for csm in tqdm(data, total=len(data), desc="Annotation crosslink-spectrum-matches..."):
                proteins_a, pep_position0_proteins_a = __get_proteins_and_positions(xl["alpha_peptide"], protein_db)
                proteins_b, pep_position0_proteins_b = __get_proteins_and_positions(xl["beta_peptide"], protein_db)
                reannoted.append(create_csm(
                    peptide_a=csm["alpha_peptide"]
                    modifications_a=csm["alpha_modifications"]
                    xl_position_peptide_a=csm["alpha_peptide_crosslink_position"]
                    proteins_a=proteins_a,
                    xl_position_proteins_a=[pos + csm["alpha_peptide_crosslink_position"] for pos in pep_position0_proteins_a],
                    pep_position_proteins_a=[pos + 1 for pos in pep_position0_proteins_a],
                    score_a=csm["alpha_score"],
                    decoy_a=csm["alpha_decoy"],
                    peptide_b=csm["beta_peptide"],
                    modifications_b=csm["beta_modifications"],
                    xl_position_peptide_b=csm["beta_peptide_crosslink_position"],
                    proteins_b=proteins_b,
                    xl_position_proteins_b=[pos + csm["beta_peptide_crosslink_position"] for pos in pep_position0_proteins_b],
                    pep_position_proteins_b=[pos + 1 for pos in pep_position0_proteins_b],
                    score_b=csm["beta_score"],
                    decoy_b=csm["beta_decoy"],
                    score=csm["score"],
                    spectrum_file=csm["spectrum_file"],
                    scan_nr=csm["scan_nr"],
                    charge=csm["charge"]
                    rt=csm["retention_time"],
                    im_cv=csm["ion_mobility"],
                    additional_information=csm["additional_information"],
                ))
        else:
            raise TypeError()
        return reannoted
    _ok = check_input(data, "data", dict)
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError()
    return create_parser_result(
        search_engine = data["search_engine"],
        csms = reannotate_positions(data["crosslink-spectrum-matches"], fasta) if data["crosslink-spectrum-matches"] is not None else None,
        crosslinks= reannotate_positions(data["crosslinks"], fasta) if data["crosslinks"] is not None else None
    )
