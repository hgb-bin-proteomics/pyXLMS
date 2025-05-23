#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
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
    proteins = list()
    positions = list()
    for id, seq in protein_db.items():
        if peptide in seq:
            proteins.append(id)
            positions.append(seq.index(peptide))
    return (proteins, positions)



def fasta_title_to_accession(str: title) -> str:
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
    if type(data) == list:
        _ok = assert_data_type_same(data)
        protein_db = dict()
        reannoted = list()
        # read fasta file
        if type(fasta) == str:
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
            for xl in data:
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
            # todo
            pass
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
