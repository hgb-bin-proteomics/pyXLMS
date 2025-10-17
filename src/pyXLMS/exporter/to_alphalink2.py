#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import re
import warnings
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

from ..data import check_input
from ..data import check_input_multi
from ..data import create_csm
from ..data import create_crosslink
from ..data import create_parser_result
from ..transform.util import assert_data_type_same

from typing import Optional
from typing import BinaryIO
from typing import Callable
from typing import Dict
from typing import Tuple
from typing import List
from typing import Any

CHAINS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")

def __get_proteins_and_positions(
    peptide: str, protein_db: Dict[str, str]
) -> Tuple[List[str], List[int]]:
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
            for match in re.finditer(peptide, seq):
                proteins.append(id)
                positions.append(match.start())
    if len(proteins) == 0:
        raise RuntimeError(f"No match found for peptide {peptide}!")
    return (proteins, positions)


def to_alphalink2(
    crosslinks: List[Dict[str, Any]],
    fasta: str | BinaryIO,
    annotated_fdr: float | List[float] = 0.01
) -> Dict[str, Any]:
    r"""Reannotates protein crosslink positions for a given fasta file.

    Reannotates the crosslink and peptide positions of the given cross-linked peptide pair and
    the specified fasta file. Takes a list of crosslink-spectrum-matches or crosslinks, or a
    parser_result as input.

    Parameters
    ----------
    data : list of dict of str, any, or dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate, or a parser_result.
    fasta : str, or file stream
        The name/path of the fasta file containing protein sequences or a file-like object/stream.
    title_to_accession : callable, or None, default = None
        A function that parses the protein accession from the fasta title/header. If None (default)
        the function ``fasta_title_to_accession`` is used.

    Returns
    -------
    list of dict of str, any, or dict of str, any
        If a list of crosslink-spectrum-matches or crosslinks was provided, a list of annotated
        crosslink-spectrum-matches or crosslinks is returned. If a parser_result was provided,
        an annotated parser_result will be returned.

    Raises
    ------
    TypeError
        If a wrong data type is provided.

    Examples
    --------
    >>> from pyXLMS.data import create_crosslink_min
    >>> from pyXLMS.transform import reannotate_positions
    >>> xls = [create_crosslink_min("ADANLDK", 7, "GNTDRHSIK", 9)]
    >>> xls = reannotate_positions(xls, "data/_fasta/Cas9_plus10.fasta")
    >>> xls[0]["alpha_proteins"]
    ["Cas9"]
    >>> xls[0]["alpha_proteins_crosslink_positions"]
    [1293]
    >>> xls[0]["beta_proteins"]
    ["Cas9"]
    >>> xls[0]["beta_proteins_crosslink_positions"]
    [48]
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input_multi(annotated_fdr, "annotated_fdr", [float, list], float)
    if len(data) == 0:
        raise ValueError()
    if "data_type" not in data[0] or data[0]["data_type"] != "crosslink":
        raise TypeError(
            "Can't annotate positions for input data. Input data has to be a list of crosslink-spectrum-matches or crosslinks "
            "or a 'parser_result'!"
        )
    _ok = assert_data_type_same(data)
    protein_db = dict()
    reannoted = list()
    # read fasta file
    i = 0
    if isinstance(fasta, str):
        with open(fasta, "r", encoding="utf-8") as f:
            for i, item in enumerate(SimpleFastaParser(f)):
                protein_db[title_to_accession(item[0])] = item[1]
        if len(protein_db) != i + 1:
            warnings.warn(
                RuntimeWarning(
                    f"Possible duplicates found in fasta file! Read {i + 1} sequences but only stored {len(protein_db)}."
                )
            )
    else:
        for i, item in enumerate(SimpleFastaParser(fasta)):
            protein_db[title_to_accession(item[0])] = item[1]
        if len(protein_db) != i + 1:
            warnings.warn(
                RuntimeWarning(
                    f"Possible duplicates found in fasta file! Read {i + 1} sequences but only stored {len(protein_db)}."
                )
            )
    # annotate crosslinks
    for xl in tqdm(data, total=len(data), desc="Annotating crosslinks..."):
        proteins_a, pep_position0_proteins_a = __get_proteins_and_positions(
            xl["alpha_peptide"], protein_db
        )
        proteins_b, pep_position0_proteins_b = __get_proteins_and_positions(
            xl["beta_peptide"], protein_db
        )
        reannoted.append(
            create_crosslink(
                peptide_a=xl["alpha_peptide"],
                xl_position_peptide_a=xl["alpha_peptide_crosslink_position"],
                proteins_a=proteins_a,
                xl_position_proteins_a=[
                    pos + xl["alpha_peptide_crosslink_position"]
                    for pos in pep_position0_proteins_a
                ],
                decoy_a=xl["alpha_decoy"],
                peptide_b=xl["beta_peptide"],
                xl_position_peptide_b=xl["beta_peptide_crosslink_position"],
                proteins_b=proteins_b,
                xl_position_proteins_b=[
                    pos + xl["beta_peptide_crosslink_position"]
                    for pos in pep_position0_proteins_b
                ],
                decoy_b=xl["beta_decoy"],
                score=xl["score"],
                additional_information=xl["additional_information"],
            )
        )
    return create_parser_result(
        search_engine=data["search_engine"], csms=new_csms, crosslinks=new_xls
    )
