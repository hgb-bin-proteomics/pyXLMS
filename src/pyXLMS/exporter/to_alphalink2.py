#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import re
import warnings
import pandas as pd
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
    peptide: str, protein_db: Dict[str, Dict[str, str]]
) -> Tuple[List[str], List[int]]:
    r"""Retrieve matching protein chains and peptide positions for a specific peptide.

    Matches the specified peptide against the given protein database and returns all protein chains
    that contain the peptides, as well as the corresponding peptide positions in those protein chains.
    Uses 0-based indexing.

    Parameters
    ----------
    peptide : str
        Unmodified peptide sequence.
    protein_db : dict of dict of str, str
        A dictionary that maps protein chain ids to their fasta entries, which are dictionaries
        that map key "header" to the sequence header and "sequence" to the sequence.

    Returns
    -------
    tuple of list of str, list of int
        List of protein chain ids, and list of peptide positions.

    Raises
    ------
    RuntimeError
        If the peptide could not be matched to any protein.

    Notes
    -----
    This function should not be called directly, it is called from ``to_alphalink2()``.

    Warnings
    --------
    Contrary to most functions in pyXLMS, this function uses 0-based indexing.
    """
    proteins = list()
    positions = list()
    for chain, item in protein_db.items():
        seq = item["sequence"] 
        if peptide in seq:
            for match in re.finditer(peptide, seq):
                proteins.append(chain)
                positions.append(match.start())
    if len(proteins) == 0:
        raise RuntimeError(f"No match found for peptide {peptide}!")
    return (proteins, positions)


def __protein_supported_by_crosslink(sequence: str, crosslinks: List[Dict[str, Any]]) -> bool:
    r"""Retrieve matching protein chains and peptide positions for a specific peptide.

    Matches the specified peptide against the given protein database and returns all protein chains
    that contain the peptides, as well as the corresponding peptide positions in those protein chains.
    Uses 0-based indexing.

    Parameters
    ----------
    sequence : str
        The sequence of the protein.
    crosslinks : list of dict of str, any
        A list of crosslinks.

    Returns
    -------
    bool
        Returns True if the protein is supported by any crosslink, otherwise False.

    Notes
    -----
    This function should not be called directly, it is called from ``to_alphalink2()``.
    """
    for crosslink in crosslinks:
        if crosslink["alpha_peptide"] in sequence:
            return True
        if crosslink["beta_peptide"] in sequence:
            return True
    return False


def to_alphalink2(
    crosslinks: List[Dict[str, Any]],
    fasta: str | BinaryIO,
    annotated_fdr: float | List[float] = 0.01,
    try_use_annotated_fdr: bool = True,
    filename_prefix: Optional[str] = None,
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
    _ok = check_input(crosslinks, "crosslinks", list, dict)
    _ok = check_input_multi(annotated_fdr, "annotated_fdr", [float, list], float)
    _ok = check_input(try_use_annotated_fdr, "try_use_annotated_fdr", bool)
    _ok = (
        check_input(filename_prefix, "filename_prefix", str)
        if filename_prefix is not None
        else True
    )
    if isinstance(annotated_fdr, list) and len(annotated_fdr) != len(crosslinks):
        raise ValueError(
            "Length of annotated_fdr does not match length of crosslinks! "
            + "When providing a list it needs to contain FDR values for every crosslink and therefore be of equal length!"
        )
    if len(crosslinks) == 0:
        raise ValueError("Provided crosslinks contain no elements!")
    if "data_type" not in crosslinks[0] or crosslinks[0]["data_type"] != "crosslink":
        raise TypeError(
            "Unsupported data type for input crosslinks! Parameter crosslinks has to be a list of crosslinks!"
        )
    _ok = assert_data_type_same(crosslinks)
    protein_db = dict()
    # read fasta file
    fasta_items = list()
    if isinstance(fasta, str):
        with open(fasta, "r", encoding="utf-8") as f:
            for item in SimpleFastaParser(f):
                fasta_items.append(item)
    else:
        for item in SimpleFastaParser(fasta):
            fasta_items.append(item)
    if len(fasta_items) > len(CHAINS):
        raise IndexError("Found more than the supported 62 proteins/chains in the fasta file! Please trim fasta file to a maximum of 62 sequences!")
    id = 0
    for item in fasta_items:
        if __protein_supported_by_crosslink(item[1], crosslinks):
            protein_db[CHAINS[id]] = {"header": item[0], "sequence": item[1]}
            id += 1
    # prepare fdr values
    fdr_values = annotated_fdr if isinstance(annotated_fdr, list) else [annotated_fdr for xl in crosslinks]
    # output
    alphalink2_txt = ""
    alphalink2_df_dict = {"residueFrom": [], "chain1": [], "residueTo": [], "chain2": [], "FDR": []}
    # export crosslinks
    for id, xl in tqdm(enumerate(crosslinks), total=len(crosslinks), desc="Exporting crosslinks to AlphaLink2..."):
        proteins_a, pep_position0_proteins_a = __get_proteins_and_positions(
            xl["alpha_peptide"], protein_db
        )
        proteins_b, pep_position0_proteins_b = __get_proteins_and_positions(
            xl["beta_peptide"], protein_db
        )
        for i in range(len(proteins_a)):
            for j in range(len(proteins_b)):
                residueFrom = pep_position0_proteins_a[i] + xl["alpha_peptide_crosslink_position"]
                chain1 = proteins_a[i]
                residueTo = pep_position0_proteins_b[j] + xl["beta_peptide_crosslink_position"]
                chain2 = proteins_b[j]
                FDR = fdr_values[id]
                if try_use_annotated_fdr:
                    if xl["additional_information"] is not None:
                        if "pyXLMS_annotated_FDR" in xl["additional_information"]:
                            if not pd.isna(xl["additional_information"]["pyXLMS_annotated_FDR"]):
                                FDR = xl["additional_information"]["pyXLMS_annotated_FDR"]
                alphalink2_df_dict["residueFrom"].append(residueFrom)
                alphalink2_df_dict["chain1"].append(chain1)
                alphalink2_df_dict["residueTo"].append(residueTo)
                alphalink2_df_dict["chain2"].append(chain2)
                alphalink2_df_dict["FDR"].append(FDR)
                alphalink2_txt += f"{residueFrom} {chain1} {residueTo} {chain2} {FDR}\n"
    # create fasta
    alphalink2_fasta = ""
    for chain, item in protein_db.items():
        alphalink2_fasta += f">{chain}|{item['header']}\n{item['sequence']}\n"
    # create pandas dataframe
    alphalink2_df = pd.DataFrame(alphalink2_df_dict)
    # export files
    exported_files = list()
    if filename_prefix is not None:
        with open(filename_prefix + "_AlphaLink2.txt", "w", encoding="utf-8") as f:
            f.write(alphalink2_txt)
            f.close()
        exported_files.append(filename_prefix + "_AlphaLink2.txt")
        with open(filename_prefix + "_AlphaLink2.fasta", "w", encoding="utf-8") as f:
            f.write(alphalink2_fasta)
            f.close()
        exported_files.append(filename_prefix + "_AlphaLink2.fasta")
    # return exported files
    return {
        "AlphaLink2 crosslinks": alphalink2_txt,
        "AlphaLink2 FASTA": alphalink2_fasta,
        "AlphaLink2 DataFrame": alphalink2_df,
        "Exported files": exported_files
    }
