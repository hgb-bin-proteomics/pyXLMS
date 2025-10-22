#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
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
from ..transform.reannotate_positions import __generate_all_sequences
from ..transform.util import assert_data_type_same

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Tuple
from typing import List
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

CHAINS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")


def __get_proteins_and_positions(
    peptide: str, protein_db: Dict[str, Dict[str, str]], error_on_no_match: bool = False
) -> Tuple[List[str], List[int]]:
    r"""Retrieve matching protein chains and peptide positions for a specific peptide.

    Matches the specified peptide against the given protein database and returns all protein chains
    that contain the peptide, as well as the corresponding peptide positions in those protein chains.
    Uses 0-based indexing.

    Parameters
    ----------
    peptide : str
        Unmodified peptide sequence.
    protein_db : dict of dict of str, str
        A dictionary that maps protein chain ids to their fasta entries, which are dictionaries
        that map key "header" to the sequence header and key "sequence" to the sequence.
    error_on_no_match : bool, default = False
        Wether an error should be raised if the peptide matches to none of the proteins.

    Returns
    -------
    tuple of list of str, list of int
        List of protein chain ids, and list of peptide positions.

    Raises
    ------
    RuntimeError
        If the peptide could not be matched to any protein. Only raised if 'error_on_no_match'
        is set to True.

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
        base_seq = item["sequence"]
        seqs = __generate_all_sequences(base_seq)
        for seq in seqs:
            if peptide in seq:
                for match in re.finditer(peptide, seq):
                    proteins.append(chain)
                    positions.append(match.start())
    if error_on_no_match and len(proteins) == 0:
        raise RuntimeError(f"No match found for peptide {peptide}!")
    return (proteins, positions)


def __protein_supported_by_crosslink(
    sequence: str, crosslinks: List[Dict[str, Any]]
) -> bool:
    r"""Check if a specific protein is supported by any of the given crosslinks.

    Checks if a specific protein that is given by its sequence is supported by any of the input
    crosslinks.

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
    seqs = __generate_all_sequences(sequence)
    for seq in seqs:
        for crosslink in crosslinks:
            if crosslink["alpha_peptide"] in seq:
                return True
            if crosslink["beta_peptide"] in seq:
                return True
    return False


def to_alphalink2(
    crosslinks: List[Dict[str, Any]],
    fasta: str | BinaryIO,
    annotated_fdr: float | List[float] = 0.01,
    try_use_annotated_fdr: bool = True,
    filename_prefix: Optional[str] = None,
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[str, Any]:
    r"""Exports a list of crosslinks to AlphaLink2 format.

    Exports a list of crosslinks to AlphaLink2 format. The tool AlphaLink2 is accessible
    via the link
    `github.com/Rappsilber-Laboratory/AlphaLink2 <https://github.com/Rappsilber-Laboratory/AlphaLink2>`_.

    Parameters
    ----------
    crosslinks : list of dict of str, any
        A list of crosslinks to export.
    fasta : str, or file stream
        The name/path of the fasta file containing protein sequences or a file-like object/stream.
        Please keep in mind that AlphaLink2 supports a maximum of 62 proteins/chains!
    annotated_fdr : float, or list of float, default = 0.01
        Value(s) to use for the "FDR" column in the AlphaLink2 crosslink table. If a single float value
        is given, all crosslinks will be annotated with that constant FDR value. If a list of floats is
        given, it should be of equal length as the crosslinks and correspondingly contain FDR values for
        all crosslinks in the same order, e.g. the FDR value for the first crosslink should be the first
        value in the list as well, the second for the second crosslink, and so forth.
    try_use_annotated_fdr : bool, default = True
        If pyXLMS annotated FDR - e.g. via ``transform.annotate_fdr()`` - should be used for the "FDR"
        column in the AlphaLink2 crosslink table, if available. This will override the values given via
        parameter ``annotated_fdr`` and should be ``False`` if you are passing a custom list of FDR values!
    filename_prefix : str, or None, default = None
        If not None, the exported data will be written to files with the specified filename prefix.
        The full list of written files can be accessed via the returned dictionary.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    dict of str, any
        Returns a dictionary with key ``AlphaLink2 crosslinks`` containing the formatted crosslink input for AlphaLink2,
        with key ``AlphaLink2 FASTA`` containing the FASTA file content for AlphaLink2,
        with key ``AlphaLink2 DataFrame`` containing the exported crosslinks as a pandas DataFrame,
        and with key ``Exported files`` containing a list of filenames of all files that were
        written to disk.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    ValueError
        If the provided crosslinks contain no elements.
    ValueError
        If the length of annotated_fdr does not match the length of the input crosslinks.
        Only applies if annotated_fdr is given as a list.
    IndexError
        If the provided crosslinks match to more than 62 proteins/chains in the FASTA file.
    RuntimeError
        If one or more of the crosslinks could not be matched to the proteins/chains in the
        FASTA file. Only is raised if 'verbose' is set to 2.

    Notes
    -----
    Please note that the legacy PDB format and therefore also AlphaLink2 only supports a maximum
    of 62 proteins/chains. If crosslinks of more than 62 proteins/chains are given, an error will
    be thrown!

    Examples
    --------
    >>> from pyXLMS.pipelines import pipeline
    >>> from pyXLMS.transform import filter_proteins
    >>> from pyXLMS.exporter import to_alphalink2
    >>> pr = pipeline(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> cas9 = filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]
    >>> _export = to_alphalink2(
    ...     cas9, fasta="data/_fasta/Cas9_plus10.fasta", filename_prefix="Cas9"
    ... )
    """
    _ok = check_input(crosslinks, "crosslinks", list, dict)
    _ok = check_input_multi(annotated_fdr, "annotated_fdr", [float, list], float)
    _ok = check_input(try_use_annotated_fdr, "try_use_annotated_fdr", bool)
    _ok = (
        check_input(filename_prefix, "filename_prefix", str)
        if filename_prefix is not None
        else True
    )
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")
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
        raise IndexError(
            "Found more than the supported 62 proteins/chains in the fasta file! Please trim fasta file to a maximum of 62 sequences!"
        )
    id = 0
    for item in fasta_items:
        if __protein_supported_by_crosslink(item[1], crosslinks):
            protein_db[CHAINS[id]] = {"header": item[0], "sequence": item[1]}
            id += 1
    # prepare fdr values
    fdr_values = (
        annotated_fdr
        if isinstance(annotated_fdr, list)
        else [annotated_fdr for xl in crosslinks]
    )
    # output
    alphalink2_txt = ""
    alphalink2_df_dict = {
        "residueFrom": [],
        "chain1": [],
        "residueTo": [],
        "chain2": [],
        "FDR": [],
    }
    # export crosslinks
    for id, xl in tqdm(
        enumerate(crosslinks),
        total=len(crosslinks),
        desc="Exporting crosslinks to AlphaLink2...",
    ):
        proteins_a, pep_position0_proteins_a = __get_proteins_and_positions(
            xl["alpha_peptide"], protein_db, error_on_no_match=False
        )
        proteins_b, pep_position0_proteins_b = __get_proteins_and_positions(
            xl["beta_peptide"], protein_db, error_on_no_match=False
        )
        # skip if crosslink matches to none of the proteins in the fasta
        if len(proteins_a) == 0 or len(proteins_b) == 0:
            if verbose == 1:
                warnings.warn(
                    RuntimeWarning(
                        (
                            f"Could not find matching proteins in FASTA file for crosslink id {id}:{xl['alpha_peptide']}-{xl['beta_peptide']}! "
                            "This warning can be ignored if this is to be expected."
                        )
                    )
                )
            if verbose == 2:
                raise RuntimeError(
                    (
                        f"Could not find matching proteins in FASTA file for crosslink id {id}:{xl['alpha_peptide']}-{xl['beta_peptide']}! "
                        "If this is to be expected please set verbose level to either 1 or 0!"
                    )
                )
            continue
        for i in range(len(proteins_a)):
            for j in range(len(proteins_b)):
                residueFrom = (
                    pep_position0_proteins_a[i] + xl["alpha_peptide_crosslink_position"]
                )
                chain1 = proteins_a[i]
                residueTo = (
                    pep_position0_proteins_b[j] + xl["beta_peptide_crosslink_position"]
                )
                chain2 = proteins_b[j]
                FDR = fdr_values[id]
                if try_use_annotated_fdr:
                    if xl["additional_information"] is not None:
                        if "pyXLMS_annotated_FDR" in xl["additional_information"]:
                            if not pd.isna(
                                xl["additional_information"]["pyXLMS_annotated_FDR"]
                            ):
                                FDR = xl["additional_information"][
                                    "pyXLMS_annotated_FDR"
                                ]
                alphalink2_df_dict["residueFrom"].append(residueFrom)
                alphalink2_df_dict["chain1"].append(chain1)
                alphalink2_df_dict["residueTo"].append(residueTo)
                alphalink2_df_dict["chain2"].append(chain2)
                alphalink2_df_dict["FDR"].append(FDR)
                alphalink2_txt += f"{residueFrom} {chain1} {residueTo} {chain2} {FDR}\n"
    # create fasta
    alphalink2_fasta = ""
    for chain in CHAINS:
        if chain in protein_db:
            item = protein_db[chain]
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
        "Exported files": exported_files,
    }
