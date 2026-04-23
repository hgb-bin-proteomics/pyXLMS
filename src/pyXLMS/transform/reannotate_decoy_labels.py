#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import copy
import warnings
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

from ..data import check_input
from ..data import create_parser_result
from .util import assert_data_type_same
from .reannotate_positions import __generate_all_sequences

from typing import BinaryIO
from typing import Callable
from typing import Dict
from typing import List
from typing import Tuple
from typing import Any


def __annotate_by_mapping(
    data: List[Dict[str, Any]], by_mapping: Dict[bool | None, bool | None]
) -> List[Dict[str, Any]]:
    r"""Reannotates decoy labels based on a given label mapping.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    by_mapping : dict of bool or None, bool or None
        A dictionary that maps possible ``alpha_decoy`` and ``beta_decoy`` values to their new values.
        For example, if decoy labels that are ``None`` should be labelled as targets, provide ``{None: False}``.

    Returns
    -------
    list of dict of str, any
        A list of reannotated crosslink-spectrum-matches or crosslinks is returned.

    Notes
    -----
    This function should not be called directly, it is called from ``reannotate_decoy_labels()``.
    """
    data_type = (
        "crosslinks"
        if data[0]["data_type"] == "crosslink"
        else "crosslink-spectrum-matches"
    )
    for _i, item in tqdm(
        enumerate(data), total=len(data), desc=f"Annotating {data_type}..."
    ):
        if item["alpha_decoy"] in by_mapping:
            item["alpha_decoy"] = by_mapping[item["alpha_decoy"]]
        if item["beta_decoy"] in by_mapping:
            item["beta_decoy"] = by_mapping[item["alpha_decoy"]]
    return data


def __annotate_by_protein_prefix(
    data: List[Dict[str, Any]], by_decoy_protein_prefix: str
) -> List[Dict[str, Any]]:
    r"""Reannotates decoy labels based on a given decoy protein prefix.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    by_decoy_protein_prefix : str
        Prefix that specifies that a protein is a decoy.

    Returns
    -------
    list of dict of str, any
        A list of reannotated crosslink-spectrum-matches or crosslinks is returned.

    Warns
    -----
    RuntimeWarning
        If one of the crosslink-spectrum-matches or crosslinks does not have assigned proteins.

    Notes
    -----
    This function should not be called directly, it is called from ``reannotate_decoy_labels()``.
    """
    data_type = (
        "crosslinks"
        if data[0]["data_type"] == "crosslink"
        else "crosslink-spectrum-matches"
    )
    for i, item in tqdm(
        enumerate(data), total=len(data), desc=f"Annotating {data_type}..."
    ):
        if item["alpha_proteins"] is not None and len(item["alpha_proteins"]) > 0:
            item["alpha_decoy"] = all(
                [
                    protein.startswith(by_decoy_protein_prefix)
                    for protein in item["alpha_proteins"]
                ]
            )
        else:
            warnings.warn(
                RuntimeWarning(
                    f"Could not annotate alpha decoy label at index={i} because alpha proteins is 'None'!"
                )
            )
        if item["beta_proteins"] is not None and len(item["beta_proteins"]) > 0:
            item["beta_decoy"] = all(
                [
                    protein.startswith(by_decoy_protein_prefix)
                    for protein in item["beta_proteins"]
                ]
            )
        else:
            warnings.warn(
                RuntimeWarning(
                    f"Could not annotate beta decoy label at index={i} because beta proteins is 'None'!"
                )
            )
    return data


def __annotate_by_protein_substring(
    data: List[Dict[str, Any]], by_decoy_protein_substring: str
) -> List[Dict[str, Any]]:
    r"""Reannotates decoy labels based on a given decoy protein substring.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    by_decoy_protein_substring : str
        Substring that specifies that a protein is a decoy.

    Returns
    -------
    list of dict of str, any
        A list of reannotated crosslink-spectrum-matches or crosslinks is returned.

    Warns
    -----
    RuntimeWarning
        If one of the crosslink-spectrum-matches or crosslinks does not have assigned proteins.

    Notes
    -----
    This function should not be called directly, it is called from ``reannotate_decoy_labels()``.
    """
    data_type = (
        "crosslinks"
        if data[0]["data_type"] == "crosslink"
        else "crosslink-spectrum-matches"
    )
    for i, item in tqdm(
        enumerate(data), total=len(data), desc=f"Annotating {data_type}..."
    ):
        if item["alpha_proteins"] is not None and len(item["alpha_proteins"]) > 0:
            item["alpha_decoy"] = all(
                [
                    by_decoy_protein_substring in protein
                    for protein in item["alpha_proteins"]
                ]
            )
        else:
            warnings.warn(
                RuntimeWarning(
                    f"Could not annotate alpha decoy label at index={i} because alpha proteins is 'None'!"
                )
            )
        if item["beta_proteins"] is not None and len(item["beta_proteins"]) > 0:
            item["beta_decoy"] = all(
                [
                    by_decoy_protein_substring in protein
                    for protein in item["beta_proteins"]
                ]
            )
        else:
            warnings.warn(
                RuntimeWarning(
                    f"Could not annotate beta decoy label at index={i} because beta proteins is 'None'!"
                )
            )
    return data


def __is_peptide_in_protein_db(peptide: str, protein_db: List[str]) -> bool:
    r"""Checks if a specific peptide is in the given protein database.

    Parameters
    ----------
    peptide : str
        Unmodified peptide sequence.
    protein_db : list of str
        A list of protein sequences.

    Returns
    -------
    bool
        Whether or not the protein database contains the peptide (``True``) or not (``False``).

    Notes
    -----
    This function should not be called directly, it is called from ``__annotate_by_fasta()``.
    """
    for base_seq in protein_db:
        seqs = __generate_all_sequences(base_seq)
        for seq in seqs:
            if peptide in seq:
                return True
    return False


def __annotate_by_fasta(
    data: List[Dict[str, Any]], fasta: str | BinaryIO, is_target: bool
) -> List[Dict[str, Any]]:
    r"""Reannotates decoy labels based on a given FASTA file.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    fasta : str, or file stream
        The name/path of the fasta file containing protein sequences or a file-like object/stream.
    is_target : bool
        If the FASTA file contains target sequences (``True``) or decoy sequences (``False``).

    Returns
    -------
    list of dict of str, any
        A list of reannotated crosslink-spectrum-matches or crosslinks is returned.

    Notes
    -----
    This function should not be called directly, it is called from ``reannotate_decoy_labels()``.
    """
    protein_db = list()
    # read fasta file
    if isinstance(fasta, str):
        with open(fasta, "r", encoding="utf-8") as f:
            for item in SimpleFastaParser(f):
                protein_db.append(item[1])
    else:
        for item in SimpleFastaParser(fasta):
            protein_db.append(item[1])
    # reannote
    data_type = (
        "crosslinks"
        if data[0]["data_type"] == "crosslink"
        else "crosslink-spectrum-matches"
    )
    for item in tqdm(data, total=len(data), desc=f"Annotating {data_type}..."):
        alpha_in_db = __is_peptide_in_protein_db(item["alpha_peptide"], protein_db)
        beta_in_db = __is_peptide_in_protein_db(item["beta_peptide"], protein_db)
        item["alpha_decoy"] = not alpha_in_db if is_target else alpha_in_db
        item["beta_decoy"] = not beta_in_db if is_target else beta_in_db
    return data


def __annotate_by_function(
    data: List[Dict[str, Any]],
    by_function: Callable[[Dict[str, Any]], Tuple[bool, bool]],
) -> List[Dict[str, Any]]:
    r"""Reannotates decoy labels based on a given function.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    by_function : callable
        A function that takes one crosslink-spectrum-match or crosslink as input and returns a tuple
        of two boolean values. The first value should be the decoy label for the alpha peptide (``True``
        if it is a decoy hit, ``False`` if it is a target hit) and the second value for the beta peptide.

    Returns
    -------
    list of dict of str, any
        A list of reannotated crosslink-spectrum-matches or crosslinks is returned.

    Notes
    -----
    This function should not be called directly, it is called from ``reannotate_decoy_labels()``.
    """
    data_type = (
        "crosslinks"
        if data[0]["data_type"] == "crosslink"
        else "crosslink-spectrum-matches"
    )
    for _i, item in tqdm(
        enumerate(data), total=len(data), desc=f"Annotating {data_type}..."
    ):
        alpha_decoy, beta_decoy = by_function(item)
        item["alpha_decoy"] = alpha_decoy
        item["beta_decoy"] = beta_decoy
    return data


def reannotate_decoy_labels(
    data: List[Dict[str, Any]] | Dict[str, Any],
    by_mapping: Dict[bool | None, bool | None] | None = None,
    by_decoy_protein_prefix: str | None = None,
    by_decoy_protein_substring: str | None = None,
    by_target_fasta: str | BinaryIO | None = None,
    by_decoy_fasta: str | BinaryIO | None = None,
    by_function: Callable[[Dict[str, Any]], Tuple[bool, bool]] | None = None,
) -> List[Dict[str, Any]] | Dict[str, Any]:
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
    _ok = (
        check_input(by_mapping, "by_mapping", dict) if by_mapping is not None else True
    )
    if by_mapping is not None:
        for k, v in by_mapping.items():
            if k not in [None, True, False] or v not in [None, True, False]:
                raise TypeError(
                    "Parameter 'by_mapping' has to be a dictionary that maps bool | None -> bool | None!"
                )
    _ok = (
        check_input(by_decoy_protein_prefix, "by_protein_prefix", str)
        if by_decoy_protein_prefix is not None
        else True
    )
    _ok = (
        check_input(by_decoy_protein_substring, "by_protein_substring", str)
        if by_decoy_protein_substring is not None
        else True
    )
    _ok = (
        check_input(by_function, "by_function", Callable)
        if by_function is not None
        else True
    )
    if [
        by_mapping,
        by_decoy_protein_prefix,
        by_decoy_protein_substring,
        by_target_fasta,
        by_decoy_fasta,
        by_function,
    ].count(None) < 5:
        raise RuntimeError(
            "Please only specify one option for reannotation, e.g. 'by_mapping' or 'by_target_fasta' but not both!"
        )
    if isinstance(data, list):
        _ok = check_input(data, "data", list, dict)
        if len(data) == 0:
            return data
        if "data_type" not in data[0]:
            raise TypeError(
                "Can't annotate positions for input data. Input data has to be a list of crosslink-spectrum-matches or crosslinks "
                "or a 'parser_result'!"
            )
        _ok = assert_data_type_same(data)
        # annotate decoy labels
        if (
            data[0]["data_type"] == "crosslink"
            or data[0]["data_type"] == "crosslink-spectrum-match"
        ):
            data_copy = copy.deepcopy(data)
            if by_mapping is not None:
                print(f"Reannotating decoy labels by mapping: {by_mapping}!")
                return __annotate_by_mapping(data_copy, by_mapping)
            if by_decoy_protein_prefix is not None:
                print(
                    f"Reannotating decoy labels by decoy protein prefix: {by_decoy_protein_prefix}!"
                )
                return __annotate_by_protein_prefix(data_copy, by_decoy_protein_prefix)
            if by_decoy_protein_substring is not None:
                print(
                    f"Reannotating decoy labels by decoy protein substring: {by_decoy_protein_substring}!"
                )
                return __annotate_by_protein_substring(
                    data_copy, by_decoy_protein_substring
                )
            if by_target_fasta is not None:
                print("Reannotating decoy labels by provided target fasta file!")
                return __annotate_by_fasta(data_copy, by_target_fasta, is_target=True)
            if by_decoy_fasta is not None:
                print("Reannotating decoy labels by provided decoy fasta file!")
                return __annotate_by_fasta(data_copy, by_decoy_fasta, is_target=False)
            if by_function is not None:
                print("Reannotating decoy labels by provided function!")
                return __annotate_by_function(data_copy, by_function)
            print(
                "No decoy label reannotation parameter provided - no decoy label reannotation has been performed!"
            )
            return data
        else:
            raise TypeError(
                f"Can't annotate decoy labels for data type {data[0]['data_type']}. Valid data types are:\n"
                "'crosslink-spectrum-match', 'crosslink', and 'parser_result'."
            )
        return data
    _ok = check_input(data, "data", dict)
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError(
            "Can't annotate decoy labels for dict. Dict has to be a valid 'parser_result'!"
        )
    new_csms = (
        reannotate_decoy_labels(
            data["crosslink-spectrum-matches"],
            by_mapping=by_mapping,
            by_decoy_protein_prefix=by_decoy_protein_prefix,
            by_decoy_protein_substring=by_decoy_protein_substring,
            by_target_fasta=by_target_fasta,
            by_decoy_fasta=by_decoy_fasta,
            by_function=by_function,
        )
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        reannotate_decoy_labels(
            data["crosslinks"],
            by_mapping=by_mapping,
            by_decoy_protein_prefix=by_decoy_protein_prefix,
            by_decoy_protein_substring=by_decoy_protein_substring,
            by_target_fasta=by_target_fasta,
            by_decoy_fasta=by_decoy_fasta,
            by_function=by_function,
        )
        if data["crosslinks"] is not None
        else None
    )
    if new_csms is not None:
        if not isinstance(new_csms, list):
            raise RuntimeError(
                "Something went wrong while reannotating decoy labels.\n"
                f"Expected data type: list. Got: {type(new_csms)}."
            )
    if new_xls is not None:
        if not isinstance(new_xls, list):
            raise RuntimeError(
                "Something went wrong while reannotating decoy labels.\n"
                f"Expected data type: list. Got: {type(new_xls)}."
            )
    return create_parser_result(
        search_engine=data["search_engine"], csms=new_csms, crosslinks=new_xls
    )
