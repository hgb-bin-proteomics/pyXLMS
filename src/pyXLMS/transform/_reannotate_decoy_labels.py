#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import copy
import warnings
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._util import check_input_multi
from ..data._parser_result import create_parser_result
from ._util import assert_csms
from ._util import assert_xls
from ._util import assert_csms_or_xls
from ._reannotate_positions import __generate_all_sequences

from typing import Optional
from typing import BinaryIO
from typing import Callable
from typing import Dict
from typing import List
from typing import Tuple


def __annotate_by_mapping(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    by_mapping: Dict[bool | None, bool | None],
) -> List[CrosslinkSpectrumMatch] | List[Crosslink]:
    r"""Reannotates decoy labels based on a given label mapping.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    by_mapping : dict of bool or None, bool or None
        A dictionary that maps possible ``alpha_decoy`` and ``beta_decoy`` values to their new values.
        For example, if decoy labels that are ``None`` should be labelled as targets, provide ``{None: False}``.

    Returns
    -------
    list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of reannotated crosslink-spectrum-matches or crosslinks is returned.

    Notes
    -----
    This function should not be called directly, it is called from ``reannotate_decoy_labels()``.
    """
    data_type = (
        "crosslinks" if isinstance(data[0], Crosslink) else "crosslink-spectrum-matches"
    )
    reannotated = list()
    for _i, item in tqdm(
        enumerate(data), total=len(data), desc=f"Annotating {data_type}..."
    ):
        alpha_decoy = item["alpha_decoy"]
        beta_decoy = item["beta_decoy"]
        if item["alpha_decoy"] in by_mapping:
            alpha_decoy = by_mapping[item["alpha_decoy"]]
        if item["beta_decoy"] in by_mapping:
            beta_decoy = by_mapping[item["beta_decoy"]]
        reannotated.append(
            item.copy_with_update(
                update={"alpha_decoy": alpha_decoy, "beta_decoy": beta_decoy}
            )
        )
    return reannotated


def __annotate_by_protein_prefix(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink], by_decoy_protein_prefix: str
) -> List[CrosslinkSpectrumMatch] | List[Crosslink]:
    r"""Reannotates decoy labels based on a given decoy protein prefix.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    by_decoy_protein_prefix : str
        Prefix that specifies that a protein is a decoy.

    Returns
    -------
    list of CrosslinkSpectrumMatch, or list of Crosslink
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
        "crosslinks" if isinstance(data[0], Crosslink) else "crosslink-spectrum-matches"
    )
    reannotated = list()
    for i, item in tqdm(
        enumerate(data), total=len(data), desc=f"Annotating {data_type}..."
    ):
        alpha_decoy = item["alpha_decoy"]
        beta_decoy = item["beta_decoy"]
        if item["alpha_proteins"] is not None and len(item["alpha_proteins"]) > 0:
            alpha_decoy = all(
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
            beta_decoy = all(
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
        reannotated.append(
            item.copy_with_update(
                update={"alpha_decoy": alpha_decoy, "beta_decoy": beta_decoy}
            )
        )
    return reannotated


def __annotate_by_protein_substring(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    by_decoy_protein_substring: str,
) -> List[CrosslinkSpectrumMatch] | List[Crosslink]:
    r"""Reannotates decoy labels based on a given decoy protein substring.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    by_decoy_protein_substring : str
        Substring that specifies that a protein is a decoy.

    Returns
    -------
    list of CrosslinkSpectrumMatch, or list of Crosslink
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
        "crosslinks" if isinstance(data[0], Crosslink) else "crosslink-spectrum-matches"
    )
    reannotated = list()
    for i, item in tqdm(
        enumerate(data), total=len(data), desc=f"Annotating {data_type}..."
    ):
        alpha_decoy = item["alpha_decoy"]
        beta_decoy = item["beta_decoy"]
        if item["alpha_proteins"] is not None and len(item["alpha_proteins"]) > 0:
            alpha_decoy = all(
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
            beta_decoy = all(
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
        reannotated.append(
            item.copy_with_update(
                update={"alpha_decoy": alpha_decoy, "beta_decoy": beta_decoy}
            )
        )
    return reannotated


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
        Whether the protein database contains the peptide (``True``) or not (``False``).

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
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    fasta: str | BinaryIO,
    is_target: bool,
) -> List[CrosslinkSpectrumMatch] | List[Crosslink]:
    r"""Reannotates decoy labels based on a given FASTA file.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    fasta : str, or file stream
        The name/path of the FASTA file containing protein sequences or a file-like object/stream.
    is_target : bool
        If the FASTA file contains target sequences (``True``) or decoy sequences (``False``).

    Returns
    -------
    list of CrosslinkSpectrumMatch, or list of Crosslink
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
        fasta.seek(0)
        for item in SimpleFastaParser(fasta):
            protein_db.append(item[1])
    # reannote
    data_type = (
        "crosslinks" if isinstance(data[0], Crosslink) else "crosslink-spectrum-matches"
    )
    reannotated = list()
    for item in tqdm(data, total=len(data), desc=f"Annotating {data_type}..."):
        alpha_in_db = __is_peptide_in_protein_db(item["alpha_peptide"], protein_db)
        beta_in_db = __is_peptide_in_protein_db(item["beta_peptide"], protein_db)
        alpha_decoy = not alpha_in_db if is_target else alpha_in_db
        beta_decoy = not beta_in_db if is_target else beta_in_db
        reannotated.append(
            item.copy_with_update(
                update={"alpha_decoy": alpha_decoy, "beta_decoy": beta_decoy}
            )
        )
    return reannotated


def __annotate_by_function(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    by_function: Callable[[CrosslinkSpectrumMatch | Crosslink], Tuple[bool, bool]],
) -> List[CrosslinkSpectrumMatch] | List[Crosslink]:
    r"""Reannotates decoy labels based on a given function.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks to annotate.
    by_function : callable
        A function that takes one crosslink-spectrum-match or crosslink as input and returns a tuple
        of two boolean values. The first value should be the decoy label for the alpha peptide (``True``
        if it is a decoy hit, ``False`` if it is a target hit) and the second value for the beta peptide.

    Returns
    -------
    list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of reannotated crosslink-spectrum-matches or crosslinks is returned.

    Notes
    -----
    This function should not be called directly, it is called from ``reannotate_decoy_labels()``.
    """
    data_type = (
        "crosslinks" if isinstance(data[0], Crosslink) else "crosslink-spectrum-matches"
    )
    reannotated = list()
    for _i, item in tqdm(
        enumerate(data), total=len(data), desc=f"Annotating {data_type}..."
    ):
        alpha_decoy, beta_decoy = by_function(item)
        reannotated.append(
            item.copy_with_update(
                update={"alpha_decoy": alpha_decoy, "beta_decoy": beta_decoy}
            )
        )
    return reannotated


def reannotate_decoy_labels(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult,
    by_mapping: Optional[Dict[bool | None, bool | None]] = None,
    by_decoy_protein_prefix: Optional[str] = None,
    by_decoy_protein_substring: Optional[str] = None,
    by_target_fasta: Optional[str | BinaryIO] = None,
    by_decoy_fasta: Optional[str | BinaryIO] = None,
    by_function: Optional[
        Callable[[CrosslinkSpectrumMatch | Crosslink], Tuple[bool, bool]]
    ]
    | None = None,
) -> List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult:
    r"""Reannotates decoy labels based on different parameters.

    Reannotates the decoy labels based on a provided mapping, a decoy protein prefix, a decoy protein substring,
    a target FASTA file, a decoy FASTA file, or a user-defined function. Takes a list of crosslink-spectrum-matches
    or crosslinks, or a parser_result as input.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        A list of crosslink-spectrum-matches or crosslinks to annotate, or a parser_result.
    by_mapping : dict of bool or None, bool or None, or None, default = None
        A dictionary that maps possible ``alpha_decoy`` and ``beta_decoy`` values to their new values.
        For example, if decoy labels that are ``None`` should be labelled as targets, provide ``{None: False}``.
    by_decoy_protein_prefix : str, or None, default = None
        Prefix that specifies that a protein is a decoy.
    by_decoy_protein_substring : str, or None, default = None
        Substring that specifies that a protein is a decoy.
    by_target_fasta : str, or file stream, default = None
        The name/path of the FASTA file containing target protein sequences or a file-like object/stream.
    by_decoy_fasta : str, or file stream, default = None
        The name/path of the FASTA file containing decoy protein sequences or a file-like object/stream.
    by_function : callable, or None, default = None
        A function that takes one crosslink-spectrum-match or crosslink as input and returns a tuple
        of two boolean values. The first value should be the decoy label for the alpha peptide (``True``
        if it is a decoy hit, ``False`` if it is a target hit) and the second value for the beta peptide.

    Returns
    -------
    list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        If a list of crosslink-spectrum-matches or crosslinks was provided, a list of reannotated
        crosslink-spectrum-matches or crosslinks is returned. If a parser_result was provided,
        an reannotated parser_result will be returned. Returns a copy of the original data to not
        modify the original data.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If parameter 'by_mapping' is not a dictionary that maps ``bool | None`` -> ``bool | None``.
    RuntimeError
        If more than one parameter for reannotation is given.

    Examples
    --------
    >>> from pyXLMS.data import create_crosslink_min
    >>> from pyXLMS.transform import reannotate_decoy_labels
    >>> xls = [create_crosslink_min("ADANLDK", 7, "GNTDRHSIK", 9)]
    >>> xls = reannotate_decoy_labels(xls, by_mapping={None: False})
    >>> xls[0]["alpha_decoy"]
    False
    >>> xls[0]["beta_decoy"]
    False
    """
    _ok = check_input_multi(data, "data", [list, ParserResult])
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
        if len(data) == 0:
            return data
        data = assert_csms_or_xls(data)
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
    new_csms = (
        assert_csms(
            reannotate_decoy_labels(
                data["crosslink-spectrum-matches"],
                by_mapping=by_mapping,
                by_decoy_protein_prefix=by_decoy_protein_prefix,
                by_decoy_protein_substring=by_decoy_protein_substring,
                by_target_fasta=by_target_fasta,
                by_decoy_fasta=by_decoy_fasta,
                by_function=by_function,
            )
        )
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        assert_xls(
            reannotate_decoy_labels(
                data["crosslinks"],
                by_mapping=by_mapping,
                by_decoy_protein_prefix=by_decoy_protein_prefix,
                by_decoy_protein_substring=by_decoy_protein_substring,
                by_target_fasta=by_target_fasta,
                by_decoy_fasta=by_decoy_fasta,
                by_function=by_function,
            )
        )
        if data["crosslinks"] is not None
        else None
    )
    return create_parser_result(
        search_engine=data["search_engine"],
        csms=new_csms,
        crosslinks=new_xls,
    )
