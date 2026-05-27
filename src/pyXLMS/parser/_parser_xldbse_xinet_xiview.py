#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
import pandas as pd
from tqdm import tqdm

from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._csm import create_csm
from ..data._crosslink import create_crosslink
from ..data._parser_result import create_parser_result
from ._util import format_sequence
from ._util import get_bool_from_value
from ._util import __serialize_pandas_series
from ._util import __parse_int, __parse_float

from typing import BinaryIO
from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def read_xinet(
    files: str | List[str] | BinaryIO,
    sep: str = ",",
    decimal: str = ".",
    verbose: Literal[0, 1, 2] = 1,
    **kwargs,
) -> ParserResult:
    r"""Read a xiNET exported result file.

    Reads a result file that was exported from xiNET in ``.csv`` (comma delimited) format
    and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the xiNET exported result file(s) or a file-like object/stream.
    sep : str, default = ","
        Seperator used in the ``.csv`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.
    **kwargs
        Any additional parameters will be passed to ``pandas.read*``.

    Returns
    -------
    ParserResult
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches or crosslinks.
    RuntimeError
        If the number of proteins does not match the number of protein crosslink positions. Only
        raised if verbose is set to ``2`` otherwise ``None`` will be used!
    KeyError
        If one of the required columns is not found.
    TypeError
        If parameter verbose was not set correctly.

    Warnings
    --------
    Because modifications could be encoded in very different forms depending on the xiNET/xiVIEW input source,
    the parsing of modifications is not supported with this parser! For that purpose we would recommend using
    the original result file from the corresponding crosslink search engine directly!

    Examples
    --------
    >>> from pyXLMS.parser import read_xinet
    >>> csms = read_xinet("data/xiview/DDX39B_LCSDA_shared_links_open_clamped.csv")
    """
    ## check input
    _ok = check_input(sep, "sep", str)
    _ok = check_input(decimal, "decimal", str)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")

    ## helper functions
    def __get_xl_position_peptide(row: pd.Series, alpha: bool) -> int:
        if alpha:
            if "LinkPos1" in row and not pd.isna(row["LinkPos1"]):  # pyright: ignore[reportGeneralTypeIssues]
                return __parse_int(row["LinkPos1"])
            if "PepPos1" in row and not pd.isna(row["PepPos1"]):  # pyright: ignore[reportGeneralTypeIssues]
                pep_pos1: int = __parse_int(str(row["PepPos1"]).split(";")[0])
                if "SeqPos1" in row and not pd.isna(row["SeqPos1"]):  # pyright: ignore[reportGeneralTypeIssues]
                    seq_pos1: int = __parse_int(str(row["SeqPos1"]).split(";")[0])
                    return seq_pos1 - pep_pos1 + 1
                if "AbsPos1" in row and not pd.isna(row["AbsPos1"]):  # pyright: ignore[reportGeneralTypeIssues]
                    abs_pos1: int = __parse_int(str(row["AbsPos1"]).split(";")[0])
                    return abs_pos1 - pep_pos1 + 1
        else:
            if "LinkPos2" in row and not pd.isna(row["LinkPos2"]):  # pyright: ignore[reportGeneralTypeIssues]
                return __parse_int(row["LinkPos2"])
            if "PepPos2" in row and not pd.isna(row["PepPos2"]):  # pyright: ignore[reportGeneralTypeIssues]
                pep_pos2: int = __parse_int(str(row["PepPos2"]).split(";")[0])
                if "SeqPos2" in row and not pd.isna(row["SeqPos2"]):  # pyright: ignore[reportGeneralTypeIssues]
                    seq_pos2: int = __parse_int(str(row["SeqPos2"]).split(";")[0])
                    return seq_pos2 - pep_pos2 + 1
                if "AbsPos2" in row and not pd.isna(row["AbsPos2"]):  # pyright: ignore[reportGeneralTypeIssues]
                    abs_pos2: int = __parse_int(str(row["AbsPos2"]).split(";")[0])
                    return abs_pos2 - pep_pos2 + 1
        raise KeyError(
            "Could not get a suitable column for the peptide crosslink position!"
        )
        return -1

    def __get_proteins(row: pd.Series, alpha: bool) -> List[str] | None:
        if alpha:
            if "Protein1" in row and not pd.isna(row["Protein1"]):  # pyright: ignore[reportGeneralTypeIssues]
                return [p.strip() for p in str(row["Protein1"]).split(";")]
        else:
            if "Protein2" in row and not pd.isna(row["Protein2"]):  # pyright: ignore[reportGeneralTypeIssues]
                return [p.strip() for p in str(row["Protein2"]).split(";")]
        return None

    def __get_xl_position_proteins(row: pd.Series, alpha: bool) -> List[int] | None:
        if alpha:
            if "SeqPos1" in row:
                if not pd.isna(row["SeqPos1"]):  # pyright: ignore[reportGeneralTypeIssues]
                    return [__parse_int(x) for x in str(row["SeqPos1"]).split(";")]
            if "AbsPos1" in row:
                if not pd.isna(row["AbsPos1"]):  # pyright: ignore[reportGeneralTypeIssues]
                    return [__parse_int(x) for x in str(row["AbsPos1"]).split(";")]
            if "PepPos1" in row and "LinkPos1" in row:
                if not pd.isna(row["PepPos1"]) and not pd.isna(row["LinkPos1"]):  # pyright: ignore[reportGeneralTypeIssues]
                    return [
                        __parse_int(x) + __parse_int(row["LinkPos1"]) - 1
                        for x in str(row["PepPos1"]).split(";")
                    ]
        else:
            if "SeqPos2" in row:
                if not pd.isna(row["SeqPos2"]):  # pyright: ignore[reportGeneralTypeIssues]
                    return [__parse_int(x) for x in str(row["SeqPos2"]).split(";")]
            if "AbsPos2" in row:
                if not pd.isna(row["AbsPos2"]):  # pyright: ignore[reportGeneralTypeIssues]
                    return [__parse_int(x) for x in str(row["AbsPos2"]).split(";")]
            if "PepPos2" in row and "LinkPos2" in row:
                if not pd.isna(row["PepPos2"]) and not pd.isna(row["LinkPos2"]):  # pyright: ignore[reportGeneralTypeIssues]
                    return [
                        __parse_int(x) + __parse_int(row["LinkPos2"]) - 1
                        for x in str(row["PepPos2"]).split(";")
                    ]
        return None

    def __get_spectrum_file(row: pd.Series, verbose: int) -> str:
        if "PeakListFileName" in row and not pd.isna(row["PeakListFileName"]):  # pyright: ignore[reportGeneralTypeIssues]
            return str(row["PeakListFileName"]).strip()
        if "RawFileName" in row and not pd.isna(row["RawFileName"]):  # pyright: ignore[reportGeneralTypeIssues]
            return str(row["RawFileName"]).strip()
        if "run" in row and not pd.isna(row["run"]):  # pyright: ignore[reportGeneralTypeIssues]
            return str(row["run"]).strip()
        if verbose == 2:
            raise KeyError(
                "Could not get a suitable column or value for the spectrum file name!"
            )
        return ""

    def __get_scan_number(row: pd.Series, id: int, verbose: int) -> int:
        if "ScanNumber" in row and not pd.isna(row["ScanNumber"]):  # pyright: ignore[reportGeneralTypeIssues]
            return __parse_int(row["ScanNumber"])
        if "Id" in row and not pd.isna(row["Id"]):  # pyright: ignore[reportGeneralTypeIssues]
            try:
                return __parse_int(row["Id"])
            except Exception as _e:
                pass
        if verbose == 2:
            raise KeyError(
                "Could not get a suitable column or value for the scan number!"
            )
        return id

    ## data structures
    csms = list()
    crosslinks = list()

    ## handle input
    if not isinstance(files, list):
        inputs = [files]
    else:
        inputs = files

    ## process data
    for input in inputs:
        data = pd.read_csv(input, sep=sep, decimal=decimal, low_memory=False, **kwargs)  # ty: ignore[no-matching-overload]
        has_csms = "ScanNumber" in data and (
            "run" in data or "RawFileName" in data or "PeakListFileName" in data
        )
        if "PepSeq1" not in data or "PepSeq2" not in data:
            raise KeyError("Could not get a suitable column for the peptide sequence!")
        data = data.dropna(axis=0, subset=["PepSeq1", "PepSeq2"])
        id = 0
        for i, row in tqdm(
            data.iterrows(),
            total=data.shape[0],
            desc=f"Reading xiNET/xiVIEW {'CSMs' if has_csms else 'Crosslinks'}...",
        ):
            peptide_a: str = format_sequence(str(row["PepSeq1"]))
            xl_position_peptide_a: int = __get_xl_position_peptide(row, alpha=True)
            proteins_a: List[str] | None = __get_proteins(row, alpha=True)
            xl_position_proteins_a: List[int] | None = (
                __get_xl_position_proteins(row, alpha=True)
                if proteins_a is not None
                else None
            )
            if proteins_a is not None and xl_position_proteins_a is not None:
                if len(proteins_a) != len(xl_position_proteins_a):
                    if verbose == 1:
                        warnings.warn(
                            RuntimeWarning(
                                f"Could not extract all proteins and protein crosslink positions for row with index {i}\n"
                                f"Extracted proteins: {proteins_a}\nExtracted protein crosslink positions: {xl_position_proteins_a}!"
                            )
                        )
                    if verbose == 2:
                        raise RuntimeError(
                            f"Could not extract all proteins and protein crosslink positions for row with index {i}\n"
                            f"Extracted proteins: {proteins_a}\nExtracted protein crosslink positions: {xl_position_proteins_a}!"
                        )
                    proteins_a = None
                    xl_position_proteins_a = None
            decoy_a: bool | None = (
                get_bool_from_value(row["Decoy1"]) if "Decoy1" in row else None
            )
            peptide_b: str = format_sequence(str(row["PepSeq2"]))
            xl_position_peptide_b: int = __get_xl_position_peptide(row, alpha=False)
            proteins_b: List[str] | None = __get_proteins(row, alpha=False)
            xl_position_proteins_b: List[int] | None = (
                __get_xl_position_proteins(row, alpha=False)
                if proteins_b is not None
                else None
            )
            if proteins_b is not None and xl_position_proteins_b is not None:
                if len(proteins_b) != len(xl_position_proteins_b):
                    if verbose == 1:
                        warnings.warn(
                            RuntimeWarning(
                                f"Could not extract all proteins and protein crosslink positions for row with index {i}\n"
                                f"Extracted proteins: {proteins_b}\nExtracted protein crosslink positions: {xl_position_proteins_b}!"
                            )
                        )
                    if verbose == 2:
                        raise RuntimeError(
                            f"Could not extract all proteins and protein crosslink positions for row with index {i}\n"
                            f"Extracted proteins: {proteins_b}\nExtracted protein crosslink positions: {xl_position_proteins_b}!"
                        )
                    proteins_b = None
                    xl_position_proteins_b = None
            decoy_b: bool | None = (
                get_bool_from_value(row["Decoy2"]) if "Decoy2" in row else None
            )
            score: float | None = (
                __parse_float(row["Score"]) if "Score" in row else None
            )
            id += 1
            if not has_csms:
                # create crosslink
                crosslink = create_crosslink(
                    peptide_a=peptide_a,
                    xl_position_peptide_a=xl_position_peptide_a,
                    proteins_a=proteins_a,
                    xl_position_proteins_a=xl_position_proteins_a,
                    decoy_a=decoy_a,
                    peptide_b=peptide_b,
                    xl_position_peptide_b=xl_position_peptide_b,
                    proteins_b=proteins_b,
                    xl_position_proteins_b=xl_position_proteins_b,
                    decoy_b=decoy_b,
                    score=score,
                    additional_information={
                        "source": __serialize_pandas_series(row),
                    },
                )
                crosslinks.append(crosslink)
            else:
                # create csm
                csm = create_csm(
                    peptide_a=peptide_a,
                    modifications_a=None,
                    xl_position_peptide_a=xl_position_peptide_a,
                    proteins_a=proteins_a,
                    xl_position_proteins_a=xl_position_proteins_a,
                    pep_position_proteins_a=[
                        xl_position_protein_a - xl_position_peptide_a + 1
                        for xl_position_protein_a in xl_position_proteins_a
                    ]
                    if xl_position_proteins_a is not None
                    else None,
                    score_a=None,
                    decoy_a=decoy_a,
                    peptide_b=peptide_b,
                    modifications_b=None,
                    xl_position_peptide_b=xl_position_peptide_b,
                    proteins_b=proteins_b,
                    xl_position_proteins_b=xl_position_proteins_b,
                    pep_position_proteins_b=[
                        xl_position_protein_b - xl_position_peptide_b + 1
                        for xl_position_protein_b in xl_position_proteins_b
                    ]
                    if xl_position_proteins_b is not None
                    else None,
                    score_b=None,
                    decoy_b=decoy_b,
                    score=score,
                    spectrum_file=__get_spectrum_file(row, verbose),
                    scan_nr=__get_scan_number(row, id, verbose),
                    charge=__parse_int(row["Charge"]) if "Charge" in row else None,
                    rt=None,
                    im_cv=None,
                    additional_information={
                        "source": __serialize_pandas_series(row),
                    },
                )
                csms.append(csm)
    ## check results
    if len(csms) + len(crosslinks) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="xiNET/xiVIEW",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )


def read_xiview(
    files: str | List[str] | BinaryIO,
    sep: str = ",",
    decimal: str = ".",
    verbose: Literal[0, 1, 2] = 1,
    **kwargs,
) -> ParserResult:
    r"""Read a xiVIEW exported result file.

    Reads a result file that was exported from xiVIEW in ``.csv`` (comma delimited) format
    and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the xiVIEW exported result file(s) or a file-like object/stream.
    sep : str, default = ","
        Seperator used in the ``.csv`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.
    **kwargs
        Any additional parameters will be passed to ``pandas.read*``.

    Returns
    -------
    ParserResult
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches or crosslinks.
    RuntimeError
        If the number of proteins does not match the number of protein crosslink positions. Only
        raised if verbose is set to ``2`` otherwise ``None`` will be used!
    KeyError
        If one of the required columns is not found.
    TypeError
        If parameter verbose was not set correctly.

    Notes
    -----
    Internally this just calls ``parser.read_xinet()`` since both formats share columns and the parser
    tries to exhaustively match all columns it can.

    Warnings
    --------
    Because modifications could be encoded in very different forms depending on the xiNET/xiVIEW input source,
    the parsing of modifications is not supported with this parser! For that purpose we would recommend using
    the original result file from the corresponding crosslink search engine directly!

    Examples
    --------
    >>> from pyXLMS.parser import read_xiview
    >>> csms = read_xiview("data/xiview/DDX39B_LCSDA_shared_links_open_clamped.csv")
    """
    return read_xinet(
        files,
        sep,
        decimal,
        verbose,
        **kwargs,
    )
