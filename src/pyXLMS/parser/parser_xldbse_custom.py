#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from tqdm import tqdm
from os.path import splitext

from ..data import check_input
from ..data import create_crosslink
from ..data import create_csm
from ..data import create_parser_result
from .util import format_sequence
from .util import get_bool_from_value

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List
from typing import Callable

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def pyxlms_modification_str_parser(modifications: str) -> Dict[int, Tuple[str, float]]:
    r"""Parse a pyXLMS modification string.

    Parses a pyXLMS modification string and returns the pyXLMS specific modification object,
    a dictionary that maps positions to their modififications.

    Parameters
    ----------
    modifications : str
        The pyXLMS modification string.

    Returns
    -------
    dict of int, tuple
        The pyXLMS specific modification object, a dictionary that maps positions (1-based)
        to their respective modifications given as tuples of modification name and modification
        delta mass.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed.

    Examples
    --------
    >>> from pyXLMS.parser import pyxlms_modification_str_parser
    >>> modification_str = "(1:[DSS|138.06808])"
    >>> pyxlms_modification_str_parser(modification_str)
    {1: ('DSS', 138.06808)}

    >>> from pyXLMS.parser import pyxlms_modification_str_parser
    >>> modification_str = "(1:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    >>> pyxlms_modification_str_parser(modification_str)
    {1: ('DSS', 138.06808), 7: ('Oxidation', 15.994915)}
    """
    parsed_modifications = dict()
    for mod in modifications.split(";"):
        pos = int(mod.split("(")[1].split(":")[0])
        desc = mod.split("[")[1].split("|")[0].strip()
        mass = float(mod.split("|")[1].split("]")[0])
        # if this is really in pyXLMS format we don't need to check
        # if pos already exists, because that is impossible
        # but if the parser is used for other formats that recreate the
        # same modification representation we should maybe check?
        if pos in parsed_modifications:
            raise RuntimeError(f"Modification at position {pos} already exists!")
        parsed_modifications[pos] = (desc, mass)
    return parsed_modifications


def __get_value(row: pd.Series, column: str) -> Any | None:
    r"""Get value from column if it exists and is not None.

    Parameters
    ----------
    row : pd.Series
        A row from a pandas DataFrame.
    column : str
        The column name to be accessed.

    Returns
    -------
    any, or None
        The column value if it exists and is not None.

    Notes
    -----
    This function should not be called directly, it is called from ``read_custom()``.
    """
    if column not in row:
        return None
    if (
        pd.isna(row[column])
        or row[column] is None
        or str(row[column]).lower().strip() in ["", "nan", "null", "none"]
    ):  # pyright: ignore [reportGeneralTypeIssues]
        return None
    return row[column]


def read_custom(
    files: str | List[str] | BinaryIO,
    column_mapping: Optional[Dict[str, str]] = None,
    parse_modifications: bool = True,
    modification_parser: Optional[Callable[[str], Dict[int, Tuple[str, float]]]] = None,
    decoy_prefix: str = "REV_",
    format: Literal["auto", "csv", "txt", "tsv", "xlsx"] = "auto",
    sep: str = ",",
    decimal: str = ".",
) -> Dict[str, Any]:
    r"""Read a custom or pyXLMS result file.

    Reads a custom or pyXLMS crosslink-spectrum-matches result file or crosslink result file in ``.csv`` or ``.xlsx`` format,
    and returns a ``parser_result``.

    The minimum required columns for a crosslink-spectrum-matches result file are:

    - "Alpha Peptide": The unmodified amino acid sequence of the first peptide.
    - "Alpha Peptide Crosslink Position": The position of the crosslinker in the sequence of the first peptide (1-based).
    - "Beta Peptide": The unmodified amino acid sequence of the second peptide.
    - "Beta Peptide Crosslink Position": The position of the crosslinker in the sequence of the second peptide (1-based).
    - "Spectrum File": Name of the spectrum file the crosslink-spectrum-match was identified in.
    - "Scan Nr": The corresponding scan number of the crosslink-spectrum-match.

    The minimum required columns for crosslink result file are:

    - "Alpha Peptide": The unmodified amino acid sequence of the first peptide.
    - "Alpha Peptide Crosslink Position": The position of the crosslinker in the sequence of the first peptide (1-based).
    - "Beta Peptide": The unmodified amino acid sequence of the second peptide.
    - "Beta Peptide Crosslink Position": The position of the crosslinker in the sequence of the second peptide (1-based).

    A full specification of columns that can be parsed can be found in the
    `docs <https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md>`_.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the result file(s) or a file-like object/stream.
    column_mapping : dict of str, str
        A dictionary that maps the result file columns to the required pyXLMS column names.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modification_parser' parameter.
    modification_parser : callable, or None
        A function that parses modification strings and returns the pyXLMS specific modifications object.
        If None, the function ``pyxlms_modification_str_parser()`` is used. If no modification columns are
        given this parameter is ignored.
    decoy_prefix : str, default = "REV\_"
        The prefix that indicates that a protein is from the decoy database.
    format : "auto", "csv", "tsv", "txt", or "xlsx", default = "auto"
        The format of the result file. ``"auto"`` is only available if the name/path to the result file is given.
    sep : str, default = ","
        Seperator used in the ``.csv`` or ``.tsv`` file. Parameter is ignored if the file is in ``.xlsx`` format.
    decimal : str, default = "."
        Character to recognize as decimal point. Parameter is ignored if the file is in ``.xlsx`` format.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    ValueError
        If the input format is not supported or cannot be inferred.
    TypeError
        If one of the values could not be parsed.
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslinks or crosslink-spectrum-matches.

    Examples
    --------
    >>> from pyXLMS.parser import read_custom
    >>> csms_from_pyxlms = read_custom("data/pyxlms/csm.txt")

    >>> from pyXLMS.parser import read_custom
    >>> crosslinks_from_pyxlms = read_custom("data/pyxlms/xl.txt")
    """
    ## check input
    _ok = (
        check_input(column_mapping, "column_mapping", dict, str)
        if column_mapping is not None
        else True
    )
    _ok = check_input(parse_modifications, "parse_modifications", bool)
    _ok = (
        check_input(modification_parser, "modification_parser", Callable)
        if modification_parser is not None
        else True
    )
    _ok = check_input(decoy_prefix, "decoy_prefix", str)
    _ok = check_input(format, "format", str)
    _ok = check_input(sep, "sep", str)
    _ok = check_input(decimal, "decimal", str)
    ## helper functions

    def get_is_decoy_value(
        row: pd.Series, decoy_prefix: str, alpha: bool
    ) -> bool | None:
        if alpha:
            if __get_value(row, "Alpha Decoy") is not None:
                return get_bool_from_value(__get_value(row, "Alpha Decoy"))
            if __get_value(row, "Alpha Proteins") is not None:
                return decoy_prefix in str(__get_value(row, "Alpha Proteins"))
            return None
        if __get_value(row, "Beta Decoy") is not None:
            return get_bool_from_value(__get_value(row, "Beta Decoy"))
        if __get_value(row, "Beta Proteins") is not None:
            return decoy_prefix in str(__get_value(row, "Beta Proteins"))
        return None

    def get_int(value: Any) -> int | None:
        if value is None:
            return None
        try:
            return int(value)
        except Exception as _e:
            pass
        raise TypeError(f"Could not parse int from value {value}!")
        return None

    def get_float(value: Any) -> float | None:
        if value is None:
            return None
        try:
            return float(value)
        except Exception as _e:
            pass
        raise TypeError(f"Could not parse float from value {value}!")
        return None

    ## set default parser
    if modification_parser is None:
        modification_parser = pyxlms_modification_str_parser

    ## data structures
    crosslinks = list()
    csms = list()

    ## handle input
    if not isinstance(files, list):
        inputs = [files]
    else:
        inputs = files

    for input in inputs:
        ## reading data
        data = None
        if format == "auto" and not isinstance(input, str):
            raise ValueError(
                "Can't detect format for file-like objects. Please specify format manually!"
            )
        # and isinstance specified for type checking
        if format == "auto" and isinstance(input, str):
            file_extension = splitext(input)[1].lower()
            if (
                file_extension == ".txt"
                or file_extension == ".tsv"
                or file_extension == ".csv"
            ):
                data = pd.read_csv(input, sep=sep, decimal=decimal, low_memory=False)
            elif file_extension == ".xlsx":
                data = pd.read_excel(input, engine="openpyxl")
            else:
                raise ValueError(
                    f"Detected file extension {file_extension} is not supported! Input file has to be a valid file with extension '.csv', '.tsv' or '.xlsx'!"
                )
        elif format in ["csv", "tsv", "txt", "xlsx"]:
            if format == "xlsx":
                data = pd.read_excel(input, engine="openpyxl")
            else:
                data = pd.read_csv(input, sep=sep, decimal=decimal, low_memory=False)
        else:
            raise ValueError(
                f"Provided input format {format} is not supported! Input format has to be of type 'csv', 'tsv' or 'xlsx'!"
            )
        if data is None:
            raise RuntimeError(
                "Something went wrong while reading the file! Please file a bug report!"
            )
        # this should be impossible, but check here for pyright
        if not isinstance(data, pd.DataFrame):
            raise RuntimeError(
                "Something went wrong while reading the file! Please file a bug report!"
            )
        ## detect input file type
        if column_mapping is not None:
            data.rename(columns=column_mapping, inplace=True)
        col_names = data.columns.values.tolist()
        is_crosslink_dataframe = "Scan Nr" not in col_names
        ## process data
        if is_crosslink_dataframe:
            for i, row in tqdm(
                data.iterrows(),
                total=data.shape[0],
                desc="Reading crosslinks...",
            ):
                # create crosslink
                crosslink = create_crosslink(
                    peptide_a=format_sequence(str(row["Alpha Peptide"])),
                    xl_position_peptide_a=int(row["Alpha Peptide Crosslink Position"]),
                    proteins_a=[
                        protein.strip()
                        if protein.strip()[: len(decoy_prefix)] != decoy_prefix
                        else protein.strip()[len(decoy_prefix) :]
                        for protein in str(__get_value(row, "Alpha Proteins")).split(
                            ";"
                        )
                    ]
                    if __get_value(row, "Alpha Proteins") is not None
                    else None,
                    xl_position_proteins_a=[
                        int(position)
                        for position in str(
                            __get_value(row, "Alpha Proteins Crosslink Positions")
                        ).split(";")
                    ]
                    if __get_value(row, "Alpha Proteins Crosslink Positions")
                    is not None
                    else None,
                    decoy_a=get_is_decoy_value(row, decoy_prefix, True),
                    peptide_b=format_sequence(str(row["Beta Peptide"])),
                    xl_position_peptide_b=int(row["Beta Peptide Crosslink Position"]),
                    proteins_b=[
                        protein.strip()
                        if protein.strip()[: len(decoy_prefix)] != decoy_prefix
                        else protein.strip()[len(decoy_prefix) :]
                        for protein in str(__get_value(row, "Beta Proteins")).split(";")
                    ]
                    if __get_value(row, "Beta Proteins") is not None
                    else None,
                    xl_position_proteins_b=[
                        int(position)
                        for position in str(
                            __get_value(row, "Beta Proteins Crosslink Positions")
                        ).split(";")
                    ]
                    if __get_value(row, "Beta Proteins Crosslink Positions") is not None
                    else None,
                    decoy_b=get_is_decoy_value(row, decoy_prefix, False),
                    score=get_float(__get_value(row, "Crosslink Score")),
                )
                crosslinks.append(crosslink)
        else:
            for i, row in tqdm(
                data.iterrows(),
                total=data.shape[0],
                desc="Reading CSMs...",
            ):
                # create csm
                csm = create_csm(
                    peptide_a=format_sequence(str(row["Alpha Peptide"])),
                    modifications_a=modification_parser(
                        str(__get_value(row, "Alpha Peptide Modifications"))
                    )
                    if parse_modifications
                    and __get_value(row, "Alpha Peptide Modifications") is not None
                    else None,
                    xl_position_peptide_a=int(row["Alpha Peptide Crosslink Position"]),
                    proteins_a=[
                        protein.strip()
                        if protein.strip()[: len(decoy_prefix)] != decoy_prefix
                        else protein.strip()[len(decoy_prefix) :]
                        for protein in str(__get_value(row, "Alpha Proteins")).split(
                            ";"
                        )
                    ]
                    if __get_value(row, "Alpha Proteins") is not None
                    else None,
                    xl_position_proteins_a=[
                        int(position)
                        for position in str(
                            __get_value(row, "Alpha Proteins Crosslink Positions")
                        ).split(";")
                    ]
                    if __get_value(row, "Alpha Proteins Crosslink Positions")
                    is not None
                    else None,
                    pep_position_proteins_a=[
                        int(position)
                        for position in str(
                            __get_value(row, "Alpha Proteins Peptide Positions")
                        ).split(";")
                    ]
                    if __get_value(row, "Alpha Proteins Peptide Positions") is not None
                    else None,
                    score_a=get_float(__get_value(row, "Alpha Score")),
                    decoy_a=get_is_decoy_value(row, decoy_prefix, True),
                    peptide_b=format_sequence(str(row["Beta Peptide"])),
                    modifications_b=modification_parser(
                        str(__get_value(row, "Beta Peptide Modifications"))
                    )
                    if parse_modifications
                    and __get_value(row, "Beta Peptide Modifications") is not None
                    else None,
                    xl_position_peptide_b=int(row["Beta Peptide Crosslink Position"]),
                    proteins_b=[
                        protein.strip()
                        if protein.strip()[: len(decoy_prefix)] != decoy_prefix
                        else protein.strip()[len(decoy_prefix) :]
                        for protein in str(__get_value(row, "Beta Proteins")).split(";")
                    ]
                    if __get_value(row, "Beta Proteins") is not None
                    else None,
                    xl_position_proteins_b=[
                        int(position)
                        for position in str(
                            __get_value(row, "Beta Proteins Crosslink Positions")
                        ).split(";")
                    ]
                    if __get_value(row, "Beta Proteins Crosslink Positions") is not None
                    else None,
                    pep_position_proteins_b=[
                        int(position)
                        for position in str(
                            __get_value(row, "Beta Proteins Peptide Positions")
                        ).split(";")
                    ]
                    if __get_value(row, "Beta Proteins Peptide Positions") is not None
                    else None,
                    score_b=get_float(__get_value(row, "Beta Score")),
                    decoy_b=get_is_decoy_value(row, decoy_prefix, False),
                    score=get_float(__get_value(row, "CSM Score")),
                    spectrum_file=str(row["Spectrum File"]).strip(),
                    scan_nr=int(row["Scan Nr"]),
                    charge=get_int(__get_value(row, "Precursor Charge")),
                    rt=get_float(__get_value(row, "Retention Time")),
                    im_cv=get_float(__get_value(row, "Ion Mobility")),
                )
                csms.append(csm)
    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="Custom",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )
