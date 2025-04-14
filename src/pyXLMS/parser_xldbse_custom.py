#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import sqlite3
import pandas as pd
from tqdm import tqdm
from os.path import splitext

from .data import check_input
from .data import create_crosslink
from .data import create_csm
from .data import create_parser_result
from .constants import MODIFICATIONS
from .parser_util import format_sequence
from .parser_util import get_bool_from_value

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
    {1: ("DSS", 138.06808)}

    >>> from pyXLMS.parser import pyxlms_modification_str_parser
    >>> modification_str = "(1:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    >>> pyxlms_modification_str_parser(modification_str)
    {1: ("DSS", 138.06808), 7: ("Oxidation", 15.994915)}
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
    ):
        return None
    return row[column]


def read_custom(
    files: str | List[str] | BinaryIO,
    column_mapping: Optional[Dict[str, str]] = None,
    modification_parser: Optional[Callable[[str], Dict[int, Tuple[str, float]]]] = None,
    format: Literal["auto", "csv", "txt", "tsv", "xlsx"] = "auto",
    sep: str = "\t",
) -> Dict[str, Any]:
    """Read a custom or pyXLMS result file.

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
    `docs <https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/md/format.md>`_.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MS Annika result file(s) or a file-like object/stream.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
    format : "auto", "csv", "tsv", "txt", "xlsx", or "pdresult", default = "auto"
        The format of the result file. ``"auto"`` is only available if the name/path to the MS Annika result file is given.
    sep : str, default = "\t"
        Seperator used in the ``.csv`` or ``.tsv`` file. Parameter is ignored if the file is in ``.xlsx`` format.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    ValueError
        If the input format is not supported or cannot be inferred.
    TypeError
        If the pdResult file is provided in the wrong format.
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslinks or crosslink-spectrum-matches.
    KeyError
        If one of the found post-translational-modifications could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_msannika
    >>> csms_from_xlsx = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")

    >>> from pyXLMS.parser import read_msannika
    >>> crosslinks_from_xlsx = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx")

    >>> from pyXLMS.parser import read_msannika
    >>> csms_from_tsv = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt")

    >>> from pyXLMS.parser import read_msannika
    >>> crosslinks_from_tsv = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt")

    >>> from pyXLMS.parser import read_msannika
    >>> csms_and_crosslinks_from_pdresult = read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult")
    """
    ## check input
    _ok = check_input(modifications, "modifications", dict, float)
    _ok = check_input(format, "format", str)
    _ok = check_input(sep, "sep", str)

    ## default parser
    if modification_parser is None:
        modification_parser = pyxlms_modification_str_parser

    ## helper functions
    def parse_modification_str(
        sequence: str,
        modification_str: str,
        modifications: Dict[str, float] = modifications,
    ) -> Dict[int, Tuple[str, float]]:
        mods = [mod.strip() for mod in modification_str.split(";")]
        parsed_mods = dict()
        for mod in mods:
            mod_type = mod.split("(")[1].split(")")[0].strip()
            mod_pos = mod.split("(")[0].strip()
            if mod_type not in modifications:
                raise KeyError(
                    f"Unable to find modification {mod_type} in the set of provided modifications. "
                    + "Please pass the full set of expected modifications to the parser."
                )
            if "Nterm" in mod_pos or "N-Term" in mod_pos:
                parsed_mods[0] = (mod_type, modifications[mod_type])
            elif "Cterm" in mod_pos or "C-Term" in mod_pos:
                parsed_mods[len(sequence)] = (mod_type, modifications[mod_type])
            else:
                parsed_mods[int(mod_pos[1:])] = (mod_type, modifications[mod_type])
        return parsed_mods

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
        data_objects = None
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
                data_objects = [pd.read_csv(input, sep=sep, low_memory=False)]
            elif file_extension == ".xlsx":
                data_objects = [pd.read_excel(input, engine="openpyxl")]
            elif file_extension == ".pdresult":
                data_objects = __read_msannika_pdresult(input)
            else:
                raise ValueError(
                    f"Detected file extension {file_extension} is not supported! Input file has to be a valid file with extension '.csv', '.tsv' or '.xlsx'!"
                )
        elif format in ["csv", "tsv", "txt", "xlsx"]:
            if format == "xlsx":
                data_objects = [pd.read_excel(input, engine="openpyxl")]
            else:
                data_objects = [pd.read_csv(input, sep=sep, low_memory=False)]
        elif format == "pdresult":
            if not isinstance(input, str):
                raise TypeError(
                    "Can't read pdResult files from a file-like object/stream. Please provide the filename/path instead!"
                )
            data_objects = __read_msannika_pdresult(input)
        else:
            raise ValueError(
                f"Provided input format {format} is not supported! Input format has to be of type 'csv', 'tsv' or 'xlsx'!"
            )
        if data_objects is None:
            raise RuntimeError(
                "Something went wrong while reading the file! Please file a bug report!"
            )
        # this should be impossible, but check here for pyright
        if not isinstance(data_objects, list):
            raise RuntimeError(
                "Something went wrong while reading the file! Please file a bug report!"
            )
        for data in data_objects:
            # this should be impossible, but check here for pyright
            if not isinstance(data, pd.DataFrame):
                raise RuntimeError(
                    "Something went wrong while reading the file! Please file a bug report!"
                )
            ## detect input file type
            col_names = data.columns.values.tolist()
            is_crosslink_dataframe = "Best CSM Score" in col_names
            ## process data
            if is_crosslink_dataframe:
                for i, row in tqdm(
                    data.iterrows(),
                    total=data.shape[0],
                    desc="Reading MS Annika crosslinks...",
                ):
                    # create crosslink
                    crosslink = create_crosslink(
                        peptide_a=format_sequence(str(row["Sequence A"])),
                        xl_position_peptide_a=int(row["Position A"]),
                        proteins_a=[
                            protein.strip()
                            for protein in str(row["Accession A"]).split(";")
                        ],
                        xl_position_proteins_a=[
                            int(position)
                            for position in str(row["In protein A"]).split(";")
                        ],
                        decoy_a=get_bool_from_value(row["Decoy"]),
                        peptide_b=format_sequence(str(row["Sequence B"])),
                        xl_position_peptide_b=int(row["Position B"]),
                        proteins_b=[
                            protein.strip()
                            for protein in str(row["Accession B"]).split(";")
                        ],
                        xl_position_proteins_b=[
                            int(position)
                            for position in str(row["In protein B"]).split(";")
                        ],
                        decoy_b=get_bool_from_value(row["Decoy"]),
                        score=float(row["Best CSM Score"]),
                    )
                    crosslinks.append(crosslink)
            else:
                for i, row in tqdm(
                    data.iterrows(),
                    total=data.shape[0],
                    desc="Reading MS Annika CSMs...",
                ):
                    # create csm
                    csm = create_csm(
                        peptide_a=format_sequence(str(row["Sequence A"])),
                        modifications_a=parse_modification_str(
                            format_sequence(str(row["Sequence A"]).strip()),
                            str(row["Modifications A"]).strip(),
                        ),
                        xl_position_peptide_a=int(row["Crosslinker Position A"]),
                        proteins_a=[
                            protein.strip()
                            for protein in str(row["Accession A"]).split(";")
                        ],
                        xl_position_proteins_a=[
                            int(position) + int(row["Crosslinker Position A"])
                            for position in str(row["A in protein"]).split(";")
                        ],
                        pep_position_proteins_a=[
                            int(position) + 1
                            for position in str(row["A in protein"]).split(";")
                        ],
                        score_a=float(row["Score Alpha"]),
                        decoy_a=not get_bool_from_value(str(row["Alpha T/D"])),
                        peptide_b=format_sequence(str(row["Sequence B"])),
                        modifications_b=parse_modification_str(
                            format_sequence(str(row["Sequence B"]).strip()),
                            str(row["Modifications B"]).strip(),
                        ),
                        xl_position_peptide_b=int(row["Crosslinker Position B"]),
                        proteins_b=[
                            protein.strip()
                            for protein in str(row["Accession B"]).split(";")
                        ],
                        xl_position_proteins_b=[
                            int(position) + int(row["Crosslinker Position B"])
                            for position in str(row["B in protein"]).split(";")
                        ],
                        pep_position_proteins_b=[
                            int(position) + 1
                            for position in str(row["B in protein"]).split(";")
                        ],
                        score_b=float(row["Score Beta"]),
                        decoy_b=not get_bool_from_value(str(row["Beta T/D"])),
                        score=float(row["Combined Score"]),
                        spectrum_file=str(row["Spectrum File"]).strip(),
                        scan_nr=int(row["First Scan"]),
                        charge=int(row["Charge"]),
                        rt=float(row["RT [min]"]) * 60.0,
                        im_cv=float(row["Compensation Voltage"]),
                    )
                    csms.append(csm)
    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="MS Annika",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )
