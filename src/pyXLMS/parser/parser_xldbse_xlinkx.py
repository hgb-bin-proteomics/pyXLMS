#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import sqlite3
import warnings
import pandas as pd
from tqdm import tqdm
from os.path import splitext

from ..data import check_input
from ..data import create_crosslink
from ..data import create_csm
from ..data import create_parser_result
from ..constants import MODIFICATIONS
from .util import format_sequence
from .util import get_bool_from_value

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __read_xlinkx_pdresult(filename: str) -> List[pd.DataFrame]:
    r"""Read an XlinkX pdResult file and convert it to standard XlinkX result tables.

    Reads an XlinkX pdResult file and converts it to standard XlinkX CSM and crosslink
    tables as pandas DataFrame objects.

    Parameters
    ----------
    filename : str
        Filename/path of the pdResult file.

    Returns
    -------
    list of pd.DataFrame
        The converted CSM and crosslink tables in XlinkX format.

    Notes
    -----
    This function should not be called directly, it is called from ``read_xlinkx()``.
    """
    conn = sqlite3.connect(filename)
    csms = pd.read_sql_query("SELECT * FROM CSMs", conn)
    dcsms = pd.read_sql_query("SELECT * FROM DecoyCSMs", conn)
    xls = pd.read_sql_query("SELECT * FROM Crosslinks", conn)
    conn.close()
    column_mapping_csms = {
        "SequenceA": "Sequence A",
        "SequenceB": "Sequence B",
        "ModificationsA": "Modifications A",
        "ModificationsB": "Modifications B",
        "CrosslinkerPositionA": "Crosslinker Position A",
        "CrosslinkerPositionB": "Crosslinker Position B",
        "ProteinAccessionA": "Protein Accession A",
        "ProteinAccessionB": "Protein Accession B",
        "LeadingProteinPositionA": "Leading Protein Position A",
        "LeadingProteinPositionB": "Leading Protein Position B",
        "IsDecoy": "Is Decoy",
        "XlinkXScore": "XlinkX Score",
        "SpectrumFileName": "Spectrum File",
        "FirstScan": "First Scan",
        "ChargeState": "Charge",
        "RetentionTime": "RT [min]",
    }
    column_mapping_xls = {
        "SequenceA": "Sequence A",
        "SequenceB": "Sequence B",
        "ModificationsA": "Modifications A",
        "ModificationsB": "Modifications B",
        "Crosslinker": "Crosslinker",
        "PositionA": "Position A",
        "PositionB": "Position B",
        "AccessionA": "Accession A",
        "AccessionB": "Accession B",
        "IsDecoy": "Is Decoy",
        "MaxXlinkXScore": "Max. XlinkX Score",
    }
    csms.rename(columns=column_mapping_csms, inplace=True)
    dcsms.rename(columns=column_mapping_csms, inplace=True)
    xls.rename(columns=column_mapping_xls, inplace=True)
    csms.drop(
        columns=list(
            set(csms.columns.values.tolist()) - set(list(column_mapping_csms.values()))
        ),
        inplace=True,
    )
    dcsms.drop(
        columns=list(
            set(dcsms.columns.values.tolist()) - set(list(column_mapping_csms.values()))
        ),
        inplace=True,
    )
    xls.drop(
        columns=list(
            set(xls.columns.values.tolist()) - set(list(column_mapping_xls.values()))
        ),
        inplace=True,
    )
    return [csms, dcsms, xls]


def read_xlinkx(
    files: str | List[str] | BinaryIO,
    decoy: Optional[bool] = None,
    parse_modifications: bool = True,
    modifications: Dict[str, float] = MODIFICATIONS,
    format: Literal["auto", "csv", "txt", "tsv", "xlsx", "pdresult"] = "auto",
    sep: str = "\t",
    decimal: str = ".",
    ignore_errors: bool = False,
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[str, Any]:
    r"""Read an XlinkX result file.

    Reads an XlinkX crosslink-spectrum-matches result file or crosslink result file in ``.csv`` or ``.xlsx`` format,
    or both from a ``.pdResult`` file from Proteome Discover, and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the XlinkX result file(s) or a file-like object/stream.
    decoy : bool, or None
        Default decoy value to use if no decoy value is found. Only used if the "Is Decoy" column is not found
        in the supplied data.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
    format : "auto", "csv", "tsv", "txt", "xlsx", or "pdresult", default = "auto"
        The format of the result file. ``"auto"`` is only available if the name/path to the XlinkX result file is given.
    sep : str, default = "\t"
        Seperator used in the ``.csv`` or ``.tsv`` file. Parameter is ignored if the file is in ``.xlsx`` or ``.pdResult`` format.
    decimal : str, default = "."
        Character to recognize as decimal point. Parameter is ignored if the file is in ``.xlsx`` or ``.pdResult`` format.
    ignore_errors : bool, default = False
        If missing crosslink positions should raise an error or not. Setting this to True will suppress the ``RuntimeError``
        for the crosslink position not being able to be parsed for at least one of the crosslinks. For these cases the crosslink
        position will be set to 100 000.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    ValueError
        If the input format is not supported or cannot be inferred.
    TypeError
        If parameter verbose was not set correctly.
    TypeError
        If the pdResult file is provided in the wrong format.
    RuntimeError
        If the crosslink position could not be parsed for at least one of the crosslinks.
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslinks or crosslink-spectrum-matches.
    KeyError
        If one of the found post-translational-modifications could not be found/mapped.

    Warnings
    --------
    XlinkX does not report if the individual peptides in a crosslink are from the target or decoy database.
    The parser assumes that both peptides from a target crosslink are from the target database, and vice versa,
    that both peptides are from the decoy database if it is a decoy crosslink. This leads to only TT and DD matches,
    which needs to be considered for FDR estimation. This applies to both crosslinks **and** crosslink-spectrum-matches.

    Examples
    --------
    >>> from pyXLMS.parser import read_xlinkx
    >>> csms_from_xlsx = read_xlinkx("data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_CSMs.xlsx")

    >>> from pyXLMS.parser import read_xlinkx
    >>> crosslinks_from_xlsx = read_xlinkx("data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_Crosslinks.xlsx")

    >>> from pyXLMS.parser import read_xlinkx
    >>> csms_from_tsv = read_xlinkx("data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_CSMs.txt")

    >>> from pyXLMS.parser import read_xlinkx
    >>> crosslinks_from_tsv = read_xlinkx("data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_Crosslinks.txt")

    >>> from pyXLMS.parser import read_xlinkx
    >>> csms_and_crosslinks_from_pdresult = read_xlinkx("data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3.pdResult")
    """
    ## check input
    _ok = check_input(decoy, "decoy", bool) if decoy is not None else True
    _ok = check_input(parse_modifications, "parse_modifications", bool)
    _ok = check_input(modifications, "modifications", dict, float)
    _ok = check_input(format, "format", str)
    _ok = check_input(sep, "sep", str)
    _ok = check_input(decimal, "decimal", str)
    _ok = check_input(ignore_errors, "ignore_errors", bool)
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")

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

    def get_crosslink_position_from_peptide_seq(
        sequence: str,
        crosslinker: str,
        modifications: str,
        ignore_errors: bool = False,
        verbose: Literal[0, 1, 2] = 1,
    ) -> int:
        seq = str(sequence).strip()
        xl = str(crosslinker).strip()
        mods = [mod.strip() for mod in str(modifications).split(";")]
        for i, aa in enumerate(seq):
            if aa == "[":
                return i + 1
        for mod in mods:
            if xl in mod:
                return int(mod.split("[")[1].split("]")[0][1:])
        if verbose == 2 or not ignore_errors:
            raise RuntimeError(
                f"Could not parse crosslink position from sequence: {seq}, or modifications {modifications}!"
            )
        if verbose == 1:
            warnings.warn(
                RuntimeWarning(
                    f"Could not parse crosslink position from sequence: {seq}, or modifications {modifications}!"
                )
            )
        return 100000

    def adjust_crosslink_position(
        position: int,
        sequence: str,
    ) -> int:
        if position == 0:
            return 1
        if position > len(sequence.strip()):
            return len(sequence.strip())
        return position

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
                data_objects = [
                    pd.read_csv(input, sep=sep, decimal=decimal, low_memory=False)
                ]
            elif file_extension == ".xlsx":
                data_objects = [pd.read_excel(input, engine="openpyxl")]
            elif file_extension == ".pdresult":
                data_objects = __read_xlinkx_pdresult(input)
            else:
                raise ValueError(
                    f"Detected file extension {file_extension} is not supported! Input file has to be a valid file with extension '.csv', '.tsv', '.xlsx', or '.pdResult'!"
                )
        elif format in ["csv", "tsv", "txt", "xlsx"]:
            if format == "xlsx":
                data_objects = [pd.read_excel(input, engine="openpyxl")]
            else:
                data_objects = [
                    pd.read_csv(input, sep=sep, decimal=decimal, low_memory=False)
                ]
        elif format == "pdresult":
            if not isinstance(input, str):
                raise TypeError(
                    "Can't read pdResult files from a file-like object/stream. Please provide the filename/path instead!"
                )
            data_objects = __read_xlinkx_pdresult(input)
        else:
            raise ValueError(
                f"Provided input format {format} is not supported! Input format has to be of type 'csv', 'tsv', 'xlsx', or 'pdresult'!"
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
            is_crosslink_dataframe = "Max. XlinkX Score" in col_names
            ## process data
            if is_crosslink_dataframe:
                for i, row in tqdm(
                    data.iterrows(),
                    total=data.shape[0],
                    desc="Reading XlinkX crosslinks...",
                ):
                    # create crosslink
                    crosslink = create_crosslink(
                        peptide_a=format_sequence(str(row["Sequence A"])),
                        xl_position_peptide_a=get_crosslink_position_from_peptide_seq(
                            str(row["Sequence A"]),
                            str(row["Crosslinker"]),
                            str(row["Modifications A"]),
                            ignore_errors,
                            verbose,
                        ),
                        proteins_a=[
                            protein.strip()
                            for protein in str(row["Accession A"]).split(";")
                        ],
                        xl_position_proteins_a=[
                            int(position)
                            for position in str(row["Position A"]).split(";")
                        ],
                        decoy_a=get_bool_from_value(row["Is Decoy"]),
                        peptide_b=format_sequence(str(row["Sequence B"])),
                        xl_position_peptide_b=get_crosslink_position_from_peptide_seq(
                            str(row["Sequence B"]),
                            str(row["Crosslinker"]),
                            str(row["Modifications B"]),
                            ignore_errors,
                            verbose,
                        ),
                        proteins_b=[
                            protein.strip()
                            for protein in str(row["Accession B"]).split(";")
                        ],
                        xl_position_proteins_b=[
                            int(position)
                            for position in str(row["Position B"]).split(";")
                        ],
                        decoy_b=get_bool_from_value(row["Is Decoy"]),
                        score=float(row["Max. XlinkX Score"]),
                    )
                    crosslinks.append(crosslink)
            else:
                for i, row in tqdm(
                    data.iterrows(),
                    total=data.shape[0],
                    desc="Reading XlinkX CSMs...",
                ):
                    # create csm
                    csm = create_csm(
                        peptide_a=format_sequence(str(row["Sequence A"])),
                        modifications_a=parse_modification_str(
                            format_sequence(str(row["Sequence A"]).strip()),
                            str(row["Modifications A"]).strip(),
                        )
                        if parse_modifications
                        else None,
                        xl_position_peptide_a=adjust_crosslink_position(
                            int(row["Crosslinker Position A"]),
                            format_sequence(str(row["Sequence A"]).strip()),
                        ),
                        proteins_a=[
                            protein.strip()
                            for protein in str(row["Protein Accession A"]).split(";")
                        ],
                        xl_position_proteins_a=[
                            int(position)
                            for position in str(
                                row["Leading Protein Position A"]
                            ).split(";")
                        ],
                        pep_position_proteins_a=[
                            int(position)
                            - adjust_crosslink_position(
                                int(row["Crosslinker Position A"]),
                                format_sequence(str(row["Sequence A"]).strip()),
                            )
                            + 1
                            for position in str(
                                row["Leading Protein Position A"]
                            ).split(";")
                        ],
                        score_a=None,
                        decoy_a=get_bool_from_value(row["Is Decoy"])
                        if "Is Decoy" in col_names
                        else decoy,
                        peptide_b=format_sequence(str(row["Sequence B"])),
                        modifications_b=parse_modification_str(
                            format_sequence(str(row["Sequence B"]).strip()),
                            str(row["Modifications B"]).strip(),
                        )
                        if parse_modifications
                        else None,
                        xl_position_peptide_b=adjust_crosslink_position(
                            int(row["Crosslinker Position B"]),
                            format_sequence(str(row["Sequence B"])),
                        ),
                        proteins_b=[
                            protein.strip()
                            for protein in str(row["Protein Accession B"]).split(";")
                        ],
                        xl_position_proteins_b=[
                            int(position)
                            for position in str(
                                row["Leading Protein Position B"]
                            ).split(";")
                        ],
                        pep_position_proteins_b=[
                            int(position)
                            - adjust_crosslink_position(
                                int(row["Crosslinker Position B"]),
                                format_sequence(str(row["Sequence B"])),
                            )
                            + 1
                            for position in str(
                                row["Leading Protein Position B"]
                            ).split(";")
                        ],
                        score_b=None,
                        decoy_b=get_bool_from_value(row["Is Decoy"])
                        if "Is Decoy" in col_names
                        else decoy,
                        score=float(row["XlinkX Score"]),
                        spectrum_file=str(row["Spectrum File"]).strip(),
                        scan_nr=int(row["First Scan"]),
                        charge=int(row["Charge"]),
                        rt=float(row["RT [min]"]) * 60.0,
                        im_cv=None,
                    )
                    csms.append(csm)
    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="XlinkX",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )
