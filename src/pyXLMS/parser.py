#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
import pandas as pd
from os.path import splitext

from .data import create_crosslink
from .data import create_csm
from .data import create_parser_result
from .constants import AMINO_ACIDS

from typing import BinaryIO
from typing import Dict
from typing import Literal
from typing import Any

def format_sequence(str: sequence, remove_non_aa = True, remove_lower = True) -> str:
    """Formats the given amino acid sequence into common represenation.

    The given amino acid sequence is re-formatted by converting all amino acids to upper case and optionally removing non-encoding and
    lower case characters.

    Parameters
    ----------
    sequence : str
        The amino acid sequence that should be formatted. Post-translational-modifications can be included in lower case but will
        be removed.
    remove_non_aa : bool, default = True
        Whether or not to remove characters that do not encode amino acids.
    remove_lower : bool, default = True
        Whether or not to remove lower case characters, this should be true if the amino acid sequence encodes post-translational-modifications
        in lower case.

    Returns
    -------
    str
        The formatted sequence.

    Examples
    --------
    >>>from pyXLMS.parser import format_sequence
    >>>format_sequence("PEP[K]TIDE")
    'PEPKTIDE'

    >>>from pyXLMS.parser import format_sequence
    >>>format_sequence("PEPKdssoTIDE")
    'PEPKTIDE'

    >>>from pyXLMS.parser import format_sequence
    >>>format_sequence("peptide", remove_lower = False)
    'PEPTIDE'
    """
    fmt_seq = ""
    for aa in sequence:
        if aa.isupper():
            if aa not in AMINO_ACIDS:
                if remove_non_aa:
                    continue
                else:
                    warnings.warn(f"The sequence {sequence} contains non-valid characters.", RuntimeWarning)
            fmt_seq += aa
        elif remove_lower:
            continue
        else:
            if aa.upper() not in AMINO_ACIDS:
                if remove_non_aa:
                    continue
                else:
                    warnings.warn(f"The sequence {sequence} contains non-valid characters.", RuntimeWarning)
            fmt_seq += aa.upper()
    return fmt_seq


def get_bool_from_value(value: Any) -> bool:
    """Parse a bool value from the given input.

    Tries to parse a boolean value from the given input object. If the object is of instance ``bool`` it will return the object, if it is of
    instance ``int`` it will return ``True`` if the object is ``1`` or ``False`` if the object is ``0``, any other number will raise a
    ``ValueError``. If the object is of instance ``str`` it will return ``True`` if the lower case version contains the letter ``t`` and
    otherwise ``False``. If the object is none of these types a ``ValueError`` will be raised.

    Parameters
    ----------
    value: Any
        The value to parse from.

    Returns
    -------
    bool
        The parsed boolean value.

    Raises
    ------
    ValueError
        If the object could not be parsed to bool.

    Examples
    --------
    >>>from pyXLMS.parser import get_bool_from_value
    >>>get_bool_from_value(0)
    False

    >>>from pyXLMS.parser import get_bool_from_value
    >>>get_bool_from_value("T")
    True
    """
    if isinstance(value, bool):
        return value
    elif isinstance(value, int):
        if value in [0, 1]:
            return bool(value)
        else:
            raise ValueError(f"Cannot parse bool value from the given input {value}.")
    elif isinstance(value, str):
        return "t" in value.lower()
    else:
        raise ValueError(f"Cannot parse bool value from the given input {value}.")
    return False


## TODO
def read_custom():
    return


## WIP
def read_msannika(
    input: str | BinaryIO,
    format: Literal["auto", "csv", "tsv", "xlsx"] = "auto",
    sep: str = "\t",
) -> Dict[str, Any]:
    """Read an MS Annika result file.

    Reads an MS Annika crosslink-spectrum-matches result file or crosslink result file in ``.csv`` or ``.xlsx`` format and returns a
    ``parser_result``.

    Parameters
    ----------
    input : str or file stream
        The name/path of the MS Annika result file or a file-like object/stream.
    format : "auto", "csv", "tsv", or "xlsx", default = "auto"
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
    RuntimeError
        If the file could not be read.
    """
    ## reading data
    data = None
    if format == "auto" and not isinstance(input, str):
        raise ValueError(
            "Can't detect format for file-like objects. Please specify format manually!"
        )
    # and isinstance specified for type checking
    if format == "auto" and isinstance(input, str):
        file_extension = splitext(input)
        if file_extension == ".tsv" or file_extension == ".csv":
            data = pd.read_csv(input, sep=sep)
        elif file_extension == ".xlsx":
            data = pd.read_excel(input, engine="openpyxl")
        else:
            raise ValueError(
                f"Detected file extension {file_extension} is not supported! Input file has to be a valid file with extension '.csv', '.tsv' or '.xlsx'!"
            )
    elif format in ["csv", "tsv", "xlsx"]:
        if format == "xlsx":
            data = pd.read_excel(input, engine="openpyxl")
        else:
            data = pd.read_csv(input, sep=sep)
    else:
        raise ValueError(
            f"Provided input format {format} is not supported! Input format has to be of type 'csv', 'tsv' or 'xlsx'!"
        )
    if data is None:
        raise RuntimeError(
            "Something went wrong while reading the file! Please file a bug report!"
        )
    ## detect input file type
    col_names = data.columns.values.tolist()
    is_crosslink_dataframe = "# CSMs" in col_names
    ## process data
    crosslinks = list()
    csms = list()
    if is_crosslink_dataframe:
        for i, row in data.iterrows():
            # create crosslink
            crosslink = create_crosslink(
                peptide_a=format_sequence(str(row["Sequence A"]).strip()),
                xl_position_peptide_a=int(row["Position A"]),
                proteins_a=[protein.strip() for protein in str(row["Accession A"]).split(";")],
                xl_position_proteins_a=[int(position) for position in str(row["In protein A"]).split(";")],
                decoy_a=get_bool_from_value(row["Decoy"]),
                peptide_b=format_sequence(str(row["Sequence B"]).strip()),
                xl_position_peptide_b=int(row["Position A"]),
                proteins_b=[protein.strip() for protein in str(row["Accession B"]).split(";")],
                xl_position_proteins_b=[int(position) for position in str(row["In protein B"]).split(";")],
                decoy_b=get_bool_from_value(row["Decoy"]),
                score=float(row["Best CSM score"]),
            )
            crosslinks.append(crosslink)
    else:
        for i, row in data.iterrows():
            # create csm
            csm = create_csm(
                peptide_a="",
                modifications_a={},
                xl_position_peptide_a=0,
                proteins_a=[""],
                xl_position_proteins_a=[0],
                pep_position_proteins_a=[0],
                score_a=0.0,
                decoy_a=False,
                peptide_b="",
                modifications_b={},
                xl_position_peptide_b=0,
                proteins_b=[""],
                xl_position_proteins_b=[0],
                pep_position_proteins_b=[0],
                score_b=0.0,
                decoy_b=False,
                score=0.0,
                spectrum_file="",
                scan_nr=0,
                charge=0,
                rt=0.0,
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
        search_engine="MS Annika",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )


## TODO
def read():
    return
