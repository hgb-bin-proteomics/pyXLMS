#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from os.path import splitext

from .data import check_input
from .data import create_crosslink
from .data import create_csm
from .data import create_parser_result
from .constants import MODIFICATIONS
from .parser_util import format_sequence
from .parser_util import get_bool_from_value

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


def read_msannika(
    files: str | List[str] | BinaryIO,
    modifications: Dict[str, float] = MODIFICATIONS,
    format: Literal["auto", "csv", "txt", "tsv", "xlsx"] = "auto",
    sep: str = "\t",
) -> Dict[str, Any]:
    """Read an MS Annika result file.

    Reads an MS Annika crosslink-spectrum-matches result file or crosslink result file in ``.csv`` or ``.xlsx`` format and returns a
    ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MS Annika result file(s) or a file-like object/stream.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
    format : "auto", "csv", "tsv", "txt", or "xlsx", default = "auto"
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
    """
    ## check input
    _ok = check_input(modifications, "modifications", dict, float)

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
            if "Nterm" in mod_pos:
                parsed_mods[0] = (mod_type, modifications[mod_type])
            elif "Cterm" in mod_pos:
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
        data = None
        if format == "auto" and not isinstance(input, str):
            raise ValueError(
                "Can't detect format for file-like objects. Please specify format manually!"
            )
        # and isinstance specified for type checking
        if format == "auto" and isinstance(input, str):
            file_extension = splitext(input)[1]
            if (
                file_extension == ".txt"
                or file_extension == ".tsv"
                or file_extension == ".csv"
            ):
                data = pd.read_csv(input, sep=sep)
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
        if is_crosslink_dataframe:
            for i, row in data.iterrows():
                # create crosslink
                crosslink = create_crosslink(
                    peptide_a=format_sequence(str(row["Sequence A"]).strip()),
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
                    peptide_b=format_sequence(str(row["Sequence B"]).strip()),
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
            for i, row in data.iterrows():
                # create csm
                csm = create_csm(
                    peptide_a=format_sequence(str(row["Sequence A"]).strip()),
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
                    peptide_b=format_sequence(str(row["Sequence B"]).strip()),
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
