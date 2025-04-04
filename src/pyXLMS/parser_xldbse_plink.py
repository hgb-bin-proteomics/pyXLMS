#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from tqdm import tqdm

from .data import check_input
from .data import create_csm
from .data import create_parser_result
from .constants import MODIFICATIONS
from .parser_util import format_sequence

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List


def read_plink(
    files: str | List[str] | BinaryIO,
    crosslinker: str,
    crosslinker_mass: Optional[float] = None,
    decoy_prefix: str = "REV__",
    modifications: Dict[str, float] = MODIFICATIONS,
    sep: str = "\t",
) -> Dict[str, Any]:
    """Read a pLink 2 result file.

    Reads a pLink 2 crosslink-spectrum-matches result file "*cross-linked_spectra.csv"
    in ``.csv`` (comma delimited) format and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MaxQuant result file(s) or a file-like object/stream.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float, or None, default = None
        Monoisotopic delta mass of the crosslink modification. If the crosslinker is
        defined in parameter "modifications" this can be omitted.
    decoy_prefix : str, default = "REV__"
        The prefix that indicates that a protein is from the decoy database.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
    sep : str, default = "\t"
        Seperator used in the ``.txt`` file.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches.
    KeyError
        If the specified crosslinker could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_maxquant
    >>> csms_from_xlsx = read_maxquant("data/maxquant/run1/crosslinkMsms.txt")
    """
    ## check input
    _ok = check_input(crosslinker, "crosslinker", str)
    _ok = (
        check_input(crosslinker_mass, "crosslinker_mass", float)
        if crosslinker_mass is not None
        else True
    )
    _ok = check_input(decoy_prefix, "decoy_prefix", str)
    _ok = check_input(modifications, "modifications", dict, float)
    _ok = check_input(sep, "sep", str)
    if crosslinker_mass is None:
        if crosslinker not in modifications:
            raise KeyError(
                "Cannot infer crosslinker mass because crosslinker is not defined in "
                "parameter 'modifications'. Please specify crosslinker mass manually!"
            )
        else:
            crosslinker_mass = modifications[crosslinker]

    ## data structures
    csms = list()

    ## handle input
    if not isinstance(files, list):
        inputs = [files]
    else:
        inputs = files

    ## process data
    for input in inputs:
        data = pd.read_csv(input, sep=sep, low_memory=False)
        xl = data.dropna(axis=0, subset=["Proteins2"])
        for i, row in tqdm(
            xl.iterrows(), total=xl.shape[0], desc="Reading MaxQuant CSMs..."
        ):
            # preprocess proteins
            protein_a = (
                str(row["Proteins1"]).split("(")[0].strip()
                if "(" in str(row["Proteins1"])
                else str(row["Proteins1"])
            )
            protein_b = (
                str(row["Proteins2"]).split("(")[0].strip()
                if "(" in str(row["Proteins2"])
                else str(row["Proteins2"])
            )
            # create csm
            csm = create_csm(
                peptide_a=format_sequence(str(row["Sequence1"])),
                modifications_a=parse_modifications_from_maxquant_sequence(
                    str(row["Modified sequence1"]),
                    int(row["Peptide index of Crosslink 1"]),
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                ),
                xl_position_peptide_a=int(row["Peptide index of Crosslink 1"]),
                proteins_a=[
                    protein_a.strip()
                    if protein_a.strip()[: len(decoy_prefix)] != decoy_prefix
                    else protein_a.strip()[len(decoy_prefix) :]
                ],
                xl_position_proteins_a=[int(row["Protein index of Crosslink 1"])],
                pep_position_proteins_a=[
                    int(row["Protein index of Crosslink 1"])
                    - int(row["Peptide index of Crosslink 1"])
                    + 1
                ],
                score_a=float(row["Partial score 1"]),
                decoy_a=decoy_prefix in str(row["Proteins1"]),
                peptide_b=format_sequence(str(row["Sequence2"])),
                modifications_b=parse_modifications_from_maxquant_sequence(
                    str(row["Modified sequence2"]),
                    int(row["Peptide index of Crosslink 2"]),
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                ),
                xl_position_peptide_b=int(row["Peptide index of Crosslink 2"]),
                proteins_b=[
                    protein_b.strip()
                    if protein_b.strip()[: len(decoy_prefix)] != decoy_prefix
                    else protein_b.strip()[len(decoy_prefix) :]
                ],
                xl_position_proteins_b=[int(row["Protein index of Crosslink 2"])],
                pep_position_proteins_b=[
                    int(row["Protein index of Crosslink 2"])
                    - int(row["Peptide index of Crosslink 2"])
                    + 1
                ],
                score_b=float(row["Partial score 2"]),
                decoy_b=decoy_prefix in str(row["Proteins2"]),
                score=float(row["Score"]),
                spectrum_file=str(row["Raw file"]).strip(),
                scan_nr=int(row["Scan number"]),
                charge=int(row["Charge"]),
                rt=None,
                im_cv=None,
                additional_information={
                    "Proteins1": str(row["Proteins1"]).strip(),
                    "Proteins2": str(row["Proteins2"]).strip(),
                    "Delta score": float(row["Delta score"]),
                },
            )
            csms.append(csm)
    ## check results
    if len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="MaxQuant",
        csms=csms,
        crosslinks=None,
    )


def read_maxlynx(
    files: str | List[str] | BinaryIO,
    crosslinker: str,
    crosslinker_mass: Optional[float] = None,
    decoy_prefix: str = "REV__",
    modifications: Dict[str, float] = MODIFICATIONS,
    sep: str = "\t",
) -> Dict[str, Any]:
    """Read a MaxLynx result file.

    Reads a MaxLynx crosslink-spectrum-matches result file "crosslinkMsms.txt" in ``.txt`` (tab delimited) format
    and returns a ``parser_result``. This is an alias for the MaxQuant reader.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MaxLynx result file(s) or a file-like object/stream.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float, or None, default = None
        Monoisotopic delta mass of the crosslink modification. If the crosslinker is
        defined in parameter "modifications" this can be omitted.
    decoy_prefix : str, default = "REV__"
        The prefix that indicates that a protein is from the decoy database.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
    sep : str, default = "\t"
        Seperator used in the ``.txt`` file.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches.
    KeyError
        If the specified crosslinker could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_maxlynx
    >>> csms_from_xlsx = read_maxlynx("data/maxquant/run1/crosslinkMsms.txt")
    """
    return read_maxquant(
        files, crosslinker, crosslinker_mass, decoy_prefix, modifications, sep
    )
