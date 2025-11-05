#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import io
import zipfile
import pandas as pd
from tqdm import tqdm

from ..data import check_input
from ..data import create_csm
from ..data import create_parser_result
from ..constants import AMINO_ACIDS
from ..constants import MODIFICATIONS
from ..constants import MEROX_MODIFICATION_MAPPING
from .util import __serialize_pandas_series

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List


MEROX_COLNAMES = [
    "Score",
    "m/z",
    "Charge",
    "M+H+",
    "Calculated Mass",
    "Deviation in ppm",
    "Peptide 1",
    "Protein 1",
    "From",
    "To",
    "Peptide2",
    "Protein 2",
    "From.1",
    "To.1",
    "Scan number",
    "is Selected in Table",
    "Candidate identifier",
    "Folder Number",
    "Retention time in sec",
    "miscellaneous",
    "best linkage position peptide 1",
    "best linkage position peptide 2",
    "All linkage positions",
    "Spectrum UUID",
    "local False discovery rate",
    "shortest distance in pdb",
    "Light/Heavy(1/2)",
    "pepScore1",
    "pepScore2",
    "xLinkScore",
    "resultId",
    "MS1intensity",
    "finalScoreComponent",
]


def __read_merox_file(
    file: str | BinaryIO, sep: str = ";", decimal: str = "."
) -> pd.DataFrame:
    r"""Helper function to read MeroX files into pandas DataFrames.

    Reads MeroX files into pandas DataFrames independent of input format. Accepts
    both ``.csv`` and ``.zhrm`` files.

    Parameters
    ----------
    file : str, or file stream
        The name/path of the MeroX result file or a file-like object/stream.
    sep : str, default = ";"
        Seperator used in the ``.csv`` or ``.zhrm`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.

    Returns
    -------
    pd.DataFrame
        The MeroX result file as a pandas DataFrame.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    # safety
    if not isinstance(file, str):
        file.seek(0)
    # this shifts file pointer
    if zipfile.is_zipfile(file):
        with zipfile.ZipFile(file, "r") as f:
            return pd.read_csv(
                io.BytesIO(f.read("Result.csv")),
                header=None,
                names=MEROX_COLNAMES,
                sep=sep,
                decimal=decimal,
                low_memory=False,
            )
    if not isinstance(file, str):
        file.seek(0)
    return pd.read_csv(file, sep=sep, decimal=decimal, low_memory=False)


def __get_merox_sequence(
    sequence: str,
    parse_modifications: bool = True,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
) -> str:
    r"""Helper funtion to read the MeroX peptide sequence.

    Reads the MeroX peptide sequence and removes prefix and suffix square brackets,
    then replaces modification symbols with their respective amino acids. If
    ``parse_modifications = True`` non-resolved symbols will raise an error,
    otherwise the symbol will be used in the sequence - even if it is not a valid
    amino acid sequence anymore.

    Parameters
    ----------
    sequence : str
        The MeroX sequence e.g. from column "Peptide 1" or "Peptide2".
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples.

    Returns
    -------
    str
        The parsed sequence.

    Raises
    ------
    KeyError
        If one of the symbols in the sequence could not be resolved/an unknown modification is encountered.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    parsed_seq = ""
    for amino_acid in sequence.lstrip("[").rstrip("]").strip():
        if amino_acid in modifications:
            parsed_seq += modifications[amino_acid]["Amino Acid"]
        elif amino_acid in AMINO_ACIDS:
            parsed_seq += amino_acid
        else:
            if parse_modifications:
                raise KeyError(
                    f"Key {amino_acid} not found in parameter 'modifications'. Are you missing a modification?"
                )
            else:
                parsed_seq += amino_acid
    return parsed_seq


def __get_merox_modifications(
    sequence: str,
    crosslink_position: int,
    crosslinker: str,
    crosslinker_mass: float,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
) -> Dict[int, Tuple[str, float]]:
    r"""Helper function to read modifications from a MeroX sequence.

    Parses post-translational-modifications from the given MeroX sequence.

    Parameters
    ----------
    sequence : str
        The MeroX sequence e.g. from column "Peptide 1" or "Peptide2".
    crosslink_position : int
        Position of the crosslinker in the sequence (1-based).
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float
        Monoisotopic delta mass of the crosslink modification.
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples.

    Returns
    -------
    dict of int, tuple of str, float
        The pyXLMS specfic modification representation of the parsed modifications.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed.
    KeyError
        If an unknown modification is encountered.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    parsed_modifications = {crosslink_position: (crosslinker, crosslinker_mass)}
    for i, amino_acid in enumerate(sequence.lstrip("[").rstrip("]").strip()):
        if amino_acid in modifications:
            if i + 1 not in parsed_modifications:
                parsed_modifications[i + 1] = modifications[amino_acid]["Modification"]
            else:
                raise RuntimeError(f"Modification at position {i + 1} already exists!")
        elif amino_acid in AMINO_ACIDS:
            pass
        else:
            raise KeyError(
                f"Key {amino_acid} not found in parameter 'modifications'. Are you missing a modification?"
            )
    return parsed_modifications


def __get_merox_position(position_str: str) -> int:
    r"""Helper function to extract the peptide crosslink position from MeroX.

    Parameters
    ----------
    position_str : str
        The position string from MeroX e.g. from column "best linkage position peptide 1".

    Returns
    -------
    int
        The parsed peptide crosslink position.

    Raises
    ------
    RuntimeError
        If the position could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    position = None
    try:
        position = int(position_str[1:])
    except Exception as _e:
        pass
    if position is None:
        raise RuntimeError(f"Could not parse position from {position_str}!")
    return position


def __get_merox_protein(proteins: str) -> List[str]:
    r"""Helper function to extract the protein accession from MeroX.

    Parses the leading protein accession from the MeroX protein string. Additional
    proteins are not parsed as they do not have corresponding protein crosslink
    positions.

    Parameters
    ----------
    proteins : str
        The protein string from MeroX e.g. from column "Protein 1".

    Returns
    -------
    list of str
        A list containing the single parsed protein accession.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    return [proteins.split("(>")[0].strip().lstrip(">")]


def __get_merox_scan_number(scan_nr_and_file: str) -> int:
    r"""Helper function to parse the scan number from MeroX.

    Parses the scan number from the MeroX scan number string. Assumes that scan
    number and spectrum file are delimted by the tilde (wave) sign.

    Parameters
    ----------
    scan_nr_and_file : str
        The scan number string from MeroX e.g. from the column "Scan number".

    Returns
    -------
    int
        The parsed scan number.

    Raises
    ------
    RuntimeError
        If the scan number could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    scan_nr = None
    try:
        scan_nr = int(scan_nr_and_file.split("~")[0])
    except Exception as _e:
        pass
    if scan_nr is None:
        raise RuntimeError(f"Could not scan number from {scan_nr_and_file}!")
    return scan_nr


def __get_merox_spectrum_file(scan_nr_and_file: str) -> str:
    r"""Helper function to parse the spectrum file from MeroX.

    Parses the spectrum file from the MeroX scan number string. Assumes that scan
    number and spectrum file are delimted by the tilde (wave) sign.

    Parameters
    ----------
    scan_nr_and_file : str
        The scan number string from MeroX e.g. from the column "Scan number".

    Returns
    -------
    str
        The parsed spectrum file.

    Raises
    ------
    RuntimeError
        If the spectrum file could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    spectrum_file = None
    try:
        spectrum_file = scan_nr_and_file.split("~")[1]
    except Exception as _e:
        pass
    if spectrum_file is None:
        raise RuntimeError(f"Could not spectrum file name from {scan_nr_and_file}!")
    return spectrum_file


def read_merox(
    files: str | List[str] | BinaryIO,
    crosslinker: str,
    crosslinker_mass: Optional[float] = None,
    decoy_prefix: str = "REV__",
    parse_modifications: bool = True,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
    sep: str = ";",
    decimal: str = ".",
) -> Dict[str, Any]:
    r"""Read a MeroX result file.

    Reads a MeroX crosslink-spectrum-matches result file in ``.csv`` or ``.zhrm`` format
    and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MeroX result file(s) or a file-like object/stream.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float, or None, default = None
        Monoisotopic delta mass of the crosslink modification. If the crosslinker is
        defined in ``constants.MODIFICATIONS`` this can be omitted.
    decoy_prefix : str, default = "REV\_\_"
        The prefix that indicates that a protein is from the decoy database.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples.
    sep : str, default = ";"
        Seperator used in the ``.csv`` or ``.zhrm`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.

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

    Notes
    -----
    Uses ``pepScore1`` as the score for the alpha peptide, ``pepScore2`` as the score of the
    beta peptide, and ``Score`` as the score of the crosslink-spectrum-match.

    Warnings
    --------
    MeroX only reports a single protein crosslink position per peptide, for ambiguous peptides
    only the crosslink position of the first matching protein is reported. All matching proteins can be
    retrieved via ``additional_information``, however not their corresponding crosslink positions. For this
    reason it is recommended to use ``transform.reannotate_positions()`` to correctly annotate all crosslink
    positions for all peptides if that is important for downstream analysis. Additionally, please note that
    target and decoy information is derived based off the protein accession and parameter ``decoy_prefix``.
    By default, MeroX only reports target matches that are above the desired FDR.

    Examples
    --------
    >>> from pyXLMS.parser import read_merox
    >>> csms_from_csv = read_merox(
    ...     "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.csv", crosslinker="DSS"
    ... )

    >>> from pyXLMS.parser import read_merox
    >>> csms_from_zhrm = read_merox(
    ...     "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.zhrm", crosslinker="DSS"
    ... )
    """
    ## check input
    _ok = check_input(crosslinker, "crosslinker", str)
    _ok = (
        check_input(crosslinker_mass, "crosslinker_mass", float)
        if crosslinker_mass is not None
        else True
    )
    _ok = check_input(decoy_prefix, "decoy_prefix", str)
    _ok = check_input(parse_modifications, "parse_modifications", bool)
    _ok = check_input(modifications, "modifications", dict, dict)
    _ok = check_input(sep, "sep", str)
    _ok = check_input(decimal, "decimal", str)
    if crosslinker_mass is None:
        if crosslinker not in MODIFICATIONS:
            if parse_modifications:
                raise KeyError(
                    "Cannot infer crosslinker mass because crosslinker is unknown. "
                    "Please specify crosslinker mass manually!"
                )
            else:
                crosslinker_mass = 0.0
        else:
            crosslinker_mass = MODIFICATIONS[crosslinker]

    ## data structures
    csms = list()

    ## handle input
    if not isinstance(files, list):
        inputs = [files]
    else:
        inputs = files

    for input in inputs:
        ## reading data
        data = __read_merox_file(input, sep=sep, decimal=decimal)
        for i, row in tqdm(
            data.iterrows(),
            total=data.shape[0],
            desc="Reading MeroX CSMs...",
        ):
            # create crosslink
            csm = create_csm(
                peptide_a=__get_merox_sequence(
                    str(row["Peptide 1"]), parse_modifications, modifications
                ),
                modifications_a=__get_merox_modifications(
                    str(row["Peptide 1"]),
                    __get_merox_position(str(row["best linkage position peptide 1"])),
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                )
                if parse_modifications
                else None,
                xl_position_peptide_a=__get_merox_position(
                    str(row["best linkage position peptide 1"])
                ),
                proteins_a=__get_merox_protein(str(row["Protein 1"])),
                xl_position_proteins_a=[
                    int(row["From"])
                    + __get_merox_position(str(row["best linkage position peptide 1"]))
                    - 1
                ],
                pep_position_proteins_a=[int(row["From"])],
                score_a=float(row["pepScore1"]),
                decoy_a=__get_merox_protein(str(row["Protein 1"]))[0].startswith(
                    decoy_prefix
                ),
                peptide_b=__get_merox_sequence(
                    str(row["Peptide2"]), parse_modifications, modifications
                ),
                modifications_b=__get_merox_modifications(
                    str(row["Peptide2"]),
                    __get_merox_position(str(row["best linkage position peptide 2"])),
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                )
                if parse_modifications
                else None,
                xl_position_peptide_b=__get_merox_position(
                    str(row["best linkage position peptide 2"])
                ),
                proteins_b=__get_merox_protein(str(row["Protein 2"])),
                xl_position_proteins_b=[
                    int(row["From.1"])
                    + __get_merox_position(str(row["best linkage position peptide 2"]))
                    - 1
                ],
                pep_position_proteins_b=[int(row["From.1"])],
                score_b=float(row["pepScore2"]),
                decoy_b=__get_merox_protein(str(row["Protein 2"]))[0].startswith(
                    decoy_prefix
                ),
                score=float(row["Score"]),
                spectrum_file=__get_merox_spectrum_file(str(row["Scan number"])),
                scan_nr=__get_merox_scan_number(str(row["Scan number"])),
                charge=int(row["Charge"]),
                rt=float(row["Retention time in sec"]),
                im_cv=None,
                additional_information={
                    "source": __serialize_pandas_series(row),
                    "xLinkScore": row["xLinkScore"],
                    "Protein 1": row["Protein 1"],
                    "Protein 2": row["Protein 2"],
                    "MS1intensity": row["MS1intensity"],
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
        search_engine="MeroX",
        csms=csms,
        crosslinks=None,
    )
