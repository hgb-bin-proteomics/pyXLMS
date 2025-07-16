#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
import pandas as pd
from tqdm import tqdm

from ..data import check_input
from ..data import create_crosslink
from ..data import create_csm
from ..data import create_parser_result
from ..constants import SCOUT_MODIFICATION_MAPPING
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


def detect_scout_filetype(
    data: pd.DataFrame,
) -> Literal["scout_csms_unfiltered", "scout_csms_filtered", "scout_xl"]:
    r"""Detects the Scout-related source of the data.

    Detects whether the input data is unfiltered crosslink-spectrum-matches, filtered crosslink-spectrum-matches,
    or crosslinks from Scout.

    Parameters
    ----------
    data : pd.DataFrame
        The input data originating from Scout.

    Returns
    -------
    str
        "scout_csms_unfiltered" if a Scout unfiltered CSMs file was read, "scout_csms_filtered" if a Scout filtered CSMs file was read,
        "scout_xl" if a Scout crosslink/residue pair result file was read.

    Raises
    ------
    ValueError
        If the data source could not be determined.

    Examples
    --------
    >>> from pyXLMS.parser import detect_scout_filetype
    >>> import pandas as pd
    >>> df1 = pd.read_csv("data/scout/Cas9_Unfiltered_CSMs.csv")
    >>> detect_scout_filetype(df1)
    'scout_csms_unfiltered'

    >>> from pyXLMS.parser import detect_scout_filetype
    >>> import pandas as pd
    >>> df2 = pd.read_csv("data/scout/Cas9_Filtered_CSMs.csv")
    >>> detect_scout_filetype(df2)
    'scout_csms_filtered'

    >>> from pyXLMS.parser import detect_scout_filetype
    >>> import pandas as pd
    >>> df3 = pd.read_csv("data/scout/Cas9_Residue_Pairs.csv")
    >>> detect_scout_filetype(df3)
    'scout_xl'
    """
    ## check input
    _ok = check_input(data, "data", pd.DataFrame)

    col_names = data.columns.values.tolist()
    if "ScanNumber" in col_names:
        return "scout_csms_unfiltered"
    if "Scan" in col_names:
        return "scout_csms_filtered"
    if "CSM count" in col_names:
        return "scout_xl"

    raise ValueError(
        "Could not infer data source, are you sure you read a Scout result file?"
    )

    return "err"


def parse_modifications_from_scout_sequence(
    seq: str,
    crosslink_position: int,
    crosslinker: str,
    crosslinker_mass: float,
    modifications: Dict[str, Tuple[str, float]] = SCOUT_MODIFICATION_MAPPING,
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[int, Tuple[str, float]]:
    r"""Parse post-translational-modifications from a Scout peptide sequence.

    Parses post-translational-modifications (PTMs) from a Scout peptide sequence,
    for example "M(+15.994900)LASAGELQKGNELALPSK".

    Parameters
    ----------
    seq : str
        The Scout sequence string.
    crosslink_position : int
        Position of the crosslinker in the sequence (1-based).
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float
        Monoisotopic delta mass of the crosslink modification.
    modifications: dict of str, float, default = ``constants.SCOUT_MODIFICATION_MAPPING``
        Mapping of modification names to modification masses.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    dict of int, tuple
        The ``pyXLMS`` specific modifications object, a dictionary that maps positions to their corresponding modifications and their
        monoisotopic masses.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed (only if ``verbose = 2``).
    KeyError
        If an unknown modification is encountered.

    Examples
    --------
    >>> from pyXLMS.parser import parse_modifications_from_scout_sequence
    >>> seq = "M(+15.994900)LASAGELQKGNELALPSK"
    >>> parse_modifications_from_scout_sequence(seq, 10, "DSS", 138.06808)
    {10: ('DSS', 138.06808), 1: ('Oxidation', 15.994915)}

    >>> from pyXLMS.parser import parse_modifications_from_scout_sequence
    >>> seq = "KIEC(+57.021460)FDSVEISGVEDR"
    >>> parse_modifications_from_scout_sequence(seq, 1, "DSS", 138.06808)
    {1: ('DSS', 138.06808), 4: ('Carbamidomethyl', 57.021464)}
    """
    # clean seq
    sequence = seq.strip()
    # init parsed modifications dict
    parsed_modifications = {crosslink_position: (crosslinker, crosslinker_mass)}
    # parse modifications from sequence
    pos = 0
    current_mod = ""
    for i, aa in enumerate(sequence):
        if aa.isupper():
            pos += 1
            current_mod = ""
        else:
            current_mod += aa
            if (i + 1 >= len(sequence)) or (sequence[i + 1].isupper()):
                mod_key = current_mod.strip("()").strip()
                if mod_key not in modifications:
                    raise KeyError(
                        f"Key {mod_key} not found in parameter 'modifications'. Are you missing a modification?"
                    )
                if pos in parsed_modifications:
                    err_str = (
                        f"Modification at position {pos} already exists!\n"
                        f"Sequence: {sequence}, Crosslink position: {crosslink_position}"
                    )
                    if verbose == 1:
                        warnings.warn(RuntimeWarning(err_str))
                    elif verbose == 2:
                        raise RuntimeError(err_str)
                    t1 = parsed_modifications[pos][0] + "," + modifications[mod_key][0]
                    t2 = parsed_modifications[pos][1] + modifications[mod_key][1]
                    parsed_modifications[pos] = (t1, t2)
                else:
                    parsed_modifications[pos] = modifications[mod_key]
    return parsed_modifications


def __read_scout_csms_unfiltered(
    data: pd.DataFrame,
    crosslinker: str,
    crosslinker_mass: float,
    parse_modifications: bool,
    modifications: Dict[str, Tuple[str, float]],
    verbose: Literal[0, 1, 2],
) -> List[Dict[str, Any]]:
    r"""Reads crosslink-spectrum-matches from a Scout unfiltered CSMs result.

    Parameters
    ----------
    data : pandas.DataFrame
        The Scout unfiltered CSMs result data.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float
        Monoisotopic delta mass of the crosslink modification.
    parse_modifications : bool
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications : dict of str, tuple
        Mapping of Scout sequence elements (e.g. ``"+15.994900"``) and modifications (e.g ``"Oxidation of Methionine"``)
        to their modifications (e.g. ``("Oxidation", 15.994915)``).
    verbose : 0, 1, or 2
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    list of dict
        The read crosslink-spectrum-matches.

    Notes
    -----
    This function should not be called directly, it is called from ``read_scout()``.
    """
    csms = list()
    xl = data.dropna(axis=0, subset=["AlphaPeptide", "BetaPeptide"])
    for i, row in tqdm(
        xl.iterrows(), total=xl.shape[0], desc="Reading Scout unfiltered CSMs..."
    ):
        csm = create_csm(
            peptide_a=format_sequence(str(row["AlphaPeptide"])),
            modifications_a=parse_modifications_from_scout_sequence(
                str(row["AlphaPeptide"]),
                int(row["AlphaPos"]) + 1,
                crosslinker,
                crosslinker_mass,
                modifications,
                verbose,
            )
            if parse_modifications
            else None,
            xl_position_peptide_a=int(row["AlphaPos"]) + 1,
            proteins_a=[
                protein.strip() for protein in str(row["AlphaMappings"]).split(";")
            ],
            xl_position_proteins_a=None,
            pep_position_proteins_a=None,
            score_a=float(row["AlphaScore"]),
            decoy_a=str(row["Class"]).strip() in ["FullDecoy", "BetaTarget"],
            peptide_b=format_sequence(str(row["BetaPeptide"])),
            modifications_b=parse_modifications_from_scout_sequence(
                str(row["BetaPeptide"]),
                int(row["BetaPos"]) + 1,
                crosslinker,
                crosslinker_mass,
                modifications,
                verbose,
            )
            if parse_modifications
            else None,
            xl_position_peptide_b=int(row["BetaPos"]) + 1,
            proteins_b=[
                protein.strip() for protein in str(row["BetaMappings"]).split(";")
            ],
            xl_position_proteins_b=None,
            pep_position_proteins_b=None,
            score_b=float(row["BetaScore"]),
            decoy_b=str(row["Class"]) in ["FullDecoy", "AlphaTarget"],
            score=float(row["XLScore"]),
            spectrum_file=str(row["FileName"]).strip(),
            scan_nr=int(row["ScanNumber"]),
            charge=int(row["Charge"]),
            rt=None,
            im_cv=None,
            additional_information={
                "ClassificationScore": float(row["ClassificationScore"]),
                "XlinkxAlpha": float(row["XlinkxAlpha"]),
                "XlinkxBeta": float(row["XlinkxBeta"]),
                "XlinkxScore": float(row["XlinkxScore"]),
                "PoissonScore": float(row["PoissonScore"]),
            },
        )
        csms.append(csm)
    return csms


def __read_scout_csms_filtered(
    data: pd.DataFrame,
    crosslinker: str,
    crosslinker_mass: float,
    parse_modifications: bool,
    modifications: Dict[str, Tuple[str, float]],
    verbose: Literal[0, 1, 2],
) -> List[Dict[str, Any]]:
    r"""Reads crosslink-spectrum-matches from a Scout filtered CSMs result.

    Parameters
    ----------
    data : pandas.DataFrame
        The Scout filtered CSMs result data.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float
        Monoisotopic delta mass of the crosslink modification.
    parse_modifications : bool
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications : dict of str, tuple
        Mapping of Scout sequence elements (e.g. ``"+15.994900"``) and modifications (e.g ``"Oxidation of Methionine"``)
        to their modifications (e.g. ``("Oxidation", 15.994915)``).
    verbose : 0, 1, or 2
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    list of dict
        The read crosslink-spectrum-matches.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed (only if ``verbose = 2``).
    KeyError
        If an unknown modification is encountered.

    Notes
    -----
    This function should not be called directly, it is called from ``read_scout()``.
    """

    ## helper functions
    def str_contains(s: str, contains: List[str]) -> bool:
        for subs in contains:
            if subs in s:
                return True
        return False

    def parse_modifications_fn(
        row: pd.Series,
        alpha: bool,
        crosslinker: str,
        crosslinker_mass: float,
        modifications: Dict[str, Tuple[str, float]] = SCOUT_MODIFICATION_MAPPING,
        verbose: Literal[0, 1, 2] = 1,
    ) -> Dict[int, Tuple[str, float]]:
        sequence = (
            str(row["Alpha peptide"]).strip()
            if alpha
            else str(row["Beta peptide"]).strip()
        )
        crosslink_position = (
            int(row["Alpha peptide position"])
            if alpha
            else int(row["Beta peptide position"])
        )
        parsed_modifications = {crosslink_position: (crosslinker, crosslinker_mass)}
        if alpha and bool(pd.isna(row["Alpha modification(s)"])):
            return parsed_modifications
        if not alpha and bool(pd.isna(row["Beta modification(s)"])):
            return parsed_modifications
        mods = (
            str(row["Alpha modification(s)"]).split(";")
            if alpha
            else str(row["Beta modification(s)"]).split(";")
        )
        for mod in mods:
            rpos = mod.split("(")[0].strip()
            mod_key = mod.split("(")[1].rstrip(")").strip()
            pos = -1
            if str_contains(
                rpos.lower(),
                [
                    "nterm",
                    "nterminal",
                    "nterminus",
                    "n-term",
                    "n-terminal",
                    "n-terminus",
                ],
            ):
                pos = 0
            elif str_contains(
                rpos.lower(),
                [
                    "cterm",
                    "cterminal",
                    "cterminus",
                    "c-term",
                    "c-terminal",
                    "c-terminus",
                ],
            ):
                pos = len(sequence)
            else:
                pos = int(rpos[1:])
            if mod_key not in modifications:
                raise KeyError(
                    f"Key {mod_key} not found in parameter 'modifications'. Are you missing a modification?"
                )
            if pos in parsed_modifications:
                err_str = (
                    f"Modification at position {pos} already exists!\n"
                    f"CSM Scan Number: {int(row['Scan'])}!\n"
                    f"Sequence: {sequence}, Crosslink position: {crosslink_position}, Modifications: {';'.join(mods)}"
                )
                if verbose == 1:
                    warnings.warn(RuntimeWarning(err_str))
                elif verbose == 2:
                    raise RuntimeError(err_str)
                t1 = parsed_modifications[pos][0] + "," + modifications[mod_key][0]
                t2 = parsed_modifications[pos][1] + modifications[mod_key][1]
                parsed_modifications[pos] = (t1, t2)
            else:
                parsed_modifications[pos] = modifications[mod_key]
        return parsed_modifications

    ## create csms
    csms = list()
    xl = data.dropna(axis=0, subset=["Alpha peptide", "Beta peptide"])
    for i, row in tqdm(
        xl.iterrows(), total=xl.shape[0], desc="Reading Scout filtered CSMs..."
    ):
        csm = create_csm(
            peptide_a=format_sequence(str(row["Alpha peptide"])),
            modifications_a=parse_modifications_fn(
                row,
                True,
                crosslinker,
                crosslinker_mass,
                modifications,
                verbose,
            )
            if parse_modifications
            else None,
            xl_position_peptide_a=int(row["Alpha peptide position"]),
            proteins_a=[
                protein.strip()
                for protein in str(row["Alpha protein mapping(s)"]).split(";")
            ],
            xl_position_proteins_a=[
                int(pos) for pos in str(row["Alpha protein(s) position(s)"]).split(";")
            ],
            pep_position_proteins_a=[
                int(pos) - int(row["Alpha peptide position"]) + 1
                for pos in str(row["Alpha protein(s) position(s)"]).split(";")
            ],
            score_a=None,
            decoy_a=get_bool_from_value(row["IsDecoy"]),
            peptide_b=format_sequence(str(row["Beta peptide"])),
            modifications_b=parse_modifications_fn(
                row,
                False,
                crosslinker,
                crosslinker_mass,
                modifications,
                verbose,
            )
            if parse_modifications
            else None,
            xl_position_peptide_b=int(row["Beta peptide position"]),
            proteins_b=[
                protein.strip()
                for protein in str(row["Beta protein mapping(s)"]).split(";")
            ],
            xl_position_proteins_b=[
                int(pos) for pos in str(row["Beta protein(s) position(s)"]).split(";")
            ],
            pep_position_proteins_b=[
                int(pos) - int(row["Beta peptide position"]) + 1
                for pos in str(row["Beta protein(s) position(s)"]).split(";")
            ],
            score_b=None,
            decoy_b=get_bool_from_value(row["IsDecoy"]),
            score=float(row["Score"]),
            spectrum_file=str(row["File"]).strip(),
            scan_nr=int(row["Scan"]),
            charge=int(row["Precursor charge"]),
            rt=None,
            im_cv=None,
        )
        csms.append(csm)
    return csms


def __read_scout_crosslinks(data: pd.DataFrame) -> List[Dict[str, Any]]:
    r"""Reads crosslinks from a Scout crosslink/residue pair result.

    Parameters
    ----------
    data : pandas.DataFrame
        The Scout crosslink/residue pair result data.

    Returns
    -------
    list of dict
        The read crosslinks.

    Notes
    -----
    This function should not be called directly, it is called from ``read_scout()``.
    """
    crosslinks = list()
    xl = data.dropna(axis=0, subset=["Alpha peptide", "Beta peptide"])
    for i, row in tqdm(
        xl.iterrows(), total=xl.shape[0], desc="Reading Scout crosslinks..."
    ):
        crosslink = create_crosslink(
            peptide_a=format_sequence(str(row["Alpha peptide"])),
            xl_position_peptide_a=int(row["Alpha peptide position"]),
            proteins_a=[
                protein.strip()
                for protein in str(row["Alpha protein mapping(s)"]).split(";")
            ],
            xl_position_proteins_a=[
                int(pos) for pos in str(row["Alpha protein(s) position(s)"]).split(";")
            ],
            decoy_a=get_bool_from_value(row["IsDecoy"]),
            peptide_b=format_sequence(str(row["Beta peptide"])),
            xl_position_peptide_b=int(row["Beta peptide position"]),
            proteins_b=[
                protein.strip()
                for protein in str(row["Beta protein mapping(s)"]).split(";")
            ],
            xl_position_proteins_b=[
                int(pos) for pos in str(row["Beta protein(s) position(s)"]).split(";")
            ],
            decoy_b=get_bool_from_value(row["IsDecoy"]),
            score=float(row["Score"]),
        )
        crosslinks.append(crosslink)
    return crosslinks


def read_scout(
    files: str | List[str] | BinaryIO,
    crosslinker: str,
    crosslinker_mass: Optional[float] = None,
    parse_modifications: bool = True,
    modifications: Dict[str, Tuple[str, float]] = SCOUT_MODIFICATION_MAPPING,
    sep: str = ",",
    decimal: str = ".",
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[str, Any]:
    r"""Read a Scout result file.

    Reads a Scout filtered or unfiltered crosslink-spectrum-matches result file or crosslink/residue pair result file in ``.csv`` format
    and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the Scout result file(s) or a file-like object/stream.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float, or None, default = None
        Monoisotopic delta mass of the crosslink modification. If the crosslinker is
        defined in parameter "modifications" this can be omitted.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications : dict of str, tuple, default = ``constants.SCOUT_MODIFICATION_MAPPING``
        Mapping of Scout sequence elements (e.g. ``"+15.994900"``) and modifications (e.g ``"Oxidation of Methionine"``)
        to their modifications (e.g. ``("Oxidation", 15.994915)``).
    sep : str, default = ","
        Seperator used in the ``.csv`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.
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
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslinks or crosslink-spectrum-matches.
    KeyError
        If the specified crosslinker could not be found/mapped.
    TypeError
        If parameter verbose was not set correctly.

    Warnings
    --------
    - When reading unfiltered crosslink-spectrum-matches, no protein crosslink positions or protein peptide positions are
      available, as these are not reported. If needed they should be annotated with ``transform.reannotate_positions()``.
    - When reading filtered crosslink-spectrum-matches, Scout does not report if the individual peptides in a crosslink are
      from the target or decoy database. The parser assumes that both peptides from a target crosslink-spectrum-match are
      from the target database, and vice versa, that both peptides are from the decoy database if it is a decoy crosslink-spectrum-match.
      This leads to only TT and DD matches, which needs to be considered for FDR estimation.
    - When reading crosslinks / residue pairs, Scout does not report if the individual peptides in a crosslink are from the
      target or decoy database. The parser assumes that both peptides from a target crosslink are from the target database,
      and vice versa, that both peptides are from the decoy database if it is a decoy crosslink. This leads to only TT and DD
      matches, which needs to be considered for FDR estimation.

    Examples
    --------
    >>> from pyXLMS.parser import read_scout
    >>> csms_unfiltered = read_scout("data/scout/Cas9_Unfiltered_CSMs.csv")

    >>> from pyXLMS.parser import read_scout
    >>> csms_filtered = read_scout("data/scout/Cas9_Filtered_CSMs.csv")

    >>> from pyXLMS.parser import read_scout
    >>> crosslinks = read_scout("data/scout/Cas9_Residue_Pairs.csv")
    """
    ## check input
    _ok = check_input(crosslinker, "crosslinker", str)
    _ok = (
        check_input(crosslinker_mass, "crosslinker_mass", float)
        if crosslinker_mass is not None
        else True
    )
    _ok = check_input(parse_modifications, "parse_modifications", bool)
    _ok = check_input(modifications, "modifications", dict, tuple)
    _ok = check_input(sep, "sep", str)
    _ok = check_input(decimal, "decimal", str)
    _ok = check_input(verbose, "verbose", int)
    if crosslinker_mass is None:
        if crosslinker not in modifications:
            if parse_modifications:
                raise KeyError(
                    "Cannot infer crosslinker mass because crosslinker is not defined in "
                    "parameter 'modifications'. Please specify crosslinker mass manually!"
                )
            else:
                crosslinker_mass = 0.0
        else:
            crosslinker_mass = modifications[crosslinker][1]
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")

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
        data = pd.read_csv(input, sep=sep, decimal=decimal, low_memory=False)
        ## detect input file type
        scout_file_type = detect_scout_filetype(data)
        ## process data
        if scout_file_type == "scout_csms_unfiltered":
            csms += __read_scout_csms_unfiltered(
                data,
                crosslinker,
                crosslinker_mass,
                parse_modifications,
                modifications,
                verbose,
            )
        elif scout_file_type == "scout_csms_filtered":
            csms += __read_scout_csms_filtered(
                data,
                crosslinker,
                crosslinker_mass,
                parse_modifications,
                modifications,
                verbose,
            )
        else:
            crosslinks += __read_scout_crosslinks(data)

    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="Scout",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )
