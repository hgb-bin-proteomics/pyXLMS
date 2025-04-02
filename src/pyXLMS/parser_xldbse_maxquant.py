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
from .parser_util import get_bool_from_value

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List


def __parse_maxquant_modification_str(
    seq: str, 
    crosslink_position: int,
    crosslinker: str,
    crosslinker_mass: float,
    modifications: Dict[str, float] = MODIFICATIONS
) -> Dict[int, Tuple[str, float]]:
    """
    Examples
    --------
    >>> from pyXLMS.parser_xldbse_maxquant import __parse_maxquant_modification_str
    >>> seq = "_VVDELVKVM(Oxidation (M))GR_"
    >>> __parse_maxquant_modification_str(seq, 2, "DSS", 138.06808)
    {2: ('DSS', 3.0), 9: ('Oxidation', 15.994915)}

    >>> from pyXLMS.parser_xldbse_maxquant import __parse_maxquant_modification_str
    >>> seq = "_VVDELVKVM(Oxidation (M))GRM(Oxidation (M))_"
    >>> __parse_maxquant_modification_str(seq, 2, "DSS", 138.06808)
    {2: ('DSS', 3.0), 9: ('Oxidation', 15.994915), 12: ('Oxidation', 15.994915)}

    >>> from pyXLMS.parser_xldbse_maxquant import __parse_maxquant_modification_str
    >>> seq = "_M(Oxidation (M))VVDELVKVM(Oxidation (M))GRM(Oxidation (M))_"
    >>> __parse_maxquant_modification_str(seq, 2, "DSS", 138.06808)
    {2: ('DSS', 3.0), 1: ('Oxidation', 15.994915), 10: ('Oxidation', 15.994915), 13: ('Oxidation', 15.994915)}
    """
    parsed_modifications = {crosslink_position: (crosslinker, crosslinker_mass)}
    ## start parse seq
    split_seq = seq.split("_")
    if len(split_seq) != 3:
        raise RuntimeError(
            f"Could not parse sequence {seq}. Is the sequence correctly formatted?"
            )
    _n_term = split_seq[0].strip() # don't use nterm mods because I don't know how they are formatted
    internal = split_seq[1].strip()
    _c_term = split_seq[2].strip() # don't use cterm mods because I don't know how they are formatted
    ## end parse seq
    is_mod = 0
    current_pos = 0
    current_mod = ""
    for aa in internal:
        if is_mod == 0:
            if aa == "(":
                is_mod += 1
            else:
                current_pos += 1
        else:
            if aa == "(":
                is_mod += 1
            elif aa == ")":
                is_mod -= 1
            else:
                current_mod += aa
            if is_mod == 0:
                if current_pos in parsed_modifications:
                    raise RuntimeError()
                else:
                    current_mod = current_mod.split()[0]
                    if current_mod not in modifications:
                        raise KeyError()
                    else:
                        parsed_modifications[current_pos] = (current_mod, modifications[current_mod])
                current_mod = ""
    return parsed_modifications


def read_maxquant(
    files: str | List[str] | BinaryIO,
    crosslinker: str,
    crosslinker_mass: Optional[float] = None,
    decoy_prefix: str = "REV__",
    modifications: Dict[str, float] = MODIFICATIONS,
    sep: str = "\t",
) -> Dict[str, Any]:
    """Read a MaxQuant result file.

    Reads a MaxQuant crosslink-spectrum-matches result file "crosslinkMsms.txt" in ``.txt`` (tab delimited) format
    and returns a ``parser_result``.

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
        If one of the found post-translational-modifications could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_maxquant
    >>> csms_from_xlsx = read_maxquant("data/maxquant/run1/crosslinkMsms.txt")
    """
    ## check input
    _ok = check_input(crosslinker, "crosslinker", str)
    _ok = check_input(crosslinker_mass, "crosslinker_mass", float) if crosslinker_mass is not None else True
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
        for i, row in tqdm(xl.iterrows(), total=xl.shape[0], desc="Reading MaxQuant CSMs..."):
            # create csm
            csm = create_csm()
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
        If one of the found post-translational-modifications could not be found/mapped.

    Examples
    --------
    >>> from pyXLMS.parser import read_maxlynx
    >>> csms_from_xlsx = read_maxlynx("data/maxquant/run1/crosslinkMsms.txt")
    """
    return read_maxquant(files, decoy_prefix, modifications, sep)
