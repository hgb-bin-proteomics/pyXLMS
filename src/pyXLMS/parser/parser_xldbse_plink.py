#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import io
import warnings
import pandas as pd
from tqdm import tqdm

from ..data import check_input
from ..data import create_csm
from ..data import create_crosslink
from ..data import create_parser_result
from ..constants import MODIFICATIONS
from .util import format_sequence
from .util import __serialize_pandas_series

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


def __parse_modifications_from_plink_modifications_str(
    seq: str,
    mod_str: Optional[str | float],
    crosslinker: str,
    modifications: Dict[str, float] = MODIFICATIONS,
    verbose: Literal[0, 1, 2] = 1,
) -> Tuple[Dict[int, Tuple[str, float]], Dict[int, Tuple[str, float]]]:
    r"""Parse post-translational-modifications from a pLink modification string.

    Parses post-translational-modifications (PTMs) from a pLink modification string,
    for example "Carbamidomethyl[C](4);Oxidation[M](23)".

    Parameters
    ----------
    seq : str
        The pLink crosslink sequence string.
    mod_str : str, float, or None
        The pLink modification value, as string or float. Can be None.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    tuple of dict of int, tuple
        The ``pyXLMS`` specific modification objects, dictionaries that map positions to their corresponding modifications and their
        monoisotopic masses. The first object (index 0) corresponds to the modifications of the first peptide, the second object (index 1)
        corresponds to the modifications of the second peptide.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed (only for ``verbose = 2``).
    KeyError
        If an unknown modification is encountered.

    Notes
    -----
    This function should not be called directly, it is called from ``read_plink()``.
    """
    modifications_a = dict()
    modifications_b = dict()
    xl_pos_a = int(seq.split("-")[0].split("(")[1].split(")")[0])
    xl_pos_b = int(seq.split("-")[1].split("(")[1].split(")")[0])
    if crosslinker in modifications:
        modifications_a[xl_pos_a] = (crosslinker, modifications[crosslinker])
        modifications_b[xl_pos_b] = (crosslinker, modifications[crosslinker])
    else:
        raise KeyError(
            f"Key {crosslinker} not found in parameter 'modifications'. Are you missing a modification?"
        )
    if mod_str is None:
        return (modifications_a, modifications_b)
    if isinstance(mod_str, float) and pd.isna(mod_str):
        return (modifications_a, modifications_b)
    mod_str = str(mod_str).strip()
    if mod_str == "nan":
        return (modifications_a, modifications_b)
    mods = mod_str.split(";")
    for mod in mods:
        mod_desc = mod.split("[")[0].strip()
        if mod_desc not in modifications:
            raise KeyError(
                f"Key {mod_desc} not found in parameter 'modifications'. Are you missing a modification?"
            )
        mod_pos = int(mod.split("(")[1].split(")")[0])
        if mod_pos > len(seq.split("-")[0]):
            mod_pos = mod_pos - len(seq.split("-")[0])
            if mod_pos in modifications_b:
                if verbose == 2:
                    raise RuntimeError(
                        f"Modification at position {mod_pos} already exists!"
                    )
                if verbose == 1:
                    warnings.warn(
                        RuntimeWarning(
                            f"Modification at position {mod_pos} already exists!"
                        )
                    )
                t1 = modifications_b[mod_pos][0] + "," + mod_desc
                t2 = modifications_b[mod_pos][1] + modifications[mod_desc]
                modifications_b[mod_pos] = (t1, t2)
            else:
                modifications_b[mod_pos] = (mod_desc, modifications[mod_desc])
        else:
            if mod_pos in modifications_a:
                if verbose == 2:
                    raise RuntimeError(
                        f"Modification at position {mod_pos} already exists!"
                    )
                if verbose == 1:
                    warnings.warn(
                        RuntimeWarning(
                            f"Modification at position {mod_pos} already exists!"
                        )
                    )
                t1 = modifications_a[mod_pos][0] + "," + mod_desc
                t2 = modifications_a[mod_pos][1] + modifications[mod_desc]
                modifications_a[mod_pos] = (t1, t2)
            else:
                modifications_a[mod_pos] = (mod_desc, modifications[mod_desc])
    return (modifications_a, modifications_b)


def __parse_proteins_and_position_from_plink(
    seq: str,
    proteins: str,
) -> Dict[str, Any]:
    r"""Parses proteins and positions from pLink results.

    Parses proteins, as well as peptide and crosslink positions from a pLink crosslink sequence
    and protein string.

    Parameters
    ----------
    seq : str
        The pLink crosslink sequence string.
    proteins : str
        The pLink proteins string.

    Returns
    -------
    dict of str, Any
        A dictionary with the following keys and information:
        ``xl_pos_a``, ``proteins_a``, ``proteins_a_xl_positions``, ``proteins_a_pep_positions``,
        ``xl_pos_b``, ``proteins_b``, ``proteins_b_xl_positions``, ``proteins_b_pep_positions``.

    Notes
    -----
    This function should not be called directly, it is called from ``read_plink()``.
    """
    xl_pos_a = int(seq.split("-")[0].split("(")[1].split(")")[0])
    xl_pos_b = int(seq.split("-")[1].split("(")[1].split(")")[0])
    # proteins a
    proteins_set_a = set()
    proteins_a = list()
    proteins_a_xl_positions = list()
    proteins_a_pep_positions = list()
    # proteins b
    proteins_set_b = set()
    proteins_b = list()
    proteins_b_xl_positions = list()
    proteins_b_pep_positions = list()
    # find unique
    proteins = proteins.strip().rstrip("/")
    for protein_pair in proteins.split("/"):
        protein_a = protein_pair.split("-")[0].strip()
        protein_b = protein_pair.split("-")[1].strip()
        proteins_set_a.add(protein_a)
        proteins_set_b.add(protein_b)
    # get proteins a
    for protein in sorted(proteins_set_a):
        acc = protein.split("(")[0]
        pos = int(protein.split("(")[1].split(")")[0])
        proteins_a.append(acc)
        proteins_a_xl_positions.append(pos)
        proteins_a_pep_positions.append(pos - xl_pos_a + 1)
    # get proteins b
    for protein in sorted(proteins_set_b):
        acc = protein.split("(")[0]
        pos = int(protein.split("(")[1].split(")")[0])
        proteins_b.append(acc)
        proteins_b_xl_positions.append(pos)
        proteins_b_pep_positions.append(pos - xl_pos_b + 1)
    return {
        "xl_pos_a": xl_pos_a,
        "proteins_a": proteins_a,
        "proteins_a_xl_positions": proteins_a_xl_positions,
        "proteins_a_pep_positions": proteins_a_pep_positions,
        "xl_pos_b": xl_pos_b,
        "proteins_b": proteins_b,
        "proteins_b_xl_positions": proteins_b_xl_positions,
        "proteins_b_pep_positions": proteins_b_pep_positions,
    }


def __read_plink_cross_linked_peptides_file(
    file: str | BinaryIO, sep: str = ",", decimal: str = "."
) -> pd.DataFrame:
    r"""Reads a pLink cross-linked peptides file into a pandas DataFrame.

    Reads a pLink cross-linked peptides file into a pandas DataFrame. Reading
    skips all lines starting with ``sep``.

    Parameters
    ----------
    file : str, or BinaryIO
        The name/path of the pLink result file or a file-like object/stream.
    sep : str, default = ","
        Seperator used in the ``.csv`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.

    Returns
    -------
    pd.DataFrame
        The parsed cross-linked peptides.

    Raises
    ------
    RuntimeError
        If the file could not be parsed.
    ValueError
        If the file is not a pLink cross-linked peptides file.

    Notes
    -----
    This function should not be called directly, it is called from ``read_plink()`` and ``detect_plink_filetype()``.
    """
    parsed_lines = None
    if isinstance(file, str):
        with open(file, "r", encoding="utf-8") as f:
            parsed_lines = f.readlines()
            f.close()
    else:
        parsed_lines = file.readlines()
        file.seek(0)
    if parsed_lines is None:
        raise RuntimeError("Something went wrong while parsing the file!")
    lines = list()
    for line in parsed_lines:
        preprocessed_line = str(line).strip()
        if not preprocessed_line.startswith(sep):
            lines.append(preprocessed_line)
    df = pd.read_csv(
        io.StringIO("\n".join(lines)), sep=sep, decimal=decimal, low_memory=False
    )
    colnames = df.columns.values.tolist()
    if not (
        "Peptide_Order" in colnames
        and "Peptide" in colnames
        and "Peptide_Mass" in colnames
        and "Modifications" in colnames
        and "Proteins" in colnames
        and "Protein_Type" in colnames
    ):
        raise ValueError(
            "The provided file seems not to be a pLink cross-linked peptides file!"
        )
    return df


def parse_spectrum_file_from_plink(title: str) -> str:
    r"""Parse the spectrum file name from a spectrum title.

    Parameters
    ----------
    title : str
        The spectrum title.

    Returns
    -------
    str
        The spectrum file name.

    Examples
    --------
    >>> from pyXLMS.parser import parse_spectrum_file_from_plink
    >>> parse_spectrum_file_from_plink(
    ...     "XLpeplib_Beveridge_QEx-HFX_DSS_R1.20588.20588.3.0.dta"
    ... )
    'XLpeplib_Beveridge_QEx-HFX_DSS_R1'
    """
    return str(title).split(".")[0].strip()


def parse_scan_nr_from_plink(title: str) -> int:
    r"""Parse the scan number from a spectrum title.

    Parameters
    ----------
    title : str
        The spectrum title.

    Returns
    -------
    int
        The scan number.

    Examples
    --------
    >>> from pyXLMS.parser import parse_scan_nr_from_plink
    >>> parse_scan_nr_from_plink(
    ...     "XLpeplib_Beveridge_QEx-HFX_DSS_R1.20588.20588.3.0.dta"
    ... )
    20588
    """
    return int(str(title).split(".")[1])


def detect_plink_filetype(
    file: str | BinaryIO, sep: str = ",", decimal: str = "."
) -> Literal["crosslinks", "crosslink-spectrum-matches"]:
    r"""Detects the pLink-related file type of the data.

    Detects whether the input data is a pLink "\*cross-linked_peptides.csv" file or
    a pLink "\*cross-linked_spectra.csv" file.

    Parameters
    ----------
    file : str, or BinaryIO
        The name/path of the pLink result file or a file-like object/stream.
    sep : str, default = ","
        Seperator used in the ``.csv`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.

    Returns
    -------
    str
        Returns "crosslinks" if ``file`` is a "\*cross-linked_peptides.csv" or
        "crosslink-spectrum-matches" if ``file`` is a "\*cross-linked_spectra.csv".

    Raises
    ------
    RuntimeError
        If the file could not be parsed.
    RuntimeError
        If the file does not contain any data.
    ValueError
        If the file does not match any of the supported pLink input files.

    Examples
    --------
    >>> from pyXLMS.parser import detect_plink_filetype
    >>> detect_plink_filetype(
    ...     "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_peptides.csv"
    ... )
    'crosslinks'

    >>> from pyXLMS.parser import detect_plink_filetype
    >>> detect_plink_filetype(
    ...     "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv"
    ... )
    'crosslink-spectrum-matches'
    """
    # parse file
    parsed_lines = None
    if isinstance(file, str):
        with open(file, "r", encoding="utf-8") as f:
            parsed_lines = f.readlines()
            f.close()
    else:
        parsed_lines = file.readlines()
        file.seek(0)
    if parsed_lines is None:
        raise RuntimeError("Something went wrong while parsing the file!")
    lines = list()
    # filter empty lines
    for line in parsed_lines:
        if len(str(line).strip()) > 0:
            lines.append(str(line).strip())
    # checks
    if len(lines) == 0:
        raise RuntimeError("The provided file does not contain any data!")
    if len(lines) == 1:
        df = pd.read_csv(file, sep=sep, decimal=decimal, low_memory=False)
        if not isinstance(file, str):
            file.seek(0)
        colnames = df.columns.values.tolist()
        if (
            "Peptide" in colnames
            and "Modifications" in colnames
            and "Linker" in colnames
            and "Proteins" in colnames
            and "Score" in colnames
            and "Title" in colnames
            and "Charge" in colnames
            and "Evalue" in colnames
            and "Alpha_Evalue" in colnames
            and "Beta_Evalue" in colnames
        ):
            return "crosslink-spectrum-matches"
        else:
            raise ValueError(
                "The provided file seems not to be one of the support pLink input files!"
            )
    if lines[1].startswith(sep):
        _ = __read_plink_cross_linked_peptides_file(file, sep=sep, decimal=decimal)
        return "crosslinks"
    df = pd.read_csv(file, sep=sep, decimal=decimal, low_memory=False)
    if not isinstance(file, str):
        file.seek(0)
    colnames = df.columns.values.tolist()
    if not (
        "Peptide" in colnames
        and "Modifications" in colnames
        and "Linker" in colnames
        and "Proteins" in colnames
        and "Score" in colnames
        and "Title" in colnames
        and "Charge" in colnames
        and "Evalue" in colnames
        and "Alpha_Evalue" in colnames
        and "Beta_Evalue" in colnames
    ):
        raise ValueError(
            "The provided file seems not to be one of the support pLink input files!"
        )
    return "crosslink-spectrum-matches"


def read_plink(
    files: str | List[str] | BinaryIO,
    spectrum_file_parser: Optional[Callable[[str], str]] = None,
    scan_nr_parser: Optional[Callable[[str], int]] = None,
    decoy_prefix: str = "REV_",
    parse_modifications: bool = True,
    modifications: Dict[str, float] = MODIFICATIONS,
    sep: str = ",",
    decimal: str = ".",
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[str, Any]:
    r"""Read a pLink result file.

    Reads a pLink crosslink-spectrum-matches result file "\*cross-linked_spectra.csv"
    in ``.csv`` (comma delimited) format or a pLink crosslinks result file
    "\*cross-linked_peptides.csv" in ``.csv`` (comma delimited) format and
    returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the pLink result file(s) or a file-like object/stream.
    spectrum_file_parser: callable, or None, default = None
        A function that parses the spectrum file name from spectrum titles. If None (default)
        the function ``parse_spectrum_file_from_plink()`` is used.
    scan_nr_parser : callable, or None, default = None
        A function that parses the scan number from spectrum titles. If None (default)
        the function ``parse_scan_nr_from_plink()`` is used.
    decoy_prefix : str, default = "REV\_"
        The prefix that indicates that a protein is from the decoy database.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses.
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
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches.
    TypeError
        If parameter verbose was not set correctly.

    Notes
    -----
    Uses ``Score`` as the score of the crosslink-spectrum-match, all other scores are ``None``.

    Warnings
    --------
    Target and decoy information is derived based off the protein accession and parameter ``decoy_prefix``.
    By default, pLink only reports target matches that are above the desired FDR.

    Examples
    --------
    >>> from pyXLMS.parser import read_plink
    >>> csms = read_plink(
    ...     "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv"
    ... )

    >>> from pyXLMS.parser import read_plink
    >>> crosslinks = read_plink(
    ...     "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_peptides.csv"
    ... )
    """
    ## check input
    _ok = (
        check_input(spectrum_file_parser, "spectrum_file_parser", Callable)
        if spectrum_file_parser is not None
        else True
    )
    _ok = (
        check_input(scan_nr_parser, "scan_nr_parser", Callable)
        if scan_nr_parser is not None
        else True
    )
    _ok = check_input(decoy_prefix, "decoy_prefix", str)
    _ok = check_input(parse_modifications, "parse_modifications", bool)
    _ok = check_input(modifications, "modifications", dict, float)
    _ok = check_input(sep, "sep", str)
    _ok = check_input(decimal, "decimal", str)
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")

    ## set default parsers
    if spectrum_file_parser is None:
        spectrum_file_parser = parse_spectrum_file_from_plink
    if scan_nr_parser is None:
        scan_nr_parser = parse_scan_nr_from_plink

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
        if detect_plink_filetype(input, sep=sep, decimal=decimal) == "crosslinks":
            data = __read_plink_cross_linked_peptides_file(
                input, sep=sep, decimal=decimal
            )
            for i, row in tqdm(
                data.iterrows(), total=data.shape[0], desc="Reading pLink crosslinks..."
            ):
                parsed_positions = __parse_proteins_and_position_from_plink(
                    seq=str(row["Peptide"]).strip(),
                    proteins=str(row["Proteins"]).strip(),
                )
                # create csm
                crosslink = create_crosslink(
                    peptide_a=format_sequence(
                        str(row["Peptide"]).split("-")[0].split("(")[0].strip()
                    ),
                    xl_position_peptide_a=parsed_positions["xl_pos_a"],
                    proteins_a=[
                        protein_a.strip()
                        if protein_a.strip()[: len(decoy_prefix)] != decoy_prefix
                        else protein_a.strip()[len(decoy_prefix) :]
                        for protein_a in parsed_positions["proteins_a"]
                    ],
                    xl_position_proteins_a=parsed_positions["proteins_a_xl_positions"],
                    decoy_a=decoy_prefix in " ".join(parsed_positions["proteins_a"]),
                    peptide_b=format_sequence(
                        str(row["Peptide"]).split("-")[1].split("(")[0].strip()
                    ),
                    xl_position_peptide_b=parsed_positions["xl_pos_b"],
                    proteins_b=[
                        protein_b.strip()
                        if protein_b.strip()[: len(decoy_prefix)] != decoy_prefix
                        else protein_b.strip()[len(decoy_prefix) :]
                        for protein_b in parsed_positions["proteins_b"]
                    ],
                    xl_position_proteins_b=parsed_positions["proteins_b_xl_positions"],
                    decoy_b=decoy_prefix in " ".join(parsed_positions["proteins_b"]),
                    score=None,
                    additional_information={"source": __serialize_pandas_series(row)},
                )
                crosslinks.append(crosslink)
        else:
            data = pd.read_csv(input, sep=sep, decimal=decimal, low_memory=False)
            for i, row in tqdm(
                data.iterrows(), total=data.shape[0], desc="Reading pLink CSMs..."
            ):
                # pre information
                parsed_modifications = (
                    __parse_modifications_from_plink_modifications_str(
                        seq=str(row["Peptide"]).strip(),
                        mod_str=row["Modifications"],  # pyright: ignore [reportArgumentType]
                        crosslinker=str(row["Linker"]).strip(),
                        modifications=modifications,
                        verbose=verbose,
                    )
                    if parse_modifications
                    else None
                )
                parsed_positions = __parse_proteins_and_position_from_plink(
                    seq=str(row["Peptide"]).strip(),
                    proteins=str(row["Proteins"]).strip(),
                )
                # create csm
                csm = create_csm(
                    peptide_a=format_sequence(
                        str(row["Peptide"]).split("-")[0].split("(")[0].strip()
                    ),
                    modifications_a=parsed_modifications[0]
                    if parsed_modifications is not None
                    else None,
                    xl_position_peptide_a=parsed_positions["xl_pos_a"],
                    proteins_a=[
                        protein_a.strip()
                        if protein_a.strip()[: len(decoy_prefix)] != decoy_prefix
                        else protein_a.strip()[len(decoy_prefix) :]
                        for protein_a in parsed_positions["proteins_a"]
                    ],
                    xl_position_proteins_a=parsed_positions["proteins_a_xl_positions"],
                    pep_position_proteins_a=parsed_positions[
                        "proteins_a_pep_positions"
                    ],
                    score_a=None,
                    decoy_a=decoy_prefix in " ".join(parsed_positions["proteins_a"]),
                    peptide_b=format_sequence(
                        str(row["Peptide"]).split("-")[1].split("(")[0].strip()
                    ),
                    modifications_b=parsed_modifications[1]
                    if parsed_modifications is not None
                    else None,
                    xl_position_peptide_b=parsed_positions["xl_pos_b"],
                    proteins_b=[
                        protein_b.strip()
                        if protein_b.strip()[: len(decoy_prefix)] != decoy_prefix
                        else protein_b.strip()[len(decoy_prefix) :]
                        for protein_b in parsed_positions["proteins_b"]
                    ],
                    xl_position_proteins_b=parsed_positions["proteins_b_xl_positions"],
                    pep_position_proteins_b=parsed_positions[
                        "proteins_b_pep_positions"
                    ],
                    score_b=None,
                    decoy_b=decoy_prefix in " ".join(parsed_positions["proteins_b"]),
                    score=float(row["Score"]),
                    spectrum_file=spectrum_file_parser(str(row["Title"]).strip()),
                    scan_nr=scan_nr_parser(str(row["Title"]).strip()),
                    charge=int(row["Charge"]),
                    rt=None,
                    im_cv=None,
                    additional_information={
                        "source": __serialize_pandas_series(row),
                        "Evalue": float(row["Evalue"]),
                        "Alpha_Evalue": float(row["Alpha_Evalue"]),
                        "Beta_Evalue": float(row["Beta_Evalue"]),
                    },
                )
                csms.append(csm)
    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="pLink",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )
