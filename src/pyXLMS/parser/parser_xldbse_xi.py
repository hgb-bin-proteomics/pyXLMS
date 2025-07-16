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
from ..constants import XI_MODIFICATION_MAPPING
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


def detect_xi_filetype(
    data: pd.DataFrame,
) -> Literal["xisearch", "xifdr_csms", "xifdr_crosslinks"]:
    r"""Detects the xi-related source (application) of the data.

    Detects whether the input data is originating from xiSearch or xiFDR, and if xiFDR which type of data is
    being read (crosslink-spectrum-matches or crosslinks).

    Parameters
    ----------
    data : pd.DataFrame
        The input data originating from xiSearch or xiFDR.

    Returns
    -------
    str
        "xisearch" if a xiSearch result file was read, "xifdr_csms" if CSMs from xiFDR were read,
        "xifdr_crosslinks" if crosslinks from xiFDR were read.

    Raises
    ------
    ValueError
        If the data source could not be determined.

    Examples
    --------
    >>> from pyXLMS.parser import detect_xi_filetype
    >>> import pandas as pd
    >>> df1 = pd.read_csv("data/xi/r1_Xi1.7.6.7.csv")
    >>> detect_xi_filetype(df1)
    'xisearch'

    >>> from pyXLMS.parser import detect_xi_filetype
    >>> import pandas as pd
    >>> df2 = pd.read_csv("data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv")
    >>> detect_xi_filetype(df2)
    'xifdr_csms'

    >>> from pyXLMS.parser import detect_xi_filetype
    >>> import pandas as pd
    >>> df3 = pd.read_csv("data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv")
    >>> detect_xi_filetype(df3)
    'xifdr_crosslinks'
    """
    ## check input
    _ok = check_input(data, "data", pd.DataFrame)

    col_names = data.columns.values.tolist()
    if "AllScore" in col_names:
        return "xisearch"
    if "LinkPos1" in col_names:
        return "xifdr_csms"
    if "ToSite" in col_names:
        return "xifdr_crosslinks"

    raise ValueError(
        "Could not infer data source, are you sure you read a xi result file?"
    )

    return "err"


def parse_peptide(sequence: str, term_char: str = ".") -> str:
    r"""Parses the peptide sequence from a sequence string including flanking amino acids.

    Parses the peptide sequence from a sequence string including flanking amino acids, for example ``"K.KKMoxKLS.S"``.
    The returned peptide sequence for this example would be ``"KKMoxKLS"``.

    Parameters
    ----------
    sequence : str
        The sequence string containing the peptide sequence and flanking amino acids.
    term_char : str (single character), default = "."
        The character used to denote N-terminal and C-terminal.

    Returns
    -------
    str
        The parsed peptide sequence without flanking amino acids.

    Raises
    ------
    RuntimeError
        If (one of) the peptide sequence(s) could not be parsed.

    Examples
    --------
    >>> from pyXLMS.parser import parse_peptide
    >>> parse_peptide("K.KKMoxKLS.S")
    'KKMoxKLS'

    >>> from pyXLMS.parser import parse_peptide
    >>> parse_peptide("-.CcmCcmPSR.T")
    'CcmCcmPSR'

    >>> from pyXLMS.parser import parse_peptide
    >>> parse_peptide("CCPSR")
    'CCPSR'
    """
    ## check input
    _ok = check_input(sequence, "sequence", str)

    # PEPTIDE
    if term_char not in sequence and len(sequence.strip()) > 1:
        return sequence.strip()
    if term_char in sequence:
        parts = [part.strip() for part in sequence.split(term_char)]
        # K.PEPTPIDE.P.EP <- wrong format
        if len(parts) > 3:
            raise RuntimeError(f"Could not parse peptide from sequence {sequence}!")
        # K.PEPTIDE.R
        if len(parts) == 3 and len(parts[1]) > 1:
            return parts[1]
        if len(parts) == 2:
            # PEPTIDE.R
            if len(parts[0]) > 1 and len(parts[1]) == 1:
                return parts[0]
            # K.PEPTIDE
            if len(parts[1]) > 1 and len(parts[0]) == 1:
                return parts[1]
    # if none of these cases match, raise error
    raise RuntimeError(f"Could not parse peptide from sequence {sequence}!")
    return "err"


def parse_modifications_from_xi_sequence(sequence: str) -> Dict[int, str]:
    r"""Parses all post-translational-modifications from a peptide sequence as reported by xiFDR.

    Parses all post-translational-modifications from a peptide sequence as reported by xiFDR. This assumes
    that amino acids are given in upper case letters and post-translational-modifications in lower case letters.
    The parsed modifications are returned as a dictionary that maps their position in the sequence (1-based) to
    their xiFDR annotation (``SYMBOLEXT``), for example ``"cm"`` or ``"ox"``.

    Parameters
    ----------
    sequence : str
        The peptide sequence as given by xiFDR.

    Returns
    -------
    dict of int, str
        Dictionary that maps modifications (values) to their respective positions in the peptide sequence (1-based)
        (keys). The modifications are given in xiFDR annotation style (``SYMBOLEXT``) which is the lower letter
        modification code, for example ``"cm"`` for carbamidomethylation.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed.

    Examples
    --------
    >>> from pyXLMS.parser import parse_modifications_from_xi_sequence
    >>> seq1 = "KIECcmFDSVEISGVEDR"
    >>> parse_modifications_from_xi_sequence(seq1)
    {4: 'cm'}

    >>> from pyXLMS.parser import parse_modifications_from_xi_sequence
    >>> seq2 = "KIECcmFDSVEMoxISGVEDR"
    >>> parse_modifications_from_xi_sequence(seq2)
    {4: 'cm', 10: 'ox'}

    >>> from pyXLMS.parser import parse_modifications_from_xi_sequence
    >>> seq3 = "KIECcmFDSVEISGVEDRMox"
    >>> parse_modifications_from_xi_sequence(seq3)
    {4: 'cm', 17: 'ox'}

    >>> from pyXLMS.parser import parse_modifications_from_xi_sequence
    >>> seq4 = "CcmKIECcmFDSVEISGVEDRMox"
    >>> parse_modifications_from_xi_sequence(seq4)
    {1: 'cm', 5: 'cm', 18: 'ox'}
    """
    ## check input
    _ok = check_input(sequence, "sequence", str)

    modifications = dict()
    pos = 0
    current_mod = ""
    for i, aa in enumerate(str(sequence).strip()):
        if aa.isupper():
            pos += 1
            current_mod = ""
        else:
            current_mod += aa
            if (i + 1 >= len(sequence)) or (sequence[i + 1].isupper()):
                if pos in modifications:
                    raise RuntimeError(
                        f"Modification at position {pos} already exists!"
                    )
                modifications[pos] = current_mod
    return modifications


def __parse_xisearch_modifications(
    row: pd.Series,
    alpha: bool,
    modifications: Dict[str, Tuple[str, float]] = XI_MODIFICATION_MAPPING,
    ignore_errors: bool = False,
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[int, Tuple[str, float]]:
    r"""Returns the corresponding modifications object for a crosslink-spectrum-match from xiSearch.

    Parameters
    ----------
    row : pandas.Series
        One row/crosslink-spectrum-match of the xiSearch result file.
    alpha : bool
        Whether to parse modifications from the alpha peptide or - if ``False`` - from the beta peptide.
    modifications : dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"cm"``) to their modifications (e.g. ``("Carbamidomethyl", 57.021464)``).
    ignore_errors : bool, default = False
        If modifications that are not given in parameter 'modifications' should raise an error or not. By default an error is
        raised if an unknown modification is encountered. If ``True`` modifications that are unknown are encoded with the xi
        shortcode (``SYMBOLEXT``) and ``float("nan")`` modification mass.
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
        If the parsed modifications and positions are not of the same length.
    RuntimeError
        If multiple modifications on the same residue are parsed (only for ``verbose = 2``).
    KeyError
        If an unknown modification is encountered.

    Notes
    -----
    This function should not be called directly, it is called from ``__read_xisearch()``.
    """

    # EXAMPLE VALUES
    # Modifications2            Mox;Mox
    # ModificationPositions2    5;7
    # helper function that changes ``SYMBOL`` to ``SYMBOLEXT``
    def preprocess_mod(mod: str) -> str:
        return "".join([c for c in mod if not c.isupper()]).strip()

    crosslinker = str(row["Crosslinker"]).strip()
    crosslinker_mass = float(row["CrosslinkerMass"])
    parsed_modifications = dict()
    # parse from Modifications
    if alpha:
        parsed_modifications[int(row["Link1"])] = (crosslinker, crosslinker_mass)
        if not pd.isna(row["Modifications1"]):  # pyright: ignore [reportGeneralTypeIssues]
            if ";" in str(row["Modifications1"]):
                mods = [
                    preprocess_mod(mod) for mod in str(row["Modifications1"]).split(";")
                ]
                positions = [
                    int(pos) for pos in str(row["ModificationPositions1"]).split(";")
                ]
                if len(mods) != len(positions):
                    err_str = "Parsed modifications and their positions are not of the same length!\n"
                    err_str += f"Parsed modifications: {row['Modifications1']}; Parsed positions: {row['ModificationPositions1']}\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    raise RuntimeError(err_str)
                for i in range(len(mods)):
                    if positions[i] in parsed_modifications:
                        err_str = (
                            f"Modification at position {positions[i]} already exists!\n"
                        )
                        err_str += (
                            f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        )
                        if verbose == 1:
                            warnings.warn(RuntimeWarning(err_str))
                        elif verbose == 2:
                            raise RuntimeError(err_str)
                        try:
                            t1 = parsed_modifications[positions[i]][0] + (
                                "," + modifications[mods[i]][0]
                            )
                            t2 = (
                                parsed_modifications[positions[i]][1]
                                + modifications[mods[i]][1]
                            )
                            parsed_modifications[positions[i]] = (t1, t2)
                        except KeyError:
                            if ignore_errors:
                                t1 = (
                                    parsed_modifications[positions[i]][0]
                                    + ","
                                    + mods[i]
                                )
                                t2 = parsed_modifications[positions[i]][1] + float(
                                    "nan"
                                )
                                parsed_modifications[positions[i]] = (t1, t2)
                            else:
                                err_str = f"Key {mods[i]} not found in parameter 'modifications'. Are you missing a modification?\n"
                                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                                raise KeyError(err_str)
                    else:
                        try:
                            parsed_modifications[positions[i]] = (
                                modifications[mods[i]][0],
                                modifications[mods[i]][1],
                            )
                        except KeyError:
                            if ignore_errors:
                                parsed_modifications[positions[i]] = (
                                    mods[i],
                                    float("nan"),
                                )
                            else:
                                err_str = f"Key {mods[i]} not found in parameter 'modifications'. Are you missing a modification?\n"
                                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                                raise KeyError(err_str)
            else:
                mod = preprocess_mod(str(row["Modifications1"]))
                pos = int(row["ModificationPositions1"])
                if pos in parsed_modifications:
                    err_str = f"Modification at position {pos} already exists!\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    if verbose == 1:
                        warnings.warn(RuntimeWarning(err_str))
                    elif verbose == 2:
                        raise RuntimeError(err_str)
                    try:
                        t1 = parsed_modifications[pos][0] + "," + modifications[mod][0]
                        t2 = parsed_modifications[pos][1] + modifications[mod][1]
                        parsed_modifications[pos] = (t1, t2)
                    except KeyError:
                        if ignore_errors:
                            t1 = parsed_modifications[pos][0] + "," + mod
                            t2 = parsed_modifications[pos][1] + float("nan")
                            parsed_modifications[pos] = (t1, t2)
                        else:
                            err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                            err_str += (
                                f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                            )
                            raise KeyError(err_str)
                else:
                    try:
                        parsed_modifications[pos] = (
                            modifications[mod][0],
                            modifications[mod][1],
                        )
                    except KeyError:
                        if ignore_errors:
                            parsed_modifications[pos] = (
                                mod,
                                float("nan"),
                            )
                        else:
                            err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                            err_str += (
                                f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                            )
                            raise KeyError(err_str)
    else:
        parsed_modifications[int(row["Link2"])] = (crosslinker, crosslinker_mass)
        if not pd.isna(row["Modifications2"]):  # pyright: ignore [reportGeneralTypeIssues]
            if ";" in str(row["Modifications2"]):
                mods = [
                    preprocess_mod(mod) for mod in str(row["Modifications2"]).split(";")
                ]
                positions = [
                    int(pos) for pos in str(row["ModificationPositions2"]).split(";")
                ]
                if len(mods) != len(positions):
                    err_str = "Parsed modifications and their positions are not of the same length!\n"
                    err_str += f"Parsed modifications: {row['Modifications2']}; Parsed positions: {row['ModificationPositions2']}\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    raise RuntimeError(err_str)
                for i in range(len(mods)):
                    if positions[i] in parsed_modifications:
                        err_str = (
                            f"Modification at position {positions[i]} already exists!\n"
                        )
                        err_str += (
                            f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        )
                        if verbose == 1:
                            warnings.warn(RuntimeWarning(err_str))
                        elif verbose == 2:
                            raise RuntimeError(err_str)
                        try:
                            t1 = parsed_modifications[positions[i]][0] + (
                                "," + modifications[mods[i]][0]
                            )
                            t2 = (
                                parsed_modifications[positions[i]][1]
                                + modifications[mods[i]][1]
                            )
                            parsed_modifications[positions[i]] = (t1, t2)
                        except KeyError:
                            if ignore_errors:
                                t1 = parsed_modifications[positions[i]][0] + (
                                    "," + mods[i]
                                )
                                t2 = parsed_modifications[positions[i]][1] + float(
                                    "nan"
                                )
                                parsed_modifications[positions[i]] = (t1, t2)
                            else:
                                err_str = f"Key {mods[i]} not found in parameter 'modifications'. Are you missing a modification?\n"
                                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                                raise KeyError(err_str)
                    else:
                        try:
                            parsed_modifications[positions[i]] = (
                                modifications[mods[i]][0],
                                modifications[mods[i]][1],
                            )
                        except KeyError:
                            if ignore_errors:
                                parsed_modifications[positions[i]] = (
                                    mods[i],
                                    float("nan"),
                                )
                            else:
                                err_str = f"Key {mods[i]} not found in parameter 'modifications'. Are you missing a modification?\n"
                                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                                raise KeyError(err_str)
            else:
                mod = preprocess_mod(str(row["Modifications2"]))
                pos = int(row["ModificationPositions2"])
                if pos in parsed_modifications:
                    err_str = f"Modification at position {pos} already exists!\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    if verbose == 1:
                        warnings.warn(RuntimeWarning(err_str))
                    elif verbose == 2:
                        raise RuntimeError(err_str)
                    try:
                        t1 = parsed_modifications[pos][0] + "," + modifications[mod][0]
                        t2 = parsed_modifications[pos][1] + modifications[mod][1]
                        parsed_modifications[pos] = (t1, t2)
                    except KeyError:
                        if ignore_errors:
                            t1 = parsed_modifications[pos][0] + "," + mod
                            t2 = parsed_modifications[pos][1] + float("nan")
                            parsed_modifications[pos] = (t1, t2)
                        else:
                            err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                            err_str += (
                                f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                            )
                            raise KeyError(err_str)
                else:
                    try:
                        parsed_modifications[pos] = (
                            modifications[mod][0],
                            modifications[mod][1],
                        )
                    except KeyError:
                        if ignore_errors:
                            parsed_modifications[pos] = (
                                mod,
                                float("nan"),
                            )
                        else:
                            err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                            err_str += (
                                f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                            )
                            raise KeyError(err_str)
    # parse from sequence (because fixed modifcations are not reported in Modifications)
    if alpha:
        modified_sequence = parse_peptide(str(row["Peptide1"]).strip())
        mods_from_sequence = parse_modifications_from_xi_sequence(modified_sequence)
        for pos, mod in mods_from_sequence.items():
            if pos in parsed_modifications:
                err_str = f"Modification at position {pos} already exists!\n"
                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                if verbose == 1:
                    warnings.warn(RuntimeWarning(err_str))
                elif verbose == 2:
                    raise RuntimeError(err_str)
                mod_mapped = None
                try:
                    mod_mapped = modifications[mod]
                except KeyError:
                    if ignore_errors:
                        mod_mapped = (mod, float("nan"))
                    else:
                        err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                        err_str += (
                            f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        )
                        raise KeyError(err_str)
                if mod_mapped is not None and isinstance(mod_mapped, tuple):
                    if mod_mapped[0] not in parsed_modifications[pos][0]:
                        parsed_modifications[pos] = (
                            parsed_modifications[pos][0] + "," + mod_mapped[0],
                            parsed_modifications[pos][1] + mod_mapped[1],
                        )
            else:
                mod_mapped = None
                try:
                    mod_mapped = modifications[mod]
                except KeyError:
                    if ignore_errors:
                        mod_mapped = (mod, float("nan"))
                    else:
                        err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                        err_str += (
                            f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        )
                        raise KeyError(err_str)
                if mod_mapped is not None and isinstance(mod_mapped, tuple):
                    parsed_modifications[pos] = mod_mapped
    else:
        modified_sequence = parse_peptide(str(row["Peptide2"]).strip())
        mods_from_sequence = parse_modifications_from_xi_sequence(modified_sequence)
        for pos, mod in mods_from_sequence.items():
            if pos in parsed_modifications:
                err_str = f"Modification at position {pos} already exists!\n"
                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                if verbose == 1:
                    warnings.warn(RuntimeWarning(err_str))
                elif verbose == 2:
                    raise RuntimeError(err_str)
                mod_mapped = None
                try:
                    mod_mapped = modifications[mod]
                except KeyError:
                    if ignore_errors:
                        mod_mapped = (mod, float("nan"))
                    else:
                        err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                        err_str += (
                            f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        )
                        raise KeyError(err_str)
                if mod_mapped is not None and isinstance(mod_mapped, tuple):
                    if mod_mapped[0] not in parsed_modifications[pos][0]:
                        parsed_modifications[pos] = (
                            parsed_modifications[pos][0] + "," + mod_mapped[0],
                            parsed_modifications[pos][1] + mod_mapped[1],
                        )
            else:
                mod_mapped = None
                try:
                    mod_mapped = modifications[mod]
                except KeyError:
                    if ignore_errors:
                        mod_mapped = (mod, float("nan"))
                    else:
                        err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                        err_str += (
                            f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        )
                        raise KeyError(err_str)
                if mod_mapped is not None and isinstance(mod_mapped, tuple):
                    parsed_modifications[pos] = mod_mapped
    return parsed_modifications


def __read_xisearch(
    data: pd.DataFrame,
    decoy_prefix: str,
    parse_modifications: bool,
    modifications: Dict[str, Tuple[str, float]],
    ignore_errors: bool,
    verbose: Literal[0, 1, 2],
) -> List[Dict[str, Any]]:
    r"""Reads a xiSearch pandas dataframe and returns a list of crosslink-spectrum-matches.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe of a xiSearch result ``.csv`` file read with pandas.
    decoy_prefix : str
        The prefix that indicates that a protein is from the decoy database.
    parse_modifications : bool
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications : dict of str, tuple
        Mapping of xi sequence elements (e.g. ``"cm"``) to their modifications (e.g. ``("Carbamidomethyl", 57.021464)``).
    ignore_errors : bool
        If modifications that are not given in parameter 'modifications' should raise an error or not. By default an error is
        raised if an unknown modification is encountered. If ``True`` modifications that are unknown are encoded with the xi
        shortcode (``SYMBOLEXT``) and ``float("nan")`` modification mass.
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
    This function should not be called directly, it is called from ``read_xi()``.
    """
    # remove monolinks
    xl = data.dropna(axis=0, subset=["BasePeptide2"])
    # create csms list
    csms = list()
    # create csms
    for i, row in tqdm(
        xl.iterrows(), total=xl.shape[0], desc="Reading xiSearch CSMs..."
    ):
        csm = create_csm(
            peptide_a=format_sequence(str(row["BasePeptide1"])),
            modifications_a=__parse_xisearch_modifications(
                row, True, modifications, ignore_errors, verbose
            )
            if parse_modifications
            else None,
            xl_position_peptide_a=int(row["Link1"]),
            proteins_a=[
                p.strip()
                if p.strip()[: len(decoy_prefix)] != decoy_prefix
                else p.strip()[len(decoy_prefix) :]
                for p in str(row["Protein1"]).split(";")
            ],
            xl_position_proteins_a=[
                int(float(p)) for p in str(row["ProteinLink1"]).split(";")
            ],
            pep_position_proteins_a=[
                int(float(p)) for p in str(row["Start1"]).split(";")
            ],
            score_a=float(row["Pep1Score"]),
            decoy_a=get_bool_from_value(int(row["Protein1decoy"])),
            peptide_b=format_sequence(str(row["BasePeptide2"])),
            modifications_b=__parse_xisearch_modifications(
                row, False, modifications, ignore_errors, verbose
            )
            if parse_modifications
            else None,
            xl_position_peptide_b=int(row["Link2"]),
            proteins_b=[
                p.strip()
                if p.strip()[: len(decoy_prefix)] != decoy_prefix
                else p.strip()[len(decoy_prefix) :]
                for p in str(row["Protein2"]).split(";")
            ],
            xl_position_proteins_b=[
                int(float(p)) for p in str(row["ProteinLink2"]).split(";")
            ],
            pep_position_proteins_b=[
                int(float(p)) for p in str(row["Start2"]).split(";")
            ],
            score_b=float(row["Pep2Score"]),
            decoy_b=get_bool_from_value(int(row["Protein2decoy"])),
            score=float(row["match score"]),
            spectrum_file=str(row["peakListFileName"]).strip(),
            scan_nr=int(row["Scan"]),
            charge=int(row["PrecoursorCharge"]),
            rt=None,
            im_cv=None,
            additional_information={
                "spectrum quality score": float(row["spectrum quality score"]),
            },
        )
        csms.append(csm)
    return csms


def __parse_xifdr_modifications(
    row: pd.Series,
    alpha: bool,
    modifications: Dict[str, Tuple[str, float]] = XI_MODIFICATION_MAPPING,
    ignore_errors: bool = False,
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[int, Tuple[str, float]]:
    r"""Returns the corresponding modifications object for a crosslink-spectrum-match from xiFDR.

    Parameters
    ----------
    row : pandas.Series
        One row/crosslink-spectrum-match of the xiFDR CSM result file.
    alpha : bool
        Whether to parse modifications from the alpha peptide or - if ``False`` - from the beta peptide.
    modifications : dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"cm"``) to their modifications (e.g. ``("Carbamidomethyl", 57.021464)``).
    ignore_errors : bool, default = False
        If modifications that are not given in parameter 'modifications' should raise an error or not. By default an error is
        raised if an unknown modification is encountered. If ``True`` modifications that are unknown are encoded with the xi
        shortcode (``SYMBOLEXT``) and ``float("nan")`` modification mass.
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

    Notes
    -----
    This function should not be called directly, it is called from ``__read_xifdr_csms()``.
    """
    crosslinker = str(row["Crosslinker"]).strip()
    crosslinker_mass = float(row["CrosslinkerModMass"])
    parsed_modifications = dict()
    if alpha:
        parsed_modifications[int(row["LinkPos1"])] = (crosslinker, crosslinker_mass)
        for pos, mod in parse_modifications_from_xi_sequence(
            str(row["PepSeq1"]).strip()
        ).items():
            if pos in parsed_modifications:
                err_str = f"Modification at position {pos} already exists!\n"
                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                if verbose == 1:
                    warnings.warn(RuntimeWarning(err_str))
                elif verbose == 2:
                    raise RuntimeError(err_str)
                try:
                    t1 = parsed_modifications[pos][0] + "," + modifications[mod][0]
                    t2 = parsed_modifications[pos][1] + modifications[mod][1]
                    parsed_modifications[pos] = (t1, t2)
                except KeyError:
                    if ignore_errors:
                        t1 = parsed_modifications[pos][0] + "," + mod
                        t2 = parsed_modifications[pos][1] + float("nan")
                        parsed_modifications[pos] = (t1, t2)
                    else:
                        err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                        err_str += (
                            f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        )
                        raise KeyError(err_str)
            try:
                parsed_modifications[pos] = (
                    modifications[mod][0],
                    modifications[mod][1],
                )
            except KeyError:
                if ignore_errors:
                    parsed_modifications[pos] = (mod, float("nan"))
                else:
                    err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    raise KeyError(err_str)
    else:
        parsed_modifications[int(row["LinkPos2"])] = (crosslinker, crosslinker_mass)
        for pos, mod in parse_modifications_from_xi_sequence(
            str(row["PepSeq2"]).strip()
        ).items():
            if pos in parsed_modifications:
                err_str = f"Modification at position {pos} already exists!\n"
                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                if verbose == 1:
                    warnings.warn(RuntimeWarning(err_str))
                elif verbose == 2:
                    raise RuntimeError(err_str)
                try:
                    t1 = parsed_modifications[pos][0] + "," + modifications[mod][0]
                    t2 = parsed_modifications[pos][1] + modifications[mod][1]
                    parsed_modifications[pos] = (t1, t2)
                except KeyError:
                    if ignore_errors:
                        t1 = parsed_modifications[pos][0] + "," + mod
                        t2 = parsed_modifications[pos][1] + float("nan")
                        parsed_modifications[pos] = (t1, t2)
                    else:
                        err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                        err_str += (
                            f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        )
                        raise KeyError(err_str)
            try:
                parsed_modifications[pos] = (
                    modifications[mod][0],
                    modifications[mod][1],
                )
            except KeyError:
                if ignore_errors:
                    parsed_modifications[pos] = (mod, float("nan"))
                else:
                    err_str = f"Key {mod} not found in parameter 'modifications'. Are you missing a modification?\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    raise KeyError(err_str)
    return parsed_modifications


def __read_xifdr_csms(
    data: pd.DataFrame,
    decoy_prefix: str,
    parse_modifications: bool,
    modifications: Dict[str, Tuple[str, float]],
    ignore_errors: bool,
    verbose: Literal[0, 1, 2],
) -> List[Dict[str, Any]]:
    r"""Reads a xiFDR CSM pandas dataframe and returns a list of crosslink-spectrum-matches.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe of a xiFDR CSM result ``.csv`` file read with pandas.
    decoy_prefix : str
        The prefix that indicates that a protein is from the decoy database.
    parse_modifications : bool
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications : dict of str, tuple
        Mapping of xi sequence elements (e.g. ``"cm"``) to their modifications (e.g. ``("Carbamidomethyl", 57.021464)``).
    ignore_errors : bool
        If modifications that are not given in parameter 'modifications' should raise an error or not. By default an error is
        raised if an unknown modification is encountered. If ``True`` modifications that are unknown are encoded with the xi
        shortcode (``SYMBOLEXT``) and ``float("nan")`` modification mass.
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
    This function should not be called directly, it is called from ``read_xi()``.
    """
    # create csms list
    csms = list()
    # create csms
    for i, row in tqdm(
        data.iterrows(), total=data.shape[0], desc="Reading xiFDR CSMs..."
    ):
        csm = create_csm(
            peptide_a=format_sequence(str(row["PepSeq1"])),
            modifications_a=__parse_xifdr_modifications(
                row, True, modifications, ignore_errors, verbose
            )
            if parse_modifications
            else None,
            xl_position_peptide_a=int(row["LinkPos1"]),
            proteins_a=[
                p.strip()
                if p.strip()[: len(decoy_prefix)] != decoy_prefix
                else p.strip()[len(decoy_prefix) :]
                for p in str(row["Protein1"]).split(";")
            ],
            xl_position_proteins_a=[
                int(p) for p in str(row["ProteinLinkPos1"]).split(";")
            ],
            pep_position_proteins_a=[int(p) for p in str(row["PepPos1"]).split(";")],
            score_a=None,
            decoy_a=get_bool_from_value(row["Decoy1"]),
            peptide_b=format_sequence(str(row["PepSeq2"])),
            modifications_b=__parse_xifdr_modifications(
                row, False, modifications, ignore_errors, verbose
            )
            if parse_modifications
            else None,
            xl_position_peptide_b=int(row["LinkPos2"]),
            proteins_b=[
                p.strip()
                if p.strip()[: len(decoy_prefix)] != decoy_prefix
                else p.strip()[len(decoy_prefix) :]
                for p in str(row["Protein2"]).split(";")
            ],
            xl_position_proteins_b=[
                int(p) for p in str(row["ProteinLinkPos2"]).split(";")
            ],
            pep_position_proteins_b=[int(p) for p in str(row["PepPos2"]).split(";")],
            score_b=None,
            decoy_b=get_bool_from_value(row["Decoy2"]),
            score=float(row["Score"]),
            spectrum_file=str(row["PeakListFileName"]).strip(),
            scan_nr=int(row["scan"]),
            charge=int(row["exp charge"]),
            rt=None,
            im_cv=None,
            additional_information=None,
        )
        csms.append(csm)
    return csms


def __read_xifdr_crosslinks(
    data: pd.DataFrame, decoy_prefix: str
) -> List[Dict[str, Any]]:
    r"""Reads a xiFDR Links pandas dataframe and returns a list of crosslinks.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe of a xiFDR Links result ``.csv`` file read with pandas.
    decoy_prefix : str
        The prefix that indicates that a protein is from the decoy database.

    Returns
    -------
    list of dict
        The read crosslinks.

    Notes
    -----
    This function should not be called directly, it is called from ``read_xi()``.
    """
    # create crosslink list
    crosslinks = list()
    # create crosslinks
    for i, row in tqdm(
        data.iterrows(), total=data.shape[0], desc="Reading xiFDR crosslinks..."
    ):
        psmid = str(row["PSMIDs"]).split(";")[0]
        s1 = psmid.split("P1_")[1].split(" ")[0]
        p1 = parse_peptide(s1)
        s2 = psmid.split("P2_")[1].split(" ")[0]
        p2 = parse_peptide(s2)
        pos1 = int(psmid.split("P2_")[1].split(" ")[1])
        pos2 = int(psmid.split("P2_")[1].split(" ")[2])
        crosslink = create_crosslink(
            peptide_a=format_sequence(p1),
            xl_position_peptide_a=pos1,
            proteins_a=[
                p.strip()
                if p.strip()[: len(decoy_prefix)] != decoy_prefix
                else p.strip()[len(decoy_prefix) :]
                for p in str(row["Protein1"]).split(";")
            ],
            xl_position_proteins_a=[int(p) for p in str(row["fromSite"]).split(";")],
            decoy_a=get_bool_from_value(row["Decoy1"]),
            peptide_b=format_sequence(p2),
            xl_position_peptide_b=pos2,
            proteins_b=[
                p.strip()
                if p.strip()[: len(decoy_prefix)] != decoy_prefix
                else p.strip()[len(decoy_prefix) :]
                for p in str(row["Protein2"]).split(";")
            ],
            xl_position_proteins_b=[int(p) for p in str(row["ToSite"]).split(";")],
            decoy_b=get_bool_from_value(row["Decoy2"]),
            score=float(row["Score"]),
            additional_information=None,
        )
        crosslinks.append(crosslink)
    return crosslinks


def read_xi(
    files: str | List[str] | BinaryIO,
    decoy_prefix: Optional[str] = "auto",
    parse_modifications: bool = True,
    modifications: Dict[str, Tuple[str, float]] = XI_MODIFICATION_MAPPING,
    sep: str = ",",
    decimal: str = ".",
    ignore_errors: bool = False,
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[str, Any]:
    r"""Read a xiSearch/xiFDR result file.

    Reads a xiSearch crosslink-spectrum-matches result file or a xiFDR crosslink-spectrum-matches
    result file or crosslink result file in ``.csv`` format and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the xiSearch/xiFDR result file(s) or a file-like object/stream.
    decoy_prefix : str, or None, default = "auto"
        The prefix that indicates that a protein is from the decoy database.
        If "auto" or None it will use the default for each xi file type.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications : dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"cm"``) to their modifications (e.g. ``("Carbamidomethyl", 57.021464)``).
        This corresponds to the ``SYMBOLEXT`` field, or the ``SYMBOL`` field minus the amino acid in the xiSearch config.
    sep : str, default = ","
        Seperator used in the ``.csv`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.
    ignore_errors : bool, default = False
        If modifications that are not given in parameter 'modifications' should raise an error or not. By default an error is
        raised if an unknown modification is encountered. If ``True`` modifications that are unknown are encoded with the xi
        shortcode (``SYMBOLEXT``) and ``float("nan")`` modification mass.
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
        If the file(s) contain no crosslinks or crosslink-spectrum-matches.
    TypeError
        If parameter verbose was not set correctly.

    Examples
    --------
    >>> from pyXLMS.parser import read_xi
    >>> csms_from_xiSearch = read_xi("data/xi/r1_Xi1.7.6.7.csv")

    >>> from pyXLMS.parser import read_xi
    >>> csms_from_xiFDR = read_xi("data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv")

    >>> from pyXLMS.parser import read_xi
    >>> crosslinks_from_xiFDR = read_xi("data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv")
    """
    ## check input
    _ok = (
        check_input(decoy_prefix, "decoy_prefix", str)
        if decoy_prefix is not None
        else True
    )
    _ok = check_input(parse_modifications, "parse_modifications", bool)
    _ok = check_input(modifications, "modifications", dict, tuple)
    _ok = check_input(sep, "sep", str)
    _ok = check_input(decimal, "decimal", str)
    _ok = check_input(ignore_errors, "ignore_errors", bool)
    _ok = check_input(verbose, "verbose", int)
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
        xi_file_type = detect_xi_filetype(data)
        ## set decoy prefix
        if decoy_prefix is None or decoy_prefix == "auto":
            decoy_prefix = "REV_" if xi_file_type == "xisearch" else "decoy:"
        ## process data
        if xi_file_type == "xifdr_csms":
            csms += __read_xifdr_csms(
                data,
                decoy_prefix,
                parse_modifications,
                modifications,
                ignore_errors,
                verbose,
            )
        elif xi_file_type == "xifdr_crosslinks":
            crosslinks += __read_xifdr_crosslinks(data, decoy_prefix)
        else:
            csms += __read_xisearch(
                data,
                decoy_prefix,
                parse_modifications,
                modifications,
                ignore_errors,
                verbose,
            )

    ## check results
    if len(crosslinks) + len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches or crosslinks were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="xiSearch/xiFDR",
        csms=csms if len(csms) > 0 else None,
        crosslinks=crosslinks if len(crosslinks) > 0 else None,
    )
