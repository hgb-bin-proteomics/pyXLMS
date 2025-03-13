#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .data import create_crosslink
from .data import create_csm
from .data import create_parser_result
from .constants import XI_MODIFICATION_MAPPING
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

def detect_xi_filetype(data: pd.DataFrame) -> Literal["xisearch", "xifdr_csms", "xifdr_crosslinks"]:
    """Detects the source application of the data.

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
    >>> from pyXLMS.parser_xi import detect_xi_filetype
    >>> import pandas as pd
    >>> df1 = pd.read_csv("data/xi/r1_Xi1.7.6.7.csv")
    >>> detect_xi_filetype(df1)
    'xisearch'

    >>> from pyXLMS.parser_xi import detect_xi_filetype
    >>> import pandas as pd
    >>> df2 = pd.read_csv("data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv")
    >>> detect_xi_filetype(df2)
    'xifdr_csms'

    >>> from pyXLMS.parser_xi import detect_xi_filetype
    >>> import pandas as pd
    >>> df3 = pd.read_csv("data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv")
    >>> detect_xi_filetype(df3)
    'xifdr_crosslinks'
    """
    col_names = data.columns.values.tolist()
    if "AllScore" in col_names:
        return "xisearch"
    if "LinkPos1" in col_names:
        return "xifdr_csms"
    if "ToSite" in col_names:
        return "xifdr_crosslinks"

    raise ValueError("Could not infer data source, are you sure you read a xi result file?")

    return "err"

def parse_modifications_from_xi_sequence(sequence: str) -> Dict[int, str]:
    """Parses all post-translational-modifications from a peptide sequence as reported by xiFDR.

    Parses all post-translational-modifications from a peptide sequence as reported by xiFDR. This assumes
    that amino acids are given in upper case letters and post-translational-modifications in lower case letters.
    The parsed modifications are returned as a dictionary that maps their position in the sequence (1-based) to
    their xiFDR annotation, for example ``"Ccm"`` or ``"Mox"``.

    Parameters
    ----------
    sequence : str
        The peptide sequence as given by xiFDR.

    Returns
    -------
    dict of int, str
        Dictionary that maps modifications (values) to their respective positions in the peptide sequence (1-based)
        (keys). The modifications are given in xiFDR annotation style which is the amino acid followed by the lower
        letter modification code, for example ``"Ccm"`` for carbamidomethylation of cysteine.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed.

    Examples
    --------
    >>> from pyXLMS.parser_xi import parse_modifications_from_xi_sequence
    >>> seq1 = "KIECcmFDSVEISGVEDR"
    >>> parse_modifications_from_xi_sequence(seq1)
    {4: 'Ccm'}

    >>> from pyXLMS.parser_xi import parse_modifications_from_xi_sequence
    >>> seq2 = "KIECcmFDSVEMoxISGVEDR"
    >>> parse_modifications_from_xi_sequence(seq2)
    {4: 'Ccm', 10: 'Mox'}

    >>> from pyXLMS.parser_xi import parse_modifications_from_xi_sequence
    >>> seq3 = "KIECcmFDSVEISGVEDRMox"
    >>> parse_modifications_from_xi_sequence(seq3)
    {4: 'Ccm', 17: 'Mox'}

    >>> from pyXLMS.parser_xi import parse_modifications_from_xi_sequence
    >>> seq4 = "CcmKIECcmFDSVEISGVEDRMox"
    >>> parse_modifications_from_xi_sequence(seq4)
    {1: 'Ccm', 5: 'Ccm', 18: 'Mox'}
    """
    modifications = dict()
    pos = 0
    current_mod = ""
    for i, aa in enumerate(str(sequence).strip()):
        if aa.isupper():
            pos += 1
            current_mod = aa
        else:
            current_mod += aa
            if i + 1 >= len(sequence):
                if pos in modifications:
                    raise RuntimeError(f"Modification at position {pos} already exists!")
                modifications[pos] = current_mod
            elif sequence[i + 1].isupper():
                if pos in modifications:
                    raise RuntimeError(f"Modification at position {pos} already exists!")
                modifications[pos] = current_mod
    return modifications

def __parse_xisearch_modifications(
    row: pd.Series,
    alpha: bool,
    modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING
) -> Dict[int, Tuple[str, float]]:
    """Returns the corresponding modifications object for a crosslink-spectrum-match from xiSearch.

    Parameters
    ----------
    row : pandas.Series
        One row/crosslink-spectrum-match of the xiSearch result file.
    alpha : bool
        Whether to parse modifications from the alpha peptide or - if ``False`` - from the beta peptide.
    modifications: dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"Ccm"``) to their modifications (e.g. ``("C", "Carbamidomethyl", 57.021464)``).

    Returns
    -------
    dict of int, tuple
        The ``pyXLMS`` specific modifications object, a dictionary that maps positions to their corresponding modifications and their
        monoisotopic masses.

    Raises
    ------
    RuntimeError
        If the parsed modifications and positions are not of the same length.
        If multiple modifications on the same residue are parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``__read_xisearch()``.
    """
    # EXAMPLE VALUES
    # Modifications2            Mox;Mox
    # ModificationPositions2    5;7
    crosslinker = str(row["Crosslinker"]).strip()
    crosslinker_mass = float(row["CrosslinkerMass"])
    modifications = dict()
    if alpha:
        modifications[int(row["Link1"])] = (crosslinker, crosslinker_mass)
        if not pd.isna(row["Modifications1"]):
            if ";" in str(row["Modifications1"]):
                mods = [mod.strip() for mod in str(row["Modifications1"]).split(";")]
                positions = [int(pos) for pos in str(row["ModificationPositions1"]).split(";")]
                if len(mods) != len(positions):
                    err_str = "Parsed modifications and their positions are not of the same length!\n"
                    err_str += f"Parsed modifications: {row['Modifications1']}; Parsed positions: {row['ModificationPositions1']}\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    raise RuntimeError(err_str)
                for i in range(len(mods)):
                    if positions[i] in modifications:
                        err_str = f"Modification at position {positions[i]} already exists!\n"
                        err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        raise RuntimeError(err_str)
                    modifications[positions[i]] = (XI_MODIFICATION_MAPPING[mods[i]][1], XI_MODIFICATION_MAPPING[mods[i]][2])
            else:
                mod = str(row["Modifications1"]).strip()
                pos = int(row["ModificationPositions1"])
                if pos in modifications:
                    err_str = f"Modification at position {pos} already exists!\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    raise RuntimeError(err_str)
                modifications[pos] = (XI_MODIFICATION_MAPPING[mod][1], XI_MODIFICATION_MAPPING[mod][2])
    else:
        modifications[int(row["Link2"])] = (crosslinker, crosslinker_mass)
        if not pd.isna(row["Modifications2"]):
            if ";" in str(row["Modifications2"]):
                mods = [mod.strip() for mod in str(row["Modifications2"]).split(";")]
                positions = [int(pos) for pos in str(row["ModificationPositions2"]).split(";")]
                if len(mods) != len(positions):
                    err_str = "Parsed modifications and their positions are not of the same length!\n"
                    err_str += f"Parsed modifications: {row['Modifications2']}; Parsed positions: {row['ModificationPositions2']}\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    raise RuntimeError(err_str)
                for i in range(len(mods)):
                    if positions[i] in modifications:
                        err_str = f"Modification at position {positions[i]} already exists!\n"
                        err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                        raise RuntimeError(err_str)
                    modifications[positions[i]] = (XI_MODIFICATION_MAPPING[mods[i]][1], XI_MODIFICATION_MAPPING[mods[i]][2])
            else:
                mod = str(row["Modifications2"]).strip()
                pos = int(row["ModificationPositions2"])
                if pos in modifications:
                    err_str = f"Modification at position {pos} already exists!\n"
                    err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                    raise RuntimeError(err_str)
                modifications[pos] = (XI_MODIFICATION_MAPPING[mod][1], XI_MODIFICATION_MAPPING[mod][2])
    return modifications

def __read_xisearch(
    data: pd.DataFrame,
    modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING
) -> List[Dict[str, Any]]:
    """Reads a xiSearch pandas dataframe and returns a list of crosslink-spectrum-matches.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe of a xiSearch result ``.csv`` file read with pandas.
    modifications: dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"Ccm"``) to their modifications (e.g. ``("C", "Carbamidomethyl", 57.021464)``).

    Returns
    -------
    list of dict
        The read crosslink-spectrum-matches.

    Notes
    -----
    This function should not be called directly, it is called from ``read_xi()``.
    """
    # remove monolinks
    xl = data.dropna(axis = 0, subset = "BasePeptide2")
    # create csms list
    csms = list()
    # create csms
    for i, row in xl.iterrows():
        csm = create_csm(
            peptide_a = format_sequence(str(row["BasePeptide1"])),
            modifications_a = __parse_xisearch_modifications(row, True, modifications),
            xl_position_peptide_a = int(row["Link1"]),
            proteins_a = [p.strip() if p.strip()[:4] != "REV_" else p.strip()[4:] for p in str(row["Protein1"]).split(";")],
            xl_position_proteins_a = [int(float(p)) for p in str(row["ProteinLink1"]).split(";")],
            pep_position_proteins_a = [int(float(p)) for p in str(row["Start1"]).split(";")],
            score_a = float(row["Pep1Score"]),
            decoy_a = get_bool_from_value(int(row["Protein1decoy"])),
            peptide_b = format_sequence(str(row["BasePeptide2"])),
            modifications_b = __parse_xisearch_modifications(row, False, modifications),
            xl_position_peptide_b = int(row["Link2"]),
            proteins_b = [p.strip() if p.strip()[:4] != "REV_" else p.strip()[4:] for p in str(row["Protein2"]).split(";")],
            xl_position_proteins_b = [int(float(p)) for p in str(row["ProteinLink2"]).split(";")],
            pep_position_proteins_b = [int(float(p)) for p in str(row["Start2"]).split(";")],
            score_b = float(row["Pep2Score"]),
            decoy_b = get_bool_from_value(int(row["Protein2decoy"])),
            score = float(row["match score"]),
            spectrum_file = str(row["peakListFileName"]).strip(),
            scan_nr = int(row["Scan"]),
            charge = int(row["PrecoursorCharge"]),
            rt = None,
            im_cv = None,
            additional_information = {
                "spectrum quality score": float(row["spectrum quality score"]),
            }
        )
        csms.append(csm)
    return csms

def __parse_xifdr_modifications(
    row: pd.Series,
    alpha: bool,
    modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING
) -> Dict[int, Tuple[str, float]]:
    """Returns the corresponding modifications object for a crosslink-spectrum-match from xiFDR.

    Parameters
    ----------
    row : pandas.Series
        One row/crosslink-spectrum-match of the xiFDR CSM result file.
    alpha : bool
        Whether to parse modifications from the alpha peptide or - if ``False`` - from the beta peptide.
    modifications: dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"Ccm"``) to their modifications (e.g. ``("C", "Carbamidomethyl", 57.021464)``).

    Returns
    -------
    dict of int, tuple
        The ``pyXLMS`` specific modifications object, a dictionary that maps positions to their corresponding modifications and their
        monoisotopic masses.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``__read_xifdr_csms()``.
    """
    crosslinker = str(row["Crosslinker"]).strip()
    crosslinker_mass = float(row["CrosslinkerModMass"])
    parsed_modifications = dict()
    if alpha:
        parsed_modifications[int(row["LinkPos1"])] = (crosslinker, crosslinker_mass)
        for pos, mod in parse_modifications_from_xi_sequence(str(row["PepSeq1"]).strip()).items():
            if pos in parsed_modifications:
                err_str = f"Modification at position {pos} already exists!\n"
                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                raise RuntimeError(err_str)
            parsed_modifications[pos] = (modifications[mod][1], modifications[mod][2])
    else:
        parsed_modifications[int(row["LinkPos2"])] = (crosslinker, crosslinker_mass)
        for pos, mod in parse_modifications_from_xi_sequence(str(row["PepSeq2"]).strip()).items():
            if pos in parsed_modifications:
                err_str = f"Modification at position {pos} already exists!\n"
                err_str += f"CSM ScanId: {row['ScanId']}; CSM Scan: {row['Scan']}"
                raise RuntimeError(err_str)
            parsed_modifications[pos] = (modifications[mod][1], modifications[mod][2])
    return parsed_modifications

def __read_xifdr_csms(
    data: pd.DataFrame,
    modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING
) -> List[Dict[str, Any]]:
    """Reads a xiFDR CSM pandas dataframe and returns a list of crosslink-spectrum-matches.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe of a xiFDR CSM result ``.csv`` file read with pandas.
    modifications: dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"Ccm"``) to their modifications (e.g. ``("C", "Carbamidomethyl", 57.021464)``).

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
    for i, row in data.iterrows():
        csm = create_csm(
            peptide_a = format_sequence(str(row["PepSeq1"])),
            modifications_a = __parse_xifdr_modifications(row, True, modifications),
            xl_position_peptide_a = int(row["LinkPos1"]),
            proteins_a = [p.strip() if p.strip()[:6] != "decoy:" else p.strip()[6:] for p in str(row["Protein1"]).split(";")],
            xl_position_proteins_a = [int(p) for p in str(row["ProteinLinkPos1"]).split(";")],
            pep_position_proteins_a = [int(p) for p in str(row["PepPos1"]).split(";")],
            score_a = None,
            decoy_a = get_bool_from_value(row["Decoy1"]),
            peptide_b = format_sequence(str(row["PepSeq2"])),
            modifications_b = __parse_xifdr_modifications(row, False, modifications),
            xl_position_peptide_b = int(row["LinkPos2"]),
            proteins_b = [p.strip() if p.strip()[:6] != "decoy:" else p.strip()[6:] for p in str(row["Protein2"]).split(";")],
            xl_position_proteins_b = [int(p) for p in str(row["ProteinLinkPos2"]).split(";")],
            pep_position_proteins_b = [int(p) for p in str(row["PepPos2"]).split(";")],
            score_b = None,
            decoy_b = get_bool_from_value(row["Decoy2"]),
            score = float(row["Score"]),
            spectrum_file = str(row["PeakListFileName"]).strip(),
            scan_nr = int(row["scan"]),
            charge = int(row["exp charge"]),
            rt = None,
            im_cv = None,
            additional_information = None,
        )
        csms.append(csm)
    return csms

def __read_xifdr_crosslinks(data: pd.DataFrame) -> List[Dict[str, Any]]:
    """Reads a xiFDR Links pandas dataframe and returns a list of crosslinks.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe of a xiFDR Links result ``.csv`` file read with pandas.

    Returns
    -------
    list of dict
        The read crosslinks.

    Raises
    ------
    RuntimeError
        If (one of) the peptide sequence(s) could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_xi()``.
    """
    # helper function
    def parse_peptide(sequence: str) -> str:
        # PEPTIDE
        if "." not in sequence and len(sequence.strip()) > 1:
            return sequence.strip()
        if "." in sequence:
            parts = [part.strip() for part in sequence.split(".")]
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
    # create crosslink list
    crosslinks = list()
    # create crosslinks
    for i, row in data.iterrows():
        psmid = str(row["PSMIDs"]).split(";")[0]
        s1 = psmid.split("P1_")[1].split(" ")[0]
        p1 = parse_peptide(s1)
        s2 = psmid.split("P2_")[1].split(" ")[0]
        p2 = parse_peptide(s2)
        pos1 = int(psmid.split("P2_")[1].split(" ")[1])
        pos2 = int(psmid.split("P2_")[1].split(" ")[2])
        crosslink = create_crosslink(
            peptide_a = format_sequence(p1),
            xl_position_peptide_a = pos1,
            proteins_a = [p.strip() if p.strip()[:6] != "decoy:" else p.strip()[6:] for p in str(row["Protein1"]).split(";")],
            xl_position_proteins_a = [int(p) for p in str(row["fromSite"]).split(";")],
            decoy_a = get_bool_from_value(row["Decoy1"]),
            peptide_b = format_sequence(p2),
            xl_position_peptide_b = pos2,
            proteins_b = [p.strip() if p.strip()[:6] != "decoy:" else p.strip()[6:] for p in str(row["Protein2"]).split(";")],
            xl_position_proteins_b = [int(p) for p in str(row["ToSite"]).split(";")],
            decoy_b = get_bool_from_value(row["Decoy2"]),
            score = float(row["Score"]),
            additional_information = None,
        )
        crosslinks.append(crosslink)
    return crosslinks

def read_xi(
    files: str | List[str] | BinaryIO,
    modifications: Dict[str, Tuple[Any]] = XI_MODIFICATION_MAPPING
) -> Dict[str, Any]:
    """Read a xiSearch/xiFDR result file.

    Reads a xiSearch crosslink-spectrum-matches result file or a xiFDR crosslink-spectrum-matches
    result file or crosslink result file in ``.csv`` format and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the xiSearch/xiFDR result file(s) or a file-like object/stream.
    modifications: dict of str, tuple, default = ``constants.XI_MODIFICATION_MAPPING``
        Mapping of xi sequence elements (e.g. ``"Ccm"``) to their modifications (e.g. ``("C", "Carbamidomethyl", 57.021464)``).

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) contain no crosslinks or crosslink-spectrum-matches.

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
    _ok = check_input(modifications, "modifications", dict, tuple)

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
        data = pd.read_csv(input)
        ## detect input file type
        xi_file_type = detect_xi_filetype(data)
        ## process data
        if xi_file_type == "xifdr_csms":
            csms += __read_xifdr_csms(data, modifications)
        elif xi_file_type == "xifdr_crosslinks":
            crosslinks += __read_xifdr_crosslinks(data, modifications)
        else:
            csms += __read_xisearch(data, modifications)

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
