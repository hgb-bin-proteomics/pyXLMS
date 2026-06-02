#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
from tqdm import tqdm
from pyteomics import mzid

from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._csm import create_csm
from ..data._parser_result import create_parser_result
from ..constants import CROSSLINKERS
from ._util import format_sequence
from ._util import get_bool_from_value
from ._util import __parse_int

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import List
from typing import Callable

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def parse_scan_nr_from_mzid(spectrum_id: str) -> int:
    r"""Parse the scan number (or spectrum index) from a 'spectrumID' of a mzIdentML file.

    Parameters
    ----------
    spectrum_id : str
        The 'spectrumID' of the mass spectrum from an mzIdentML file read with ``pyteomics``.

    Returns
    -------
    int
        The scan number or spectrum index.

    Notes
    -----
    This function tries to parse the scan number from the 'spectrumID' but does fall back
    to using the spectrum index if the scan number is not available!

    Examples
    --------
    >>> from pyXLMS.parser import parse_scan_nr_from_mzid
    >>> parse_scan_nr_from_mzid("scan=5321")
    5321

    >>> from pyXLMS.parser import parse_scan_nr_from_mzid
    >>> parse_scan_nr_from_mzid("index=1")
    RuntimeWarning: Could not parse scan number from spectrum - using index instead!
    Exception while parsing scan number: list index out of range
    1
    """
    try:
        return __parse_int(str(spectrum_id).split("scan=")[1].split(",")[0])
    except Exception as e:
        warnings.warn(
            RuntimeWarning(
                "Could not parse scan number from spectrum - using index instead!\n"
                f"Exception while parsing scan number: {e}"
            )
        )
    return __parse_int(str(spectrum_id).split("index=")[1].split(",")[0])


def read_mzid(
    files: str | List[str] | BinaryIO,
    scan_nr_parser: Optional[Callable[[str], int]] = None,
    decoy: Optional[bool] = None,
    crosslinkers: Dict[str, float] = CROSSLINKERS,
    verbose: Literal[0, 1, 2] = 1,
) -> ParserResult:
    r"""Read a mzIdentML (mzid) file.

    Reads crosslink-spectrum-matches from a mzIdentML (mzid) file and
    returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the mzIdentML (mzid) file(s) or a file-like object/stream.
    scan_nr_parser : callable, or None, default = None
        A function that parses the scan number from mzid spectrumIDs. If None (default)
        the function ``parse_scan_nr_from_mzid()`` is used.
    decoy : bool, or None, default = None
        Whether the mzid file contains decoy CSMs (``True``) or target CSMs (``False``).
        If None (default) the decoy label is tried to be inferred from the mzIdentML file.
    crosslinkers: dict of str, float, default = ``constants.CROSSLINKERS``
        Mapping of crosslinker names to crosslinker delta masses.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    ParserResult
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches.
    RuntimeError
        If there are warnings while reading the mzIdentML file (only for ``verbose = 2``).
    TypeError
        If parameter verbose was not set correctly.
    TypeError
        If one of the values necessary to create a crosslink-spectrum-match could not be parsed
        correctly.

    Notes
    -----
    This parser only guarantees minimal data because some information might not be available from the mzIdentML file.
    The guaranteed available data is:

    - ``alpha_peptide``
    - ``alpha_peptide_crosslink_position``
    - ``beta_peptide``
    - ``beta_peptide_crosslink_position``
    - ``spectrum_file``
    - ``scan_nr``

    Data that is parsed if available:

    - ``alpha_proteins``
    - ``alpha_proteins_crosslink_positions``
    - ``alpha_proteins_peptide_positions``
    - ``alpha_decoy``
    - ``beta_proteins``
    - ``beta_proteins_crosslink_positions``
    - ``beta_proteins_peptide_positions``
    - ``beta_decoy``

    You can retroactively check which data is available using ``transform.get_available_keys()``!

    Examples
    --------
    >>> from pyXLMS.parser import read_mzid
    >>> csms = read_mzid("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.mzid")
    """
    ## check input
    _ok = (
        check_input(scan_nr_parser, "scan_nr_parser", Callable)
        if scan_nr_parser is not None
        else True
    )
    _ok = check_input(decoy, "decoy", bool) if decoy is not None else True
    _ok = check_input(crosslinkers, "crosslinkers", dict, float)
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")

    ## set default parsers
    if scan_nr_parser is None:
        scan_nr_parser = parse_scan_nr_from_mzid

    ## helper functions
    def check_str(value: str | None) -> str:
        if value is None:
            raise TypeError("Expected str value but None was given!")
        if type(value) is str:
            return value
        raise TypeError(f"Expected str value but {type(value)} was given!")
        return "err"

    def check_int(value: int | None) -> int:
        if value is None:
            raise TypeError("Expected int value but None was given!")
        if type(value) is int:
            return value
        raise TypeError(f"Expected int value but {type(value)} was given!")
        return -1

    def is_xl_mod(modification: Dict[Any, Any]) -> bool:
        if "name" in modification:
            if str(mod["name"]).strip().upper() in crosslinkers:
                return True
        if "crosslink donor" in modification:
            return True
        if "cross-link donor" in modification:
            return True
        if "crosslink acceptor" in modification:
            return True
        if "cross-link acceptor" in modification:
            return True
        if "crosslink receiver" in modification:
            return True
        if "cross-link receiver" in modification:
            return True
        if "search modification id ref" in modification:
            if "crosslink_donor" in "search modification id ref":
                return True
            if "crosslink_acceptor" in "search modification id ref":
                return True
        return False

    def get_proteins_and_positions(
        pep_evidence_list: List[Dict[Any, Any]], pos: int
    ) -> Dict[str, Any]:
        proteins = list()
        xl_position_proteins = list()
        pep_position_proteins = list()
        decoy = None
        for pep_evidence in pep_evidence_list:
            if "start" in pep_evidence:
                try:
                    start = __parse_int(pep_evidence["start"])
                    # positions are 1-indexed in mzIdentML, so if start is smaller than 1
                    # the mzIdentML is incorrect or it's maybe a decoy?
                    if start > 0:
                        xl_position_proteins.append(start + pos - 1)
                        pep_position_proteins.append(start)
                except Exception as _e:
                    pass
            if "accession" in pep_evidence:
                accession = str(pep_evidence["accession"]).strip()
                if len(accession) > 0:
                    proteins.append(accession)
            if "isDecoy" in pep_evidence:
                parsed_decoy = None
                try:
                    parsed_decoy = get_bool_from_value(pep_evidence["isDecoy"])
                except Exception as _e:
                    pass
                if parsed_decoy is not None:
                    if decoy is None:
                        decoy = parsed_decoy
                    else:
                        # if any of the peptides are target, we classify as target
                        decoy = decoy and parsed_decoy
        if len(proteins) > 0 and len(proteins) == len(xl_position_proteins) == len(
            pep_position_proteins
        ):
            return {
                "proteins": proteins,
                "xl": xl_position_proteins,
                "pep": pep_position_proteins,
                "decoy": decoy,
            }
        return {"proteins": None, "xl": None, "pep": None, "decoy": decoy}

    ## data structures
    csms = list()

    ## handle input
    if not isinstance(files, list):
        inputs = [files]
    else:
        inputs = files

    ## process data
    for input in inputs:
        # read all items with pyteomics
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter("always")
            pyteomics_mzid = mzid.MzIdentML(input)
            items = [item for item in pyteomics_mzid]
            pyteomics_mzid.close()
        if verbose > 0 and len(wl) > 0:
            for w in wl:
                warnings.warn(w.message)
        if verbose == 2 and len(wl) > 0:
            raise RuntimeError("Reading mzIdentML file raised warnings!")
        # iterate over all items
        for item in tqdm(
            items, total=len(items), desc="Reading mzIdentML identifications..."
        ):
            # set up empty variables that are needed for a minimal CSM
            csm_id: str | None = None
            scan: int | None = None
            filename: str | None = None
            peptide_a: str | None = None
            pos_a: int | None = None
            peptide_b: str | None = None
            pos_b: int | None = None
            # optional fields
            proteins_a: List[str] | None = None
            xl_position_proteins_a: List[int] | None = None
            pep_position_proteins_a: List[int] | None = None
            decoy_a: bool | None = decoy
            proteins_b: List[str] | None = None
            xl_position_proteins_b: List[int] | None = None
            pep_position_proteins_b: List[int] | None = None
            decoy_b: bool | None = decoy
            # set scan
            if "spectrumID" in item:
                scan = scan_nr_parser(item["spectrumID"])
            # set spectrum file name
            if "location" in item:
                filename = str(item["location"]).strip()
            # check if any identifications for the spectrum
            if "SpectrumIdentificationItem" in item:
                for subitem in item["SpectrumIdentificationItem"]:
                    # we only consider rank 1 CSMs
                    if "rank" in subitem:
                        if __parse_int(subitem["rank"]) > 1:
                            continue
                    # check if item is a CSM
                    if (
                        "cross-link spectrum identification item" in subitem
                        or "crosslink spectrum identification item" in subitem
                    ):
                        parsed_csm_id = (
                            str(subitem["cross-link spectrum identification item"])
                            if "cross-link spectrum identification item" in subitem
                            else str(subitem["crosslink spectrum identification item"])
                        )
                        # if csm_id is not set yet, we parse item as alpha peptide
                        if csm_id is None:
                            csm_id = parsed_csm_id
                            if "PeptideSequence" in subitem:
                                peptide_a = format_sequence(subitem["PeptideSequence"])
                            # we only parse crosslink position from modifications
                            if "Modification" in subitem:
                                for mod in subitem["Modification"]:
                                    if is_xl_mod(mod):
                                        if "location" in mod:
                                            pos_a = __parse_int(mod["location"])
                            if "PeptideEvidenceRef" in subitem:
                                if pos_a is not None:
                                    proteins_and_positions = get_proteins_and_positions(
                                        subitem["PeptideEvidenceRef"], pos_a
                                    )
                                    proteins_a = proteins_and_positions["proteins"]
                                    xl_position_proteins_a = proteins_and_positions[
                                        "xl"
                                    ]
                                    pep_position_proteins_a = proteins_and_positions[
                                        "pep"
                                    ]
                                    decoy_a = (
                                        proteins_and_positions["decoy"]
                                        if decoy is None
                                        else decoy
                                    )
                        # if csm_id is already set, we check if csm_ids of items are equal,
                        # if yes we parse the item as the beta peptide
                        elif csm_id == parsed_csm_id:
                            if "PeptideSequence" in subitem:
                                peptide_b = format_sequence(subitem["PeptideSequence"])
                            if "Modification" in subitem:
                                for mod in subitem["Modification"]:
                                    if is_xl_mod(mod):
                                        if "location" in mod:
                                            pos_b = __parse_int(mod["location"])
                            if "PeptideEvidenceRef" in subitem:
                                if pos_b is not None:
                                    proteins_and_positions = get_proteins_and_positions(
                                        subitem["PeptideEvidenceRef"], pos_b
                                    )
                                    proteins_b = proteins_and_positions["proteins"]
                                    xl_position_proteins_b = proteins_and_positions[
                                        "xl"
                                    ]
                                    pep_position_proteins_b = proteins_and_positions[
                                        "pep"
                                    ]
                                    decoy_b = (
                                        proteins_and_positions["decoy"]
                                        if decoy is None
                                        else decoy
                                    )
            # if and only if all minimal CSM values are parsed, we create a CSM
            if None not in [csm_id, scan, filename, peptide_a, pos_a, peptide_b, pos_b]:
                csm = create_csm(
                    peptide_a=check_str(peptide_a),
                    modifications_a=None,
                    xl_position_peptide_a=check_int(pos_a),
                    proteins_a=proteins_a,
                    xl_position_proteins_a=xl_position_proteins_a,
                    pep_position_proteins_a=pep_position_proteins_a,
                    score_a=None,
                    decoy_a=decoy_a,
                    peptide_b=check_str(peptide_b),
                    modifications_b=None,
                    xl_position_peptide_b=check_int(pos_b),
                    proteins_b=proteins_b,
                    xl_position_proteins_b=xl_position_proteins_b,
                    pep_position_proteins_b=pep_position_proteins_b,
                    score_b=None,
                    decoy_b=decoy_b,
                    score=None,
                    spectrum_file=check_str(filename),
                    scan_nr=check_int(scan),
                    charge=None,
                    rt=None,
                    im_cv=None,
                )
                csms.append(csm)
    ## check results
    if len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="mzIdentML",
        csms=csms,
        crosslinks=None,
    )
