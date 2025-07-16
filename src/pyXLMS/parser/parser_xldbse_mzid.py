#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
from tqdm import tqdm
from pyteomics import mzid

from ..data import check_input
from ..data import create_csm
from ..data import create_parser_result
from ..constants import CROSSLINKERS
from .util import format_sequence

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
    r"""Parse the scan number from a 'spectrumID' of a mzIdentML file.

    Parameters
    ----------
    title : str
        The 'spectrumID' of the mass spectrum from an mzIdentML file read with ``pyteomics``.

    Returns
    -------
    int
        The scan number.

    Examples
    --------
    >>> from pyXLMS.parser import parse_scan_nr_from_mzid
    >>> parse_scan_nr_from_mzid("scan=5321")
    5321
    """
    return int(str(spectrum_id).split("scan=")[1].split(",")[0])


def read_mzid(
    files: str | List[str] | BinaryIO,
    scan_nr_parser: Optional[Callable[[str], int]] = None,
    decoy: Optional[bool] = None,
    crosslinkers: Dict[str, float] = CROSSLINKERS,
    verbose: Literal[0, 1, 2] = 1,
) -> Dict[str, Any]:
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
    crosslinkers: dict of str, float, default = ``constants.CROSSLINKERS``
        Mapping of crosslinker names to crosslinker delta masses.
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
    RuntimeError
        If parser is used with ``verbose = 2``.
    RuntimeError
        If there are warnings while reading the mzIdentML file (only for ``verbose = 2``).
    TypeError
        If parameter verbose was not set correctly.
    TypeError
        If one of the values necessary to create a crosslink-spectrum-match could not be parsed
        correctly.

    Notes
    -----
    This parser is experimental, as I don't know if the mzIdentML structure is consistent accross different
    crosslink search engines. This parser was tested with mzIdentML files from MS Annika and XlinkX.

    Warnings
    --------
    This parser only parses minimal data because most information is not available from the mzIdentML file.
    The available data is:

    - ``alpha_peptide``
    - ``alpha_peptide_crosslink_position``
    - ``beta_peptide``
    - ``beta_peptide_crosslink_position``
    - ``spectrum_file``
    - ``scan_nr``

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

    ## warning message
    if verbose == 1:
        warnings.warn(
            UserWarning(
                "Please be aware that mzIdentML parsing is currently an experimental feature!\n"
                "Please check the documentation for parser.read_mzid for more information!"
            )
        )
    if verbose == 2:
        raise RuntimeError(
            "Please be aware that mzIdentML parsing is currently an experimental feature!\n"
            "Please check the documentation for parser.read_mzid for more information!"
        )

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
        for item in tqdm(items):
            # set up empty variables that are needed for a minimal CSM
            csm_id = None
            scan = None
            filename = None
            peptide_a = None
            pos_a = None
            peptide_b = None
            pos_b = None
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
                        if int(subitem["rank"]) > 1:
                            continue
                    # check if item is a CSM
                    if "cross-link spectrum identification item" in subitem:
                        # if csm_id is not set yet, we parse item as alpha peptide
                        if csm_id is None:
                            csm_id = int(
                                float(
                                    subitem["cross-link spectrum identification item"]
                                )
                            )
                            if "PeptideSequence" in subitem:
                                peptide_a = format_sequence(subitem["PeptideSequence"])
                            # we only parse crosslink position from modifications
                            if "Modification" in subitem:
                                for mod in subitem["Modification"]:
                                    if "name" in mod:
                                        if (
                                            str(mod["name"]).strip().upper()
                                            in crosslinkers
                                        ):
                                            if "location" in mod:
                                                pos_a = int(mod["location"])
                        # if csm_id is already set, we check if csm_ids of items are equal,
                        # if yes we parse the item as the beta peptide
                        elif csm_id == int(
                            float(subitem["cross-link spectrum identification item"])
                        ):
                            if "PeptideSequence" in subitem:
                                peptide_b = format_sequence(subitem["PeptideSequence"])
                            if "Modification" in subitem:
                                for mod in subitem["Modification"]:
                                    if "name" in mod:
                                        if (
                                            str(mod["name"]).strip().upper()
                                            in crosslinkers
                                        ):
                                            if "location" in mod:
                                                pos_b = int(mod["location"])
            # if and only if all minimal CSM values are parsed, we create a CSM
            if None not in [csm_id, scan, filename, peptide_a, pos_a, peptide_b, pos_b]:
                csm = create_csm(
                    peptide_a=check_str(peptide_a),
                    modifications_a=None,
                    xl_position_peptide_a=check_int(pos_a),
                    proteins_a=None,
                    xl_position_proteins_a=None,
                    pep_position_proteins_a=None,
                    score_a=None,
                    decoy_a=decoy,
                    peptide_b=check_str(peptide_b),
                    modifications_b=None,
                    xl_position_peptide_b=check_int(pos_b),
                    proteins_b=None,
                    xl_position_proteins_b=None,
                    pep_position_proteins_b=None,
                    score_b=None,
                    decoy_b=decoy,
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
