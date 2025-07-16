#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

__all__ = [
    "read_xi",
    "read_mzid",
    "read_plink",
    "read_scout",
    "read_xlinkx",
    "read_custom",
    "read_msannika",
    "read_maxquant",
    "read_maxlynx",
    "detect_xi_filetype",
    "parse_peptide",
    "parse_modifications_from_xi_sequence",
    "parse_scan_nr_from_mzid",
    "parse_scan_nr_from_plink",
    "parse_spectrum_file_from_plink",
    "detect_scout_filetype",
    "parse_modifications_from_scout_sequence",
    "pyxlms_modification_str_parser",
    "parse_modifications_from_maxquant_sequence",
    "read",
]

# READERS
from .parser_xldbse_xi import read_xi
from .parser_xldbse_mzid import read_mzid
from .parser_xldbse_plink import read_plink
from .parser_xldbse_scout import read_scout
from .parser_xldbse_xlinkx import read_xlinkx
from .parser_xldbse_custom import read_custom
from .parser_xldbse_msannika import read_msannika
from .parser_xldbse_maxquant import read_maxquant
from .parser_xldbse_maxquant import read_maxlynx

# UTILITY
from .parser_xldbse_xi import detect_xi_filetype
from .parser_xldbse_xi import parse_peptide
from .parser_xldbse_xi import parse_modifications_from_xi_sequence
from .parser_xldbse_mzid import parse_scan_nr_from_mzid
from .parser_xldbse_plink import parse_scan_nr_from_plink
from .parser_xldbse_plink import parse_spectrum_file_from_plink
from .parser_xldbse_scout import detect_scout_filetype
from .parser_xldbse_scout import parse_modifications_from_scout_sequence
from .parser_xldbse_custom import pyxlms_modification_str_parser
from .parser_xldbse_maxquant import parse_modifications_from_maxquant_sequence

from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def read(
    files: str | List[str] | BinaryIO,
    engine: Literal[
        "Custom",
        "MaxQuant",
        "MaxLynx",
        "MS Annika",
        "mzIdentML",
        "pLink",
        "Scout",
        "xiSearch/xiFDR",
        "XlinkX",
    ],
    crosslinker: str,
    parse_modifications: bool = True,
    ignore_errors: bool = False,
    verbose: Literal[0, 1, 2] = 1,
    **kwargs,
) -> Dict[str, Any]:
    r"""Read a crosslink result file.

    Reads a crosslink or crosslink-spectrum-match result file from any of the supported crosslink search engines or formats.
    Currently supports results files from MaxLynx/MaxQuant, MS Annika, pLink 2 and pLink 3, Scout, xiSearch and xiFDR,
    XlinkX, and the mzIdentML format. Additionally supports parsing from custom ``.csv`` files in pyXLMS format, see more
    about the custom format in ``parser.read_custom()`` and in here:
    `docs <https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md>`_.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the result file(s) or a file-like object/stream.
    engine : "Custom", "MaxQuant", "MaxLynx", "MS Annika", "mzIdentML", "pLink", "Scout", "xiSearch/xiFDR", or "XlinkX"
        Crosslink search engine or format of the result file.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter for every parser. Defaults are selected
        for every parser if 'modifications' is not passed via ``**kwargs``.
    ignore_errors : bool, default = False
        Ignore errors when mapping modifications. Used in ``parser.read_xi()`` and ``parser.read_xlinkx()``.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.
    **kwargs
        Any additional parameters will be passed to the specific parsers.

    Returns
    -------
    dict
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    ValueError
        If the value entered for parameter ``engine`` is not supported.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> csms_from_xiSearch = read("data/xi/r1_Xi1.7.6.7.csv", engine="xiSearch/xiFDR", crosslinker="DSS")

    >>> from pyXLMS.parser import read
    >>> csms_from_MaxQuant = read("data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS")
    """
    supported = [
        "Custom",
        "MaxQuant",
        "MaxLynx",
        "MS Annika",
        "mzIdentML",
        "pLink",
        "Scout",
        "xiSearch/xiFDR",
        "XlinkX",
    ]
    ff = engine.lower().strip()

    if ff in ["custom", "pyxlms"]:
        return read_custom(files, parse_modifications=parse_modifications, **kwargs)
    if ff in ["maxquant", "max quant"]:
        return read_maxquant(
            files,
            crosslinker=crosslinker,
            parse_modifications=parse_modifications,
            **kwargs,
        )
    if ff in ["maxlynx", "max lynx"]:
        return read_maxlynx(
            files,
            crosslinker=crosslinker,
            parse_modifications=parse_modifications,
            **kwargs,
        )
    if ff in ["ms annika", "msannika"]:
        return read_msannika(
            files, parse_modifications=parse_modifications, verbose=verbose, **kwargs
        )
    if ff in ["mzidentml", "mzid"]:
        return read_mzid(files, verbose=verbose, **kwargs)
    if ff in ["plink", "plink2", "plink3", "plink 2", "plink 3"]:
        return read_plink(
            files, parse_modifications=parse_modifications, verbose=verbose, **kwargs
        )
    if ff in ["scout"]:
        return read_scout(
            files,
            crosslinker=crosslinker,
            parse_modifications=parse_modifications,
            verbose=verbose,
            **kwargs,
        )
    if ff in ["xisearch/xifdr", "xisearch", "xifdr", "xi search", "xi fdr", "xi"]:
        return read_xi(
            files,
            parse_modifications=parse_modifications,
            ignore_errors=ignore_errors,
            verbose=verbose,
            **kwargs,
        )
    if ff in ["xlinkx", "x link x"]:
        return read_xlinkx(
            files,
            parse_modifications=parse_modifications,
            ignore_errors=ignore_errors,
            verbose=verbose,
            **kwargs,
        )

    err_str = (
        f"{engine} is not a supported crosslink search engine or format! Valid options are:\n"
        + ", ".join(supported)
    )
    raise ValueError(err_str)

    return {"err": err_str}
