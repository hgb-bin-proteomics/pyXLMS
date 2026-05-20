#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com
__all__ = [
    "read_xi",
    "read_mzid",
    "read_plink",
    "read_scout",
    "read_xlinkx",
    "read_custom",
    "read_merox",
    "read_msannika",
    "read_maxquant",
    "read_maxlynx",
    "detect_xi_filetype",
    "parse_peptide",
    "parse_modifications_from_xi_sequence",
    "parse_scan_nr_from_mzid",
    "parse_scan_nr_from_plink",
    "parse_spectrum_file_from_plink",
    "detect_plink_filetype",
    "detect_scout_filetype",
    "parse_modifications_from_scout_sequence",
    "pyxlms_modification_str_parser",
    "parse_modifications_from_maxquant_sequence",
    "read",
    "read_xinet",
    "read_xiview",
]
# READERS
from ._parser import read
from ._parser._parser_xldbse_xi import read_xi
from ._parser._parser_xldbse_mzid import read_mzid
from ._parser._parser_xldbse_plink import read_plink
from ._parser._parser_xldbse_scout import read_scout
from ._parser._parser_xldbse_xlinkx import read_xlinkx
from ._parser._parser_xldbse_custom import read_custom
from ._parser._parser_xldbse_merox import read_merox
from ._parser._parser_xldbse_msannika import read_msannika
from ._parser._parser_xldbse_maxquant import read_maxquant
from ._parser._parser_xldbse_xinet_xiview import read_xinet
from ._parser._parser_xldbse_maxquant import read_maxlynx
from ._parser._parser_xldbse_xinet_xiview import read_xiview

# UTILITY
from ._parser._parser_xldbse_xi import detect_xi_filetype
from ._parser._parser_xldbse_xi import parse_peptide
from ._parser._parser_xldbse_xi import parse_modifications_from_xi_sequence
from ._parser._parser_xldbse_mzid import parse_scan_nr_from_mzid
from ._parser._parser_xldbse_plink import parse_scan_nr_from_plink
from ._parser._parser_xldbse_plink import parse_spectrum_file_from_plink
from ._parser._parser_xldbse_plink import detect_plink_filetype
from ._parser._parser_xldbse_scout import detect_scout_filetype
from ._parser._parser_xldbse_scout import parse_modifications_from_scout_sequence
from ._parser._parser_xldbse_custom import pyxlms_modification_str_parser
from ._parser._parser_xldbse_maxquant import parse_modifications_from_maxquant_sequence
