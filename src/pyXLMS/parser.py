#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


# READERS
from .parser_xldbse_xi import read_xi
from .parser_xldbse_plink import read_plink
from .parser_xldbse_scout import read_scout
from .parser_xldbse_xlinkx import read_xlinkx
from .parser_xldbse_custom import read_custom
from .parser_xldbse_msannika import read_msannika
from .parser_xldbse_maxquant import read_maxquant
from .parser_xldbse_maxquant import read_maxlynx

# UTILITY
from .parser_xldbse_xi import detect_xi_filetype  # noqa: F401
from .parser_xldbse_xi import parse_peptide  # noqa: F401
from .parser_xldbse_xi import parse_modifications_from_xi_sequence  # noqa: F401
from .parser_xldbse_plink import parse_scan_nr_from_plink  # noqa: F401
from .parser_xldbse_plink import parse_spectrum_file_from_plink  # noqa: F401
from .parser_xldbse_scout import detect_scout_filetype  # noqa: F401
from .parser_xldbse_scout import parse_modifications_from_scout_sequence  # noqa: F401
from .parser_xldbse_custom import pyxlms_modification_str_parser  # noqa: F401
from .parser_xldbse_maxquant import parse_modifications_from_maxquant_sequence  # noqa: F401


## TODO
def read(file: str, dbse: str):
    if dbse == "MS Annika":
        return read_msannika(file)
    if dbse == "Xi":
        return read_xi(file)
    if dbse == "MaxQuant":
        return read_maxquant(file, "DSSO")
    if dbse == "MaxLynx":
        return read_maxlynx(file, "DSSO")
    if dbse == "pLink":
        return read_plink(file)
    if dbse == "XlinkX":
        return read_xlinkx(file)
    if dbse == "Scout":
        return read_scout(file, "DSSO")
    if dbse == "Custom":
        return read_custom(file)
    return
