#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


# READERS
from .parser_xi import read_xi
from .parser_msannika import read_msannika
from .parser_xldbse_maxquant import read_maxquant
from .parser_xldbse_maxquant import read_maxlynx

# UTILITY
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
    return
