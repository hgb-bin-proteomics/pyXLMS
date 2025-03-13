#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from .parser_msannika import read_msannika
from .parser_xi import read_xi


## TODO
def read(file: str, dbse: str):
    if dbse == "MS Annika":
        return read_msannika(file)
    if dbse == "Xi":
        return read_xi(file)
    return
