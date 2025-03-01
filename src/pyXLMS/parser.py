#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from .parser_msannika import read_msannika

## TODO
def read(file: str):
    if file:
        return read_msannika(file)
    return
