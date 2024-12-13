#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS import parser as p
    
    assert p.format_sequence("PEP[K]TIDE") == "PEPKTIDE"
    assert p.format_sequence("PEPKdssoTIDE") == "PEPKTIDE"
    assert p.format_sequence("peptide", remove_lower = False) == "PEPTIDE"
    
