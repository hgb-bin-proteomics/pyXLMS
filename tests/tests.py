#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

def test1():
    from pyXLMS import data
    crosslink = parser.create_crosslink("PEPTIDE", 1, "EDITPEP", 2)
    assert crosslink["alpha_peptide"] == "EDITPEP"
    assert crosslink["alpha_crosslink_position"] == 2
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_crosslink_position"] == 1

def test2():
    from pyXLMS import data
    crosslink = parser.create_crosslink("PEPTIDE", 2, "PEPTIDE", 1)
    assert crosslink["alpha_peptide"] == "PEPTIDE"
    assert crosslink["alpha_crosslink_position"] == 1
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_crosslink_position"] == 2
