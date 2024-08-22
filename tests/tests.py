#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

def test1():
    from pyXLMS import data
    x = 1
    assert data.check_input(x, "x", int)

def test2():
    from pyXLMS import data
    x = [1, 2, 3]
    assert data.check_input(x, "x", list)

def test3():
    from pyXLMS import data
    x = [1, 2, 3]
    assert data.check_input(x, "x", list, int)

def test4():
    from pyXLMS import data
    x = 1
    with pytest.raises(TypeError, match = f"x must be {str}!"):
        i = data.check_input(x, "x", str)

def test5():
    from pyXLMS import data
    x = [1, 2, 3]
    with pytest.raises(TypeError, match = f"List values of x must be {str}"):
        i = data.check_input(x, "x", list, str)

def test6():
    from pyXLMS import data
    crosslink = data.create_crosslink("PEPTIDE", 1, ["PROTEIN"], [1],
                                      "EDITPEP", 3, ["NIETORP", "PROTEIN"], [5, 2],
                                      170.3)
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["alpha_peptide"] == "EDITPEP"
    assert crosslink["alpha_peptide_crosslink_position"] == 3
    assert len(crosslink["alpha_proteins"]) == 2
    assert crosslink["alpha_proteins"][0] == "NIETORP"
    assert len(crosslink["alpha_proteins_crosslink_positions"]) == 2
    assert crosslink["alpha_proteins_crosslink_positions"][0] == 5
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 1
    assert len(crosslink["beta_proteins"]) == 1
    assert crosslink["beta_proteins"][0] == "PROTEIN"
    assert len(crosslink["beta_proteins_crosslink_positions"]) == 1
    assert crosslink["beta_proteins_crosslink_positions"][0] == 1
    assert crosslink["score"] >= 170.25 and crosslink["score"] <= 170.35

def test7():
    from pyXLMS import data
    crosslink = data.create_crosslink("PEPTIDE", 3, ["PROTEIN"], [3],
                                      "PEPTIDE", 1, ["PROTEIN"], [1],
                                      170.3)
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["alpha_peptide"] == "PEPTIDE"
    assert crosslink["alpha_peptide_crosslink_position"] == 1
    assert len(crosslink["alpha_proteins"]) == 1
    assert crosslink["alpha_proteins"][0] == "PROTEIN"
    assert len(crosslink["alpha_proteins_crosslink_positions"]) == 1
    assert crosslink["alpha_proteins_crosslink_positions"][0] == 1
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 3
    assert len(crosslink["beta_proteins"]) == 1
    assert crosslink["beta_proteins"][0] == "PROTEIN"
    assert len(crosslink["beta_proteins_crosslink_positions"]) == 1
    assert crosslink["beta_proteins_crosslink_positions"][0] == 3
    assert crosslink["score"] >= 170.25 and crosslink["score"] <= 170.35

def test8():
    from pyXLMS import data
    with pytest.raises(TypeError, match = f"xl_position_peptide_a must be {int}!"):
        crosslink = data.create_crosslink("PEPTIDE", "3", ["PROTEIN"], [3],
                                          "PEPTIDE", 1, ["PROTEIN"], [1],
                                          170.3)

def test9():
    from pyXLMS import data
    with pytest.raises(TypeError, match = f"List values of xl_position_proteins_a must be {int}"):
        crosslink = data.create_crosslink("PEPTIDE", 3, ["PROTEIN"], ["3"],
                                          "PEPTIDE", 1, ["PROTEIN"], [1],
                                          170.3)

def test10():
    from pyXLMS import data
    with pytest.raises(ValueError, match = "Crosslink position has to be given for every protein! Length of proteins_a and xl_position_proteins_a has to match!"):
        crosslink = data.create_crosslink("PEPTIDE", 3, ["PROTEIN"], [3, 4],
                                          "PEPTIDE", 1, ["PROTEIN"], [1],
                                          170.3)

def test11():
    from pyXLMS import data
    with pytest.raises(ValueError, match = "Crosslink position has to be given for every protein! Length of proteins_b and xl_position_proteins_b has to match!"):
        crosslink = data.create_crosslink("PEPTIDE", 3, ["PROTEIN"], [3],
                                          "PEPTIDE", 1, ["PROTEIN"], [1, 2],
                                          170.3)
