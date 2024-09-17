#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

def test1():
    from pyXLMS import transform
    modifications = {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)}
    modifications_str = "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    assert transform.modifications_to_str(modifications) == modifications_str

def test2():
    from pyXLMS import transform
    modifications = {1: ("Oxidation", 15.994915)}
    modifications_str = "(1:[Oxidation|15.994915])"
    assert transform.modifications_to_str(modifications) == modifications_str

def test3():
    from pyXLMS import transform
    modifications = dict()
    modifications_str = ""
    assert transform.modifications_to_str(modifications) == modifications_str

def test4():
    from pyXLMS import transform
    assert transform.__cc([1, 2, 3]) == "1;2;3"

def test5():
    from pyXLMS import transform
    assert transform.__cc([1, 2, 3], ",") == "1,2,3"

def test6():
    from pyXLMS import transform
    assert transform.__cc([]) == ""

def test7():
    import pandas as pd
    from pyXLMS import data, transform
    c1 = data.create_crosslink("PEPTIDE", 3, ["PROTEINA"], [5], False,
                               "PEPTIDE"[::-1], 5, ["PROTEINA"], [5], True,
                               70.3)
    c2 = data.create_crosslink("PEPTIDEB", 3, ["PROTEINB"], [5], False,
                               "PEPTIDEA", 5, ["PROTEINA", "PROTEINC"], [1,2], False,
                               123.7)
    crosslinks = [c1, c2]
    df = transform.__crosslinks_to_dataframe(crosslinks)

    assert df.shape[0] == 2
    assert df.shape[1] == 11
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[0, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[0, "Alpha Proteins Crosslink Positions"] == "5"
    assert df.loc[0, "Alpha Decoy"] == True
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[0, "Beta Proteins"] == "PROTEINA"
    assert df.loc[0, "Beta Proteins Crosslink Positions"] == "5"
    assert df.loc[0, "Beta Decoy"] == False
    assert df.loc[0, "Crosslink Score"] > 70.0 and df.loc[0, "Crosslink Score"] < 71.0
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[1, "Alpha Proteins"] == "PROTEINA;PROTEINC"
    assert df.loc[1, "Alpha Proteins Crosslink Positions"] == "1;2"
    assert df.loc[1, "Alpha Decoy"] == False
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[1, "Beta Proteins"] == "PROTEINB"
    assert df.loc[1, "Beta Proteins Crosslink Positions"] == "5"
    assert df.loc[1, "Beta Decoy"] == False
    assert df.loc[1, "Crosslink Score"] > 123.0 and df.loc[0, "Crosslink Score"] < 124.0

def test8():
    import pandas as pd
    from pyXLMS import data, transform
    c1 = data.create_crosslink("PEPTIDE", 3, ["PROTEINA"], [5], False,
                               "PEPTIDE"[::-1], 5, ["PROTEINA"], [5], True,
                               70.3)
    c2 = data.create_crosslink("PEPTIDEB", 3, ["PROTEINB"], [5], False,
                               "PEPTIDEA", 5, ["PROTEINA", "PROTEINC"], [1,2], False,
                               123.7)
    crosslinks = [c1, c2]
    df = transform.to_dataframe(crosslinks)

    assert df.shape[0] == 2
    assert df.shape[1] == 11
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[0, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[0, "Alpha Proteins Crosslink Positions"] == "5"
    assert df.loc[0, "Alpha Decoy"] == True
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[0, "Beta Proteins"] == "PROTEINA"
    assert df.loc[0, "Beta Proteins Crosslink Positions"] == "5"
    assert df.loc[0, "Beta Decoy"] == False
    assert df.loc[0, "Crosslink Score"] > 70.0 and df.loc[0, "Crosslink Score"] < 71.0
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[1, "Alpha Proteins"] == "PROTEINA;PROTEINC"
    assert df.loc[1, "Alpha Proteins Crosslink Positions"] == "1;2"
    assert df.loc[1, "Alpha Decoy"] == False
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[1, "Beta Proteins"] == "PROTEINB"
    assert df.loc[1, "Beta Proteins Crosslink Positions"] == "5"
    assert df.loc[1, "Beta Decoy"] == False
    assert df.loc[1, "Crosslink Score"] > 123.0 and df.loc[0, "Crosslink Score"] < 124.0

def test9():
    # todo
    assert True

def test10():
    # todo
    assert True

def test11():
    from pyXLMS import transform
    data = [{"data_type": "peptide-spectrum-match"}]
    with pytest.raises(TypeError, match = "The given data object is not supported!"):
        df = transform.to_dataframe(data)

def test12():
    from pyXLMS import transform
    data = [{"data-type": "peptide-spectrum-match"}]
    with pytest.raises(TypeError, match = "The given data object is not supported!"):
        df = transform.to_dataframe(data)

def test12():
    from pyXLMS import transform
    data = []
    with pytest.raises(ValueError, match = "Parameter data has to be at least of length one!"):
        df = transform.to_dataframe(data)
