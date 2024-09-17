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
    import pandas as pd
    from pyXLMS import data, transform
    c1 = data.create_csm("PEPTIDE", {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)},
                         3, ["PROTEINA"], [5], [1], 70.3, False,
                         "PEPTIDE"[::-1], {1: ("Oxidation", 15.994915)},
                         5, ["PROTEINA"], [4], [2], 20.4, True,
                         score = 70.3,
                         spectrum_file = "MS_EXP1",
                         scan_nr = 1,
                         charge = 4,
                         rt = 12.8,
                         im_cv = -50.0)
    c2 = data.create_csm("PEPTIDEB", {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)},
                         3, ["PROTEINB", "PROTEINC"], [5, 3], [1, 2], 71.3, False,
                         "PEPTIDEA", {},
                         5, ["PROTEINA"], [4], [2], 21.4, False,
                         score = 71.3,
                         spectrum_file = "MS_EXP1",
                         scan_nr = 2,
                         charge = 3,
                         rt = 12.9,
                         im_cv = -70.0)
    csms = [c1, c2]
    df = transform.__csms_to_dataframe(csms)

    assert df.shape[0] == 2
    assert df.shape[1] == 22
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Modifications"] == "(1:[Oxidation|15.994915])"
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[0, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[0, "Alpha Proteins Crosslink Positions"] == "4"
    assert df.loc[0, "Alpha Proteins Peptide Positions"] == "2"
    assert df.loc[0, "Alpha Score"] > 20.0 and df.loc[0, "Alpha Score"] < 21.0
    assert df.loc[0, "Alpha Decoy"] == True
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert df.loc[0, "Beta Peptide Modifications"] == "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[0, "Beta Proteins"] == "PROTEINA"
    assert df.loc[0, "Beta Proteins Crosslink Positions"] == "5"
    assert df.loc[0, "Beta Proteins Peptide Positions"] == "1"
    assert df.loc[0, "Beta Score"] > 70.0 and df.loc[0, "Alpha Score"] < 71.0
    assert df.loc[0, "Beta Decoy"] == False
    assert df.loc[0, "CSM Score"] > 70.0 and df.loc[0, "CSM Score"] < 71.0
    assert df.loc[0, "Spectrum File"] == "MS_EXP1"
    assert df.loc[0, "Scan Nr"] == 1
    assert df.loc[0, "Precursor Charge"] == 4
    assert df.loc[0, "Retention Time"] > 12.0 and df.loc[0, "Retention Time"] < 13.0
    assert df.loc[0, "Ion Mobility"] > -51.0 and df.loc[0, "Ion Mobility"] < -49.0
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Modifications"] == ""
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[1, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[1, "Alpha Proteins Crosslink Positions"] == "4"
    assert df.loc[1, "Alpha Proteins Peptide Positions"] == "2"
    assert df.loc[1, "Alpha Score"] > 21.0 and df.loc[1, "Alpha Score"] < 22.0
    assert df.loc[1, "Alpha Decoy"] == False
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert df.loc[1, "Beta Peptide Modifications"] == "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[1, "Beta Proteins"] == "PROTEINB;PROTEINC"
    assert df.loc[1, "Beta Proteins Crosslink Positions"] == "5;3"
    assert df.loc[1, "Beta Proteins Peptide Positions"] == "1;2"
    assert df.loc[1, "Beta Score"] > 71.0 and df.loc[1, "Alpha Score"] < 72.0
    assert df.loc[1, "Beta Decoy"] == False
    assert df.loc[1, "CSM Score"] > 71.0 and df.loc[1, "CSM Score"] < 72.0
    assert df.loc[1, "Spectrum File"] == "MS_EXP1"
    assert df.loc[1, "Scan Nr"] == 2
    assert df.loc[1, "Precursor Charge"] == 3
    assert df.loc[1, "Retention Time"] > 12.0 and df.loc[1, "Retention Time"] < 13.0
    assert df.loc[1, "Ion Mobility"] > -71.0 and df.loc[1, "Ion Mobility"] < -69.0

def test10():
    import pandas as pd
    from pyXLMS import data, transform
    c1 = data.create_csm("PEPTIDE", {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)},
                         3, ["PROTEINA"], [5], [1], 70.3, False,
                         "PEPTIDE"[::-1], {1: ("Oxidation", 15.994915)},
                         5, ["PROTEINA"], [4], [2], 20.4, True,
                         score = 70.3,
                         spectrum_file = "MS_EXP1",
                         scan_nr = 1,
                         charge = 4,
                         rt = 12.8,
                         im_cv = -50.0)
    c2 = data.create_csm("PEPTIDEB", {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)},
                         3, ["PROTEINB", "PROTEINC"], [5, 3], [1, 2], 71.3, False,
                         "PEPTIDEA", {},
                         5, ["PROTEINA"], [4], [2], 21.4, False,
                         score = 71.3,
                         spectrum_file = "MS_EXP1",
                         scan_nr = 2,
                         charge = 3,
                         rt = 12.9,
                         im_cv = -70.0)
    csms = [c1, c2]
    df = transform.to_dataframe(csms)

    assert df.shape[0] == 2
    assert df.shape[1] == 22
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Modifications"] == "(1:[Oxidation|15.994915])"
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[0, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[0, "Alpha Proteins Crosslink Positions"] == "4"
    assert df.loc[0, "Alpha Proteins Peptide Positions"] == "2"
    assert df.loc[0, "Alpha Score"] > 20.0 and df.loc[0, "Alpha Score"] < 21.0
    assert df.loc[0, "Alpha Decoy"] == True
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert df.loc[0, "Beta Peptide Modifications"] == "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[0, "Beta Proteins"] == "PROTEINA"
    assert df.loc[0, "Beta Proteins Crosslink Positions"] == "5"
    assert df.loc[0, "Beta Proteins Peptide Positions"] == "1"
    assert df.loc[0, "Beta Score"] > 70.0 and df.loc[0, "Alpha Score"] < 71.0
    assert df.loc[0, "Beta Decoy"] == False
    assert df.loc[0, "CSM Score"] > 70.0 and df.loc[0, "CSM Score"] < 71.0
    assert df.loc[0, "Spectrum File"] == "MS_EXP1"
    assert df.loc[0, "Scan Nr"] == 1
    assert df.loc[0, "Precursor Charge"] == 4
    assert df.loc[0, "Retention Time"] > 12.0 and df.loc[0, "Retention Time"] < 13.0
    assert df.loc[0, "Ion Mobility"] > -51.0 and df.loc[0, "Ion Mobility"] < -49.0
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Modifications"] == ""
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[1, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[1, "Alpha Proteins Crosslink Positions"] == "4"
    assert df.loc[1, "Alpha Proteins Peptide Positions"] == "2"
    assert df.loc[1, "Alpha Score"] > 21.0 and df.loc[1, "Alpha Score"] < 22.0
    assert df.loc[1, "Alpha Decoy"] == False
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert df.loc[1, "Beta Peptide Modifications"] == "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[1, "Beta Proteins"] == "PROTEINB;PROTEINC"
    assert df.loc[1, "Beta Proteins Crosslink Positions"] == "5;3"
    assert df.loc[1, "Beta Proteins Peptide Positions"] == "1;2"
    assert df.loc[1, "Beta Score"] > 71.0 and df.loc[1, "Alpha Score"] < 72.0
    assert df.loc[1, "Beta Decoy"] == False
    assert df.loc[1, "CSM Score"] > 71.0 and df.loc[1, "CSM Score"] < 72.0
    assert df.loc[1, "Spectrum File"] == "MS_EXP1"
    assert df.loc[1, "Scan Nr"] == 2
    assert df.loc[1, "Precursor Charge"] == 3
    assert df.loc[1, "Retention Time"] > 12.0 and df.loc[1, "Retention Time"] < 13.0
    assert df.loc[1, "Ion Mobility"] > -71.0 and df.loc[1, "Ion Mobility"] < -69.0

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
