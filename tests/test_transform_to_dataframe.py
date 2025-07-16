#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.transform.to_dataframe import __cc

    assert __cc([1, 2, 3]) == "1;2;3"


def test6():
    from pyXLMS.transform.to_dataframe import __cc

    assert __cc([1, 2, 3], ",") == "1,2,3"


def test7():
    from pyXLMS.transform.to_dataframe import __cc

    assert __cc([]) == ""


def test8():
    from pyXLMS.transform.to_dataframe import __cc

    assert __cc(None) is None


def test9():
    from pyXLMS import data
    from pyXLMS.transform.to_dataframe import __crosslinks_to_dataframe

    c1 = data.create_crosslink(
        "PEPTIDE",
        3,
        ["PROTEINA"],
        [5],
        False,
        "PEPTIDE"[::-1],
        5,
        ["PROTEINA"],
        [5],
        True,
        70.3,
    )
    c2 = data.create_crosslink(
        "PEPTIDEB",
        3,
        ["PROTEINB"],
        [5],
        False,
        "PEPTIDEA",
        5,
        ["PROTEINA", "PROTEINC"],
        [1, 2],
        False,
        123.7,
    )
    crosslinks = [c1, c2]
    df = __crosslinks_to_dataframe(crosslinks)

    assert df.shape[0] == 2
    assert df.shape[1] == 13
    assert df.loc[0, "Completeness"] == "full"
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[0, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[0, "Alpha Proteins Crosslink Positions"] == "5"
    # assert df.loc[0, "Alpha Decoy"] == True
    assert df.loc[0, "Alpha Decoy"]
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[0, "Beta Proteins"] == "PROTEINA"
    assert df.loc[0, "Beta Proteins Crosslink Positions"] == "5"
    # assert df.loc[0, "Beta Decoy"] == False
    assert not df.loc[0, "Beta Decoy"]
    assert df.loc[0, "Crosslink Type"] == "intra"
    assert df.loc[0, "Crosslink Score"] == pytest.approx(70.3)
    assert df.loc[1, "Completeness"] == "full"
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[1, "Alpha Proteins"] == "PROTEINA;PROTEINC"
    assert df.loc[1, "Alpha Proteins Crosslink Positions"] == "1;2"
    # assert df.loc[1, "Alpha Decoy"] == False
    assert not df.loc[1, "Alpha Decoy"]
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[1, "Beta Proteins"] == "PROTEINB"
    assert df.loc[1, "Beta Proteins Crosslink Positions"] == "5"
    # assert df.loc[1, "Beta Decoy"] == False
    assert not df.loc[1, "Beta Decoy"]
    assert df.loc[1, "Crosslink Type"] == "inter"
    assert df.loc[1, "Crosslink Score"] == pytest.approx(123.7)


def test10():
    from pyXLMS import data, transform

    c1 = data.create_crosslink(
        "PEPTIDE",
        3,
        ["PROTEINA"],
        [5],
        False,
        "PEPTIDE"[::-1],
        5,
        ["PROTEINA"],
        [5],
        True,
        70.3,
    )
    c2 = data.create_crosslink(
        "PEPTIDEB",
        3,
        ["PROTEINB"],
        [5],
        False,
        "PEPTIDEA",
        5,
        ["PROTEINA", "PROTEINC"],
        [1, 2],
        False,
        123.7,
    )
    crosslinks = [c1, c2]
    df = transform.to_dataframe(crosslinks)

    assert df.shape[0] == 2
    assert df.shape[1] == 13
    assert df.loc[0, "Completeness"] == "full"
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[0, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[0, "Alpha Proteins Crosslink Positions"] == "5"
    # assert df.loc[0, "Alpha Decoy"] == True
    assert df.loc[0, "Alpha Decoy"]
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[0, "Beta Proteins"] == "PROTEINA"
    assert df.loc[0, "Beta Proteins Crosslink Positions"] == "5"
    # assert df.loc[0, "Beta Decoy"] == False
    assert not df.loc[0, "Beta Decoy"]
    assert df.loc[0, "Crosslink Type"] == "intra"
    assert df.loc[0, "Crosslink Score"] == pytest.approx(70.3)
    assert df.loc[1, "Completeness"] == "full"
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[1, "Alpha Proteins"] == "PROTEINA;PROTEINC"
    assert df.loc[1, "Alpha Proteins Crosslink Positions"] == "1;2"
    # assert df.loc[1, "Alpha Decoy"] == False
    assert not df.loc[1, "Alpha Decoy"]
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[1, "Beta Proteins"] == "PROTEINB"
    assert df.loc[1, "Beta Proteins Crosslink Positions"] == "5"
    # assert df.loc[1, "Beta Decoy"] == False
    assert not df.loc[1, "Beta Decoy"]
    assert df.loc[1, "Crosslink Type"] == "inter"
    assert df.loc[1, "Crosslink Score"] == pytest.approx(123.7)


def test11():
    from pyXLMS import data, transform
    import pandas as pd

    c1 = data.create_crosslink(
        "PEPTIDE",
        3,
        None,
        None,
        None,
        "PEPTIDE"[::-1],
        5,
        None,
        None,
        None,
        None,
    )
    c2 = data.create_crosslink(
        "PEPTIDEB",
        3,
        None,
        None,
        None,
        "PEPTIDEA",
        5,
        None,
        None,
        None,
        None,
    )
    crosslinks = [c1, c2]
    df = transform.to_dataframe(crosslinks)

    assert df.shape[0] == 2
    assert df.shape[1] == 13
    assert df.loc[0, "Completeness"] == "partial"
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert pd.isna(df.loc[0, "Alpha Proteins"])
    assert pd.isna(df.loc[0, "Alpha Proteins Crosslink Positions"])
    assert pd.isna(df.loc[0, "Alpha Decoy"])
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert pd.isna(df.loc[0, "Beta Proteins"])
    assert pd.isna(df.loc[0, "Beta Proteins Crosslink Positions"])
    assert pd.isna(df.loc[0, "Beta Decoy"])
    assert df.loc[0, "Crosslink Type"] == "inter"
    assert pd.isna(df.loc[0, "Crosslink Score"])
    assert df.loc[1, "Completeness"] == "partial"
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert pd.isna(df.loc[1, "Alpha Proteins"])
    assert pd.isna(df.loc[1, "Alpha Proteins Crosslink Positions"])
    assert pd.isna(df.loc[1, "Alpha Decoy"])
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert pd.isna(df.loc[1, "Beta Proteins"])
    assert pd.isna(df.loc[1, "Beta Proteins Crosslink Positions"])
    assert pd.isna(df.loc[1, "Beta Decoy"])
    assert df.loc[1, "Crosslink Type"] == "inter"
    assert pd.isna(df.loc[1, "Crosslink Score"])


def test12():
    from pyXLMS import data
    from pyXLMS.transform.to_dataframe import __csms_to_dataframe

    c1 = data.create_csm(
        "PEPTIDE",
        {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)},
        3,
        ["PROTEINA"],
        [3],
        [1],
        70.3,
        False,
        "PEPTIDE"[::-1],
        {1: ("Oxidation", 15.994915)},
        5,
        ["PROTEINA"],
        [6],
        [2],
        20.4,
        True,
        score=70.3,
        spectrum_file="MS_EXP1",
        scan_nr=1,
        charge=4,
        rt=12.8,
        im_cv=-50.0,
    )
    c2 = data.create_csm(
        "PEPTIDEB",
        {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)},
        3,
        ["PROTEINB", "PROTEINC"],
        [3, 4],
        [1, 2],
        71.3,
        False,
        "PEPTIDEA",
        {},
        5,
        ["PROTEINA"],
        [6],
        [2],
        21.4,
        False,
        score=71.3,
        spectrum_file="MS_EXP1",
        scan_nr=2,
        charge=3,
        rt=12.9,
        im_cv=-70.0,
    )
    csms = [c1, c2]
    df = __csms_to_dataframe(csms)

    assert df.shape[0] == 2
    assert df.shape[1] == 24
    assert df.loc[0, "Completeness"] == "full"
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Modifications"] == "(1:[Oxidation|15.994915])"
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[0, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[0, "Alpha Proteins Crosslink Positions"] == "6"
    assert df.loc[0, "Alpha Proteins Peptide Positions"] == "2"
    assert df.loc[0, "Alpha Score"] == pytest.approx(20.4)
    # assert df.loc[0, "Alpha Decoy"] == True
    assert df.loc[0, "Alpha Decoy"]
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert (
        df.loc[0, "Beta Peptide Modifications"]
        == "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    )
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[0, "Beta Proteins"] == "PROTEINA"
    assert df.loc[0, "Beta Proteins Crosslink Positions"] == "3"
    assert df.loc[0, "Beta Proteins Peptide Positions"] == "1"
    assert df.loc[0, "Beta Score"] == pytest.approx(70.3)
    # assert df.loc[0, "Beta Decoy"] == False
    assert not df.loc[0, "Beta Decoy"]
    assert df.loc[0, "Crosslink Type"] == "intra"
    assert df.loc[0, "CSM Score"] == pytest.approx(70.3)
    assert df.loc[0, "Spectrum File"] == "MS_EXP1"
    assert df.loc[0, "Scan Nr"] == 1
    assert df.loc[0, "Precursor Charge"] == 4
    assert df.loc[0, "Retention Time"] == pytest.approx(12.8)
    assert df.loc[0, "Ion Mobility"] == pytest.approx(-50)
    assert df.loc[1, "Completeness"] == "full"
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Modifications"] == ""
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[1, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[1, "Alpha Proteins Crosslink Positions"] == "6"
    assert df.loc[1, "Alpha Proteins Peptide Positions"] == "2"
    assert df.loc[1, "Alpha Score"] == pytest.approx(21.4)
    # assert df.loc[1, "Alpha Decoy"] == False
    assert not df.loc[1, "Alpha Decoy"]
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert (
        df.loc[1, "Beta Peptide Modifications"]
        == "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    )
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[1, "Beta Proteins"] == "PROTEINB;PROTEINC"
    assert df.loc[1, "Beta Proteins Crosslink Positions"] == "3;4"
    assert df.loc[1, "Beta Proteins Peptide Positions"] == "1;2"
    assert df.loc[1, "Beta Score"] == pytest.approx(71.3)
    # assert df.loc[1, "Beta Decoy"] == False
    assert not df.loc[1, "Beta Decoy"]
    assert df.loc[1, "Crosslink Type"] == "inter"
    assert df.loc[1, "CSM Score"] == pytest.approx(71.3)
    assert df.loc[1, "Spectrum File"] == "MS_EXP1"
    assert df.loc[1, "Scan Nr"] == 2
    assert df.loc[1, "Precursor Charge"] == 3
    assert df.loc[1, "Retention Time"] == pytest.approx(12.9)
    assert df.loc[1, "Ion Mobility"] == pytest.approx(-70)


def test13():
    from pyXLMS import data, transform

    c1 = data.create_csm(
        "PEPTIDE",
        {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)},
        3,
        ["PROTEINA"],
        [3],
        [1],
        70.3,
        False,
        "PEPTIDE"[::-1],
        {1: ("Oxidation", 15.994915)},
        5,
        ["PROTEINA"],
        [6],
        [2],
        20.4,
        True,
        score=70.3,
        spectrum_file="MS_EXP1",
        scan_nr=1,
        charge=4,
        rt=12.8,
        im_cv=-50.0,
    )
    c2 = data.create_csm(
        "PEPTIDEB",
        {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)},
        3,
        ["PROTEINB", "PROTEINC"],
        [3, 4],
        [1, 2],
        71.3,
        False,
        "PEPTIDEA",
        {},
        5,
        ["PROTEINA"],
        [6],
        [2],
        21.4,
        False,
        score=71.3,
        spectrum_file="MS_EXP1",
        scan_nr=2,
        charge=3,
        rt=12.9,
        im_cv=-70.0,
    )
    csms = [c1, c2]
    df = transform.to_dataframe(csms)

    assert df.shape[0] == 2
    assert df.shape[1] == 24
    assert df.loc[0, "Completeness"] == "full"
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert df.loc[0, "Alpha Peptide Modifications"] == "(1:[Oxidation|15.994915])"
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[0, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[0, "Alpha Proteins Crosslink Positions"] == "6"
    assert df.loc[0, "Alpha Proteins Peptide Positions"] == "2"
    assert df.loc[0, "Alpha Score"] == pytest.approx(20.4)
    # assert df.loc[0, "Alpha Decoy"] == True
    assert df.loc[0, "Alpha Decoy"]
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert (
        df.loc[0, "Beta Peptide Modifications"]
        == "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    )
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[0, "Beta Proteins"] == "PROTEINA"
    assert df.loc[0, "Beta Proteins Crosslink Positions"] == "3"
    assert df.loc[0, "Beta Proteins Peptide Positions"] == "1"
    assert df.loc[0, "Beta Score"] == pytest.approx(70.3)
    # assert df.loc[0, "Beta Decoy"] == False
    assert not df.loc[0, "Beta Decoy"]
    assert df.loc[0, "Crosslink Type"] == "intra"
    assert df.loc[0, "CSM Score"] == pytest.approx(70.3)
    assert df.loc[0, "Spectrum File"] == "MS_EXP1"
    assert df.loc[0, "Scan Nr"] == 1
    assert df.loc[0, "Precursor Charge"] == 4
    assert df.loc[0, "Retention Time"] == pytest.approx(12.8)
    assert df.loc[0, "Ion Mobility"] == pytest.approx(-50)
    assert df.loc[1, "Completeness"] == "full"
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert df.loc[1, "Alpha Peptide Modifications"] == ""
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert df.loc[1, "Alpha Proteins"] == "PROTEINA"
    assert df.loc[1, "Alpha Proteins Crosslink Positions"] == "6"
    assert df.loc[1, "Alpha Proteins Peptide Positions"] == "2"
    assert df.loc[1, "Alpha Score"] == pytest.approx(21.4)
    # assert df.loc[1, "Alpha Decoy"] == False
    assert not df.loc[1, "Alpha Decoy"]
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert (
        df.loc[1, "Beta Peptide Modifications"]
        == "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    )
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert df.loc[1, "Beta Proteins"] == "PROTEINB;PROTEINC"
    assert df.loc[1, "Beta Proteins Crosslink Positions"] == "3;4"
    assert df.loc[1, "Beta Proteins Peptide Positions"] == "1;2"
    assert df.loc[1, "Beta Score"] == pytest.approx(71.3)
    # assert df.loc[1, "Beta Decoy"] == False
    assert not df.loc[1, "Beta Decoy"]
    assert df.loc[1, "Crosslink Type"] == "inter"
    assert df.loc[1, "CSM Score"] == pytest.approx(71.3)
    assert df.loc[1, "Spectrum File"] == "MS_EXP1"
    assert df.loc[1, "Scan Nr"] == 2
    assert df.loc[1, "Precursor Charge"] == 3
    assert df.loc[1, "Retention Time"] == pytest.approx(12.9)
    assert df.loc[1, "Ion Mobility"] == pytest.approx(-70.0)


def test14():
    from pyXLMS import data, transform
    import pandas as pd

    c1 = data.create_csm(
        "PEPTIDE",
        None,
        3,
        None,
        None,
        None,
        None,
        None,
        "PEPTIDE"[::-1],
        None,
        5,
        None,
        None,
        None,
        None,
        None,
        score=None,
        spectrum_file="MS_EXP1",
        scan_nr=1,
        charge=None,
        rt=None,
        im_cv=None,
    )
    c2 = data.create_csm(
        "PEPTIDEB",
        None,
        3,
        None,
        None,
        None,
        None,
        None,
        "PEPTIDEA",
        None,
        5,
        None,
        None,
        None,
        None,
        None,
        score=None,
        spectrum_file="MS_EXP1",
        scan_nr=2,
        charge=None,
        rt=None,
        im_cv=None,
    )
    csms = [c1, c2]
    df = transform.to_dataframe(csms)

    assert df.shape[0] == 2
    assert df.shape[1] == 24
    assert df.loc[0, "Completeness"] == "partial"
    assert df.loc[0, "Alpha Peptide"] == "PEPTIDE"[::-1]
    assert pd.isna(df.loc[0, "Alpha Peptide Modifications"])
    assert df.loc[0, "Alpha Peptide Crosslink Position"] == 5
    assert pd.isna(df.loc[0, "Alpha Proteins"])
    assert pd.isna(df.loc[0, "Alpha Proteins Crosslink Positions"])
    assert pd.isna(df.loc[0, "Alpha Proteins Peptide Positions"])
    assert pd.isna(df.loc[0, "Alpha Score"])
    assert pd.isna(df.loc[0, "Alpha Decoy"])
    assert df.loc[0, "Beta Peptide"] == "PEPTIDE"
    assert pd.isna(df.loc[0, "Beta Peptide Modifications"])
    assert df.loc[0, "Beta Peptide Crosslink Position"] == 3
    assert pd.isna(df.loc[0, "Beta Proteins"])
    assert pd.isna(df.loc[0, "Beta Proteins Crosslink Positions"])
    assert pd.isna(df.loc[0, "Beta Proteins Peptide Positions"])
    assert pd.isna(df.loc[0, "Beta Score"])
    assert pd.isna(df.loc[0, "Beta Decoy"])
    assert df.loc[0, "Crosslink Type"] == "inter"
    assert pd.isna(df.loc[0, "CSM Score"])
    assert df.loc[0, "Spectrum File"] == "MS_EXP1"
    assert df.loc[0, "Scan Nr"] == 1
    assert pd.isna(df.loc[0, "Precursor Charge"])
    assert pd.isna(df.loc[0, "Retention Time"])
    assert pd.isna(df.loc[0, "Ion Mobility"])
    assert df.loc[1, "Completeness"] == "partial"
    assert df.loc[1, "Alpha Peptide"] == "PEPTIDEA"
    assert pd.isna(df.loc[1, "Alpha Peptide Modifications"])
    assert df.loc[1, "Alpha Peptide Crosslink Position"] == 5
    assert pd.isna(df.loc[1, "Alpha Proteins"])
    assert pd.isna(df.loc[1, "Alpha Proteins Crosslink Positions"])
    assert pd.isna(df.loc[1, "Alpha Proteins Peptide Positions"])
    assert pd.isna(df.loc[1, "Alpha Score"])
    assert pd.isna(df.loc[1, "Alpha Decoy"])
    assert df.loc[1, "Beta Peptide"] == "PEPTIDEB"
    assert pd.isna(df.loc[1, "Beta Peptide Modifications"])
    assert df.loc[1, "Beta Peptide Crosslink Position"] == 3
    assert pd.isna(df.loc[1, "Beta Proteins"])
    assert pd.isna(df.loc[1, "Beta Proteins Crosslink Positions"])
    assert pd.isna(df.loc[1, "Beta Proteins Peptide Positions"])
    assert pd.isna(df.loc[1, "Beta Score"])
    assert pd.isna(df.loc[1, "Beta Decoy"])
    assert df.loc[1, "Crosslink Type"] == "inter"
    assert pd.isna(df.loc[1, "CSM Score"])
    assert df.loc[1, "Spectrum File"] == "MS_EXP1"
    assert df.loc[1, "Scan Nr"] == 2
    assert pd.isna(df.loc[1, "Precursor Charge"])
    assert pd.isna(df.loc[1, "Retention Time"])
    assert pd.isna(df.loc[1, "Ion Mobility"])


def test15():
    from pyXLMS import transform

    data = [{"data_type": "peptide-spectrum-match"}]
    with pytest.raises(TypeError, match="The given data object is not supported!"):
        _df = transform.to_dataframe(data)


def test16():
    from pyXLMS import transform

    data = [{"data-type": "peptide-spectrum-match"}]
    with pytest.raises(TypeError, match="The given data object is not supported!"):
        _df = transform.to_dataframe(data)


def test17():
    from pyXLMS import transform

    data = []
    with pytest.raises(
        ValueError, match="Parameter data has to be at least of length one!"
    ):
        _df = transform.to_dataframe(data)
