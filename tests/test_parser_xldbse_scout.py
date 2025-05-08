#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

SCOUT_CSMS_9 = "data/scout/Cas9_Unfiltered_CSMs.csv"
SCOUT_CSMS_F_9 = "data/scout/Cas9_Filtered_CSMs.csv"
SCOUT_XL_9 = "data/scout/Cas9_Residue_Pairs.csv"
SCOUT_CSMS_10 = "data/scout/Cas10_Unfiltered_CSMs.csv"
SCOUT_CSMS_F_10 = "data/scout/Cas10_Filtered_CSMs.csv"
SCOUT_XL_10 = "data/scout/Cas10_Residue_Pairs.csv"
SCOUT_MM = "data/scout/Multi_Mod.csv"


def test1():
    from pyXLMS.parser import detect_scout_filetype
    import pandas as pd

    csms = pd.read_csv(SCOUT_CSMS_9)
    csms_f = pd.read_csv(SCOUT_CSMS_F_9)
    xls = pd.read_csv(SCOUT_XL_9)
    err = pd.read_excel("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")

    assert detect_scout_filetype(csms) == "scout_csms_unfiltered"
    assert detect_scout_filetype(csms_f) == "scout_csms_filtered"
    assert detect_scout_filetype(xls) == "scout_xl"

    with pytest.raises(
        ValueError,
        match="Could not infer data source, are you sure you read a Scout result file?",
    ):
        _r = detect_scout_filetype(err)


def test2():
    from pyXLMS.parser import parse_modifications_from_scout_sequence as pms

    seq = "M(+15.994900)LASAGELQKGNELALPSK"
    assert pms(seq, 10, "DSS", 138.06808) == {10: ('DSS', 138.06808), 1: ('Oxidation', 15.994915)}

    seq = "KIEC(+57.021460)FDSVEISGVEDR"
    assert pms(seq, 1, "DSS", 138.06808) == {1: ('DSS', 138.06808), 4: ('Carbamidomethyl', 57.021464)}


def test3():
    from pyXLMS import parser as p

    with pytest.raises(TypeError, match="Verbose level has to be one of 0, 1, or 2!"):
        _r = p.read_scout(SCOUT_CSMS_9, verbose=3)


def test4():
    from pyXLMS import parser as p

    pr = p.read_scout(SCOUT_CSMS_9, crosslinker="DSSO", verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 1697


def test5():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    pr = p.read_scout(SCOUT_CSMS_10, crosslinker="DSSO", verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 1696

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "MLASAGELQKGNELALPSK"
    assert mts(csm["alpha_modifications"]) == "(1:[Oxidation|15.994915]);(10:[DSSO|158.00376])"
    assert csm["alpha_peptide_crosslink_position"] == 10
    assert csm["alpha_proteins"] == ["Cas10", "Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.405408)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "MLASAGELQKGNELALPSK"
    assert mts(csm["beta_modifications"]) == "(10:[DSSO|158.00376])"
    assert csm["beta_peptide_crosslink_position"] == 10
    assert csm["beta_proteins"] == ["Cas10", "Cas9"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.390379)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(0.390379)
    assert csm["spectrum_file"] == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    assert csm["scan_nr"] == 21781
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[1668]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "CFQWQRNMRKVR"
    assert mts(csm["alpha_modifications"]) == "(1:[Carbamidomethyl|57.021464]);(8:[Oxidation|15.994915]);(10:[DSSO|158.00376])"
    assert csm["alpha_peptide_crosslink_position"] == 10
    assert csm["alpha_proteins"] == ["spTRFL_HUMAN_"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.01438)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "IYEGEKK"
    assert mts(csm["beta_modifications"]) == "(6:[DSSO|158.00376])"
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["spSUMO1_HUMAN_"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.0231)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.01438)
    assert csm["spectrum_file"] == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    assert csm["scan_nr"] == 28673
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[1685]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "HTKLFDK"
    assert mts(csm["alpha_modifications"]) == "(3:[DSSO|158.00376])"
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] == ["spPEPA_PIG_"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.031)
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "SGANGTKTSEENGGKGLDDAK"
    assert mts(csm["beta_modifications"]) == "(15:[DSSO|158.00376])"
    assert csm["beta_peptide_crosslink_position"] == 15
    assert csm["beta_proteins"] == ["spSODC_HUMAN_"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.036907)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.031)
    assert csm["spectrum_file"] == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    assert csm["scan_nr"] == 30723
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[1689]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KLVDSTDK"
    assert mts(csm["alpha_modifications"]) == "(1:[DSSO|158.00376])"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas10", "Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.125987)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "SSSYHKSSSYRVSM"
    assert mts(csm["beta_modifications"]) == "(6:[DSSO|158.00376])"
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["spK1C10_HUMAN_"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.01305)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.01305)
    assert csm["spectrum_file"] == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    assert csm["scan_nr"] == 31150
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
