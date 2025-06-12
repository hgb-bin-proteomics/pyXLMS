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
    from pyXLMS import parser as p

    pr = p.read_scout(
        SCOUT_CSMS_9,
        crosslinker="DSSO",
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 1697


def test2():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    pr = p.read_scout(
        SCOUT_CSMS_10,
        crosslinker="DSSO",
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
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
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 10
    assert csm["alpha_proteins"] == ["Cas10", "Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.405408, abs=0.00001)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "MLASAGELQKGNELALPSK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 10
    assert csm["beta_proteins"] == ["Cas10", "Cas9"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.390379, abs=0.00001)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(0.390379, abs=0.00001)
    assert (
        csm["spectrum_file"]
        == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    )
    assert csm["scan_nr"] == 21781
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[1668]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "CFQWQRNMRKVR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 10
    assert csm["alpha_proteins"] == ["spTRFL_HUMAN_"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.01438, abs=0.00001)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "IYEGEKK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["spSUMO1_HUMAN_"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.0231, abs=0.00001)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.01438, abs=0.00001)
    assert (
        csm["spectrum_file"]
        == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    )
    assert csm["scan_nr"] == 28673
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[1685]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "HTKLFDK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] == ["spPEPA_PIG_"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.031, abs=0.00001)
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "SGANGTKTSEENGGKGLDDAK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 15
    assert csm["beta_proteins"] == ["spSODC_HUMAN_"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.036907, abs=0.00001)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.031, abs=0.00001)
    assert (
        csm["spectrum_file"]
        == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    )
    assert csm["scan_nr"] == 30723
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[1689]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KLVDSTDK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas10", "Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.125987, abs=0.00001)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "SSSYHKSSSYRVSM"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["spK1C10_HUMAN_"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.01305, abs=0.00001)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.01305, abs=0.00001)
    assert (
        csm["spectrum_file"]
        == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    )
    assert csm["scan_nr"] == 31150
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "AKIQDKEGIPPDQQR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 6
    assert csm["alpha_proteins"] == ["spRS27A_HUMAN_"]
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] == pytest.approx(0.005108, abs=0.00001)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "ICSVPPGRVKRMNR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 10
    assert csm["beta_proteins"] == ["spTRFL_HUMAN_"]
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] == pytest.approx(0.024458, abs=0.00001)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.005108, abs=0.00001)
    assert (
        csm["spectrum_file"]
        == "C:\\Users\\P42587\\Downloads\\scout\\XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    )
    assert csm["scan_nr"] == 36420
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test3():
    from pyXLMS import parser as p

    pr = p.read_scout(
        SCOUT_CSMS_F_9,
        crosslinker="DSSO",
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 1306


def test4():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    pr = p.read_scout(
        SCOUT_CSMS_F_10,
        crosslinker="DSSO",
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 1305

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "MLASAGELQKGNELALPSK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 10
    assert csm["alpha_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [1237, 1226]
    assert csm["alpha_proteins_peptide_positions"] == [1228, 1217]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "MLASAGELQKGNELALPSK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 10
    assert csm["beta_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [1237, 1226]
    assert csm["beta_proteins_peptide_positions"] == [1228, 1217]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(0.390379, abs=0.00001)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    assert csm["scan_nr"] == 21781
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[91]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KIECFDSVEISGVEDR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [586, 575]
    assert csm["alpha_proteins_peptide_positions"] == [586, 575]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [699, 688]
    assert csm["beta_proteins_peptide_positions"] == [693, 682]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(0.278628, abs=0.00001)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    assert csm["scan_nr"] == 29092
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "DKQSGK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [688, 677]
    assert csm["alpha_proteins_peptide_positions"] == [687, 676]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "SVKELLGITIMER"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 3
    assert csm["beta_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [1176, 1165]
    assert csm["beta_proteins_peptide_positions"] == [1174, 1163]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(0.097123, abs=0.00001)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    assert csm["scan_nr"] == 20401
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test5():
    from pyXLMS import parser as p

    pr = p.read_scout(
        SCOUT_XL_9,
        crosslinker="DSSO",
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is None
    assert pr["crosslinks"] is not None

    crosslinks = pr["crosslinks"]
    assert len(crosslinks) == 200


def test6():
    from pyXLMS import parser as p

    pr = p.read_scout(
        SCOUT_XL_10,
        crosslinker="DSSO",
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is None
    assert pr["crosslinks"] is not None

    crosslinks = pr["crosslinks"]
    assert len(crosslinks) == 200

    xl = crosslinks[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "MLASAGELQKGNELALPSK"
    assert xl["alpha_peptide_crosslink_position"] == 10
    assert xl["alpha_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [1237, 1226]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "MLASAGELQKGNELALPSK"
    assert xl["beta_peptide_crosslink_position"] == 10
    assert xl["beta_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [1237, 1226]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(0.999998, abs=0.00001)

    xl = crosslinks[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "LPKYSLFELENGR"
    assert xl["alpha_peptide_crosslink_position"] == 3
    assert xl["alpha_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [1215, 1204]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "NPIDFLEAKGYK"
    assert xl["beta_peptide_crosslink_position"] == 9
    assert xl["beta_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [1200, 1189]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(0.999892, abs=0.00001)


def test7():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    pr = p.read_scout(
        SCOUT_MM,
        crosslinker="DSSO",
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 2266

    csm = csms[2241]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "ESILPKR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 6
    assert csm["alpha_proteins"] == ["sp|Cas10|Cas10", "sp|Cas9|Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [1128, 1117]
    assert csm["alpha_proteins_peptide_positions"] == [1123, 1112]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "ICSVPPGRVKRMNR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 10
    assert csm["beta_proteins"] == ["sp|spTRFL_HUMAN_|spTRFL_HUMAN_"]
    assert csm["beta_proteins_crosslink_positions"] == [1374]
    assert csm["beta_proteins_peptide_positions"] == [1365]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.037393, abs=0.00001)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw"
    assert csm["scan_nr"] == 17503
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test8():
    from pyXLMS import parser as p

    pr = p.read_scout(
        [SCOUT_MM, SCOUT_XL_9, SCOUT_XL_10, SCOUT_CSMS_10],
        crosslinker="DSSO",
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "full"
    assert pr["search_engine"] == "Scout"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is not None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 2266 + 1696

    crosslinks = pr["crosslinks"]
    assert len(crosslinks) == 200 + 200
