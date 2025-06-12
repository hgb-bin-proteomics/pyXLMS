#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

XISEARCH = "data/xi/r1_Xi1.7.6.7.csv"
XIFDR_CSMS = "data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv"
XIFDR_LINKS = "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv"


def test1():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    pr = p.read_xi(XISEARCH, parse_modifications=False, modifications={}, verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "xiSearch/xiFDR"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 4648

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "SDKNR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [866]
    assert csm["alpha_proteins_peptide_positions"] == [864]
    assert csm["alpha_score"] == pytest.approx(0.596154)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "SDKNR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 3
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [866]
    assert csm["beta_proteins_peptide_positions"] == [864]
    assert csm["beta_score"] == pytest.approx(0.596154)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(8.758549)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 2561
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[29]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "CCPSR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["KRA3_SHEEP"]
    assert csm["alpha_proteins_crosslink_positions"] == [1]
    assert csm["alpha_proteins_peptide_positions"] == [1]
    assert csm["alpha_score"] == pytest.approx(0.240385)
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "SKKLK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 3
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [213]
    assert csm["beta_proteins_peptide_positions"] == [211]
    assert csm["beta_score"] == pytest.approx(0.355769)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(1.468129)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 6482
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[184]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "EKVYFQFK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [394.0]
    assert csm["alpha_proteins_peptide_positions"] == [393]
    assert csm["alpha_score"] == pytest.approx(0.214286)
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "MLRNER"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["RETBP_HUMAN"]
    assert csm["beta_proteins_crosslink_positions"] == [1]
    assert csm["beta_proteins_peptide_positions"] == [1]
    assert csm["beta_score"] == pytest.approx(0.138462)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(-0.014849)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 8316
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "EYNSAMKR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 7
    assert csm["alpha_proteins"] == ["MYG_HUMAN"]
    assert csm["alpha_proteins_crosslink_positions"] == [13]
    assert csm["alpha_proteins_peptide_positions"] == [7]
    assert csm["alpha_score"] == pytest.approx(0.0)
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "LKGSPEDNEQKQLFVEQHK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [1250]
    assert csm["beta_proteins_peptide_positions"] == [1249]
    assert csm["beta_score"] == pytest.approx(0.126068)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(0.48559)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 27092
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test2():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    pr = p.read_xi(XIFDR_CSMS, parse_modifications=False, modifications={}, verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "xiSearch/xiFDR"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 413

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KIECFDSVEISGVEDR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [575]
    assert csm["alpha_proteins_peptide_positions"] == [575]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KIECFDSVEISGVEDR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [575]
    assert csm["beta_proteins_peptide_positions"] == [575]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(27.268)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 19140
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[332]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NPIDFLEAKGYK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 9
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [1189]
    assert csm["alpha_proteins_peptide_positions"] == [1181]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "QWYKNK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 4
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [488]
    assert csm["beta_proteins_peptide_positions"] == [485]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(10.071)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 21757
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "QITKHVAQILDSR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 4
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [933]
    assert csm["alpha_proteins_peptide_positions"] == [930]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "QQLPEKYK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [350]
    assert csm["beta_proteins_peptide_positions"] == [345]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(6.808)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 12190
    assert csm["charge"] == 5
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test3():
    from pyXLMS import parser as p

    pr = p.read_xi(XIFDR_LINKS, parse_modifications=False, modifications={}, verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "xiSearch/xiFDR"
    assert pr["crosslink-spectrum-matches"] is None
    assert pr["crosslinks"] is not None

    xls = pr["crosslinks"]
    assert len(xls) == 227

    xl = xls[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "VVDELVKVMGR"
    assert xl["alpha_peptide_crosslink_position"] == 7
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [753]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "VVDELVKVMGR"
    assert xl["beta_peptide_crosslink_position"] == 7
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [753]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(40.679)

    xl = xls[30]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "FDNLTKAER"
    assert xl["alpha_peptide_crosslink_position"] == 6
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [952]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "YDENDKLIR"
    assert xl["beta_peptide_crosslink_position"] == 6
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [906]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(25.616)

    xl = xls[219]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "NPIDFLEAKGYK"
    assert xl["alpha_peptide_crosslink_position"] == 9
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [488]
    assert xl["alpha_decoy"]
    assert xl["beta_peptide"] == "QWYKNK"
    assert xl["beta_peptide_crosslink_position"] == 4
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [1189]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(10.071)

    xl = xls[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "LSKSR"
    assert xl["alpha_peptide_crosslink_position"] == 3
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [222]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "LSKSR"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [222]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(9.619)


def test4():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    pr = p.read_xi(
        [XIFDR_LINKS, XIFDR_CSMS],
        parse_modifications=False,
        modifications={},
        verbose=0,
    )
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "full"
    assert pr["search_engine"] == "xiSearch/xiFDR"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is not None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 413

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KIECFDSVEISGVEDR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [575]
    assert csm["alpha_proteins_peptide_positions"] == [575]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KIECFDSVEISGVEDR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [575]
    assert csm["beta_proteins_peptide_positions"] == [575]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(27.268)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 19140
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[332]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NPIDFLEAKGYK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 9
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [1189]
    assert csm["alpha_proteins_peptide_positions"] == [1181]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "QWYKNK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 4
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [488]
    assert csm["beta_proteins_peptide_positions"] == [485]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(10.071)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 21757
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "QITKHVAQILDSR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 4
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [933]
    assert csm["alpha_proteins_peptide_positions"] == [930]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "QQLPEKYK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [350]
    assert csm["beta_proteins_peptide_positions"] == [345]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(6.808)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 12190
    assert csm["charge"] == 5
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    xls = pr["crosslinks"]
    assert len(xls) == 227

    xl = xls[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "VVDELVKVMGR"
    assert xl["alpha_peptide_crosslink_position"] == 7
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [753]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "VVDELVKVMGR"
    assert xl["beta_peptide_crosslink_position"] == 7
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [753]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(40.679)

    xl = xls[30]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "FDNLTKAER"
    assert xl["alpha_peptide_crosslink_position"] == 6
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [952]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "YDENDKLIR"
    assert xl["beta_peptide_crosslink_position"] == 6
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [906]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(25.616)

    xl = xls[219]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "NPIDFLEAKGYK"
    assert xl["alpha_peptide_crosslink_position"] == 9
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [488]
    assert xl["alpha_decoy"]
    assert xl["beta_peptide"] == "QWYKNK"
    assert xl["beta_peptide_crosslink_position"] == 4
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [1189]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(10.071)

    xl = xls[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "LSKSR"
    assert xl["alpha_peptide_crosslink_position"] == 3
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [222]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "LSKSR"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [222]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(9.619)
