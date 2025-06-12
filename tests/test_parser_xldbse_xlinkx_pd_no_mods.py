#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


PD = "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3.pdResult"


def test1():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(PD, parse_modifications=False, modifications={})
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 841 + 9

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 262

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GQKNSR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [779]
    assert csm["alpha_proteins_peptide_positions"] == [777]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "GYKEVK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 3
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [1192]
    assert csm["beta_proteins_peptide_positions"] == [1190]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(125.007057923049)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 4313
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(16.0316734313965 * 60.0)
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "LNDWTLDSASKK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 11
    assert csm["alpha_proteins"] == ["decoy_98"]
    assert csm["alpha_proteins_crosslink_positions"] == [235]
    assert csm["alpha_proteins_peptide_positions"] == [225]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "NKDSR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["decoy_2"]
    assert csm["beta_proteins_crosslink_positions"] == [505]
    assert csm["beta_proteins_peptide_positions"] == [504]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(4.15420412699015)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 43290
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(88.1825561523437 * 60.0)
    assert csm["ion_mobility"] is None

    xl = crosslinks[1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "EKIEK"
    assert xl["alpha_peptide_crosslink_position"] == 2
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [443]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "KVTVK"
    assert xl["beta_peptide_crosslink_position"] == 1
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [562]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(119.059334966149)

    xl = crosslinks[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "LNDWTLDSASKK"
    assert xl["alpha_peptide_crosslink_position"] == 11
    assert xl["alpha_proteins"] == ["decoy_98"]
    assert xl["alpha_proteins_crosslink_positions"] == [235]
    assert xl["alpha_decoy"]
    assert xl["beta_peptide"] == "NKDSR"
    assert xl["beta_peptide_crosslink_position"] == 2
    assert xl["beta_proteins"] == ["decoy_2"]
    assert xl["beta_proteins_crosslink_positions"] == [505]
    assert xl["beta_decoy"]
    assert xl["crosslink_type"] == "inter"
    assert xl["score"] == pytest.approx(4.15420412699015)


def test2():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx([PD, PD], parse_modifications=False, modifications={})
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == (841 + 9) * 2

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 262 * 2

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GQKNSR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [779]
    assert csm["alpha_proteins_peptide_positions"] == [777]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "GYKEVK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 3
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [1192]
    assert csm["beta_proteins_peptide_positions"] == [1190]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(125.007057923049)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 4313
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(16.0316734313965 * 60.0)
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "LNDWTLDSASKK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 11
    assert csm["alpha_proteins"] == ["decoy_98"]
    assert csm["alpha_proteins_crosslink_positions"] == [235]
    assert csm["alpha_proteins_peptide_positions"] == [225]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "NKDSR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["decoy_2"]
    assert csm["beta_proteins_crosslink_positions"] == [505]
    assert csm["beta_proteins_peptide_positions"] == [504]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(4.15420412699015)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 43290
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(88.1825561523437 * 60.0)
    assert csm["ion_mobility"] is None

    xl = crosslinks[1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "EKIEK"
    assert xl["alpha_peptide_crosslink_position"] == 2
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [443]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "KVTVK"
    assert xl["beta_peptide_crosslink_position"] == 1
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [562]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(119.059334966149)

    xl = crosslinks[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "LNDWTLDSASKK"
    assert xl["alpha_peptide_crosslink_position"] == 11
    assert xl["alpha_proteins"] == ["decoy_98"]
    assert xl["alpha_proteins_crosslink_positions"] == [235]
    assert xl["alpha_decoy"]
    assert xl["beta_peptide"] == "NKDSR"
    assert xl["beta_peptide_crosslink_position"] == 2
    assert xl["beta_proteins"] == ["decoy_2"]
    assert xl["beta_proteins_crosslink_positions"] == [505]
    assert xl["beta_decoy"]
    assert xl["crosslink_type"] == "inter"
    assert xl["score"] == pytest.approx(4.15420412699015)
