#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


PD = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult"


def test1():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika(PD, parse_modifications=False, modifications={})
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    csm = csms[822]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KIECFDSVEISGVEDR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [575]
    assert csm["alpha_proteins_peptide_positions"] == [575]
    assert csm["alpha_score"] == pytest.approx(376.15, abs=0.01)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87, abs=0.01)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15, abs=0.01)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995, abs=0.01)
    assert csm["ion_mobility"] == pytest.approx(0.0, abs=0.01)

    xl = crosslinks[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "GQKNSR"
    assert xl["alpha_peptide_crosslink_position"] == 3
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [779]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "GQKNSR"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [779]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(119.83, abs=0.01)
