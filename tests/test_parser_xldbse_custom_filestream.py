#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    with open("data/pyxlms/csm.txt", "r", encoding="utf-8") as f:
        parser_result = p.read_custom(f, format="csv")
        f.close()

    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 2

    first_csm = csms[0]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "KPEPTIDE"
    assert mts(first_csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [13]
    assert first_csm["alpha_proteins_peptide_positions"] == [13]
    assert first_csm["alpha_score"] == pytest.approx(87.53)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "PEPKTIDE"
    assert mts(first_csm["beta_modifications"]) == "(4:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 4
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [17]
    assert first_csm["beta_proteins_peptide_positions"] == [14]
    assert first_csm["beta_score"] == pytest.approx(100.3)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(87.53)
    assert first_csm["spectrum_file"] == "S1.raw"
    assert first_csm["scan_nr"] == 1
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(14.3)
    assert first_csm["ion_mobility"] == pytest.approx(50.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "EKTIDEM"
    assert (
        mts(last_csm["alpha_modifications"])
        == "(2:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    )
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [33, 21]
    assert last_csm["alpha_proteins_peptide_positions"] == [32, 20]
    assert last_csm["alpha_score"] == pytest.approx(5.3)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "PEKPIDE"
    assert mts(last_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 3
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [28]
    assert last_csm["beta_proteins_peptide_positions"] == [26]
    assert last_csm["beta_score"] == pytest.approx(34.89)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(5.4)
    assert last_csm["spectrum_file"] == "S1.raw"
    assert last_csm["scan_nr"] == 2
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(37.332)
    assert last_csm["ion_mobility"] == pytest.approx(-70.0)
