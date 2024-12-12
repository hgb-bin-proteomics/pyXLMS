#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


XL_TSV = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt"
XL_XLSX = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx"
CSMS_TSV = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt"
CSMS_XLSX = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"


def test1():
    from pyXLMS import parser as p

    parser_result = p.read_msannika(XL_TSV)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslink[299]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink-type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink-type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)
