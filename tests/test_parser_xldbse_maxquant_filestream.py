#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

MAXQUANT1 = "data/maxquant/run1/crosslinkMsms.txt"


def test1():
    from pyXLMS.parser import read_maxquant
    from pyXLMS.transform import modifications_to_str as mts

    with open(MAXQUANT1, "r", encoding="utf-8") as f:
        pr = read_maxquant(f, crosslinker="DSS")
        f.close()

    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MaxQuant"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 730

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GQKNSR"
    assert mts(csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [779]
    assert csm["alpha_proteins_peptide_positions"] == [777]
    assert csm["alpha_score"] == pytest.approx(46.617672)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "GQKNSR"
    assert mts(csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 3
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [779]
    assert csm["beta_proteins_peptide_positions"] == [777]
    assert csm["beta_score"] == pytest.approx(46.617672)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(46.618)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 2257
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[593]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "DFLELKANMAGQADAGFDGPK"
    assert mts(csm["alpha_modifications"]) == "(6:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 6
    assert csm["alpha_proteins"] == ["sp|MYG_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [20]
    assert csm["alpha_proteins_peptide_positions"] == [15]
    assert csm["alpha_score"] == pytest.approx(16.773198)
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "YVPKIK"
    assert mts(csm["beta_modifications"]) == "(4:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 4
    assert csm["beta_proteins"] == ["sp|MYG_HUMAN|"]
    assert csm["beta_proteins_crosslink_positions"] == [55]
    assert csm["beta_proteins_peptide_positions"] == [52]
    assert csm["beta_score"] == pytest.approx(22.066486)
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(13.746)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 17682
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[723]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "MDGTEELLVKLNR"
    assert (
        mts(csm["alpha_modifications"])
        == "(1:[Oxidation|15.994915]);(10:[DSS|138.06808])"
    )
    assert csm["alpha_peptide_crosslink_position"] == 10
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [396]
    assert csm["alpha_proteins_peptide_positions"] == [387]
    assert csm["alpha_score"] == pytest.approx(121.292258)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "MDGTEELLVKLNR"
    assert mts(csm["beta_modifications"]) == "(10:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 10
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [396]
    assert csm["beta_proteins_peptide_positions"] == [387]
    assert csm["beta_score"] == pytest.approx(140.526596)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(229.48)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 23071
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "MLLDIKTR"
    assert mts(csm["alpha_modifications"]) == "(6:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 6
    assert csm["alpha_proteins"] == ["sp|K1C15_SHEEP|"]
    assert csm["alpha_proteins_crosslink_positions"] == [394]
    assert csm["alpha_proteins_peptide_positions"] == [389]
    assert csm["alpha_score"] == pytest.approx(26.076371)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TEVQTGGFSKESILPK"
    assert mts(csm["beta_modifications"]) == "(10:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 10
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [1111]
    assert csm["beta_proteins_peptide_positions"] == [1102]
    assert csm["beta_score"] == pytest.approx(11.355476)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(23.711)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 23784
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None
