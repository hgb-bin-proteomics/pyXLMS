#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

MAXQUANT1 = "data/maxquant/run1/crosslinkMsms.txt"
MAXQUANT2 = "data/maxquant/run2/crosslinkMsms.txt"


def test1():
    from pyXLMS.parser import parse_modifications_from_maxquant_sequence

    xl = "DSS"
    xl_mass = 138.06808
    seq1 = "_VVDELVKVM(Oxidation (M))GR_"
    seq2 = "_VVDELVKVM(Oxidation (M))GRM(Oxidation (M))_"
    seq3 = "_M(Oxidation (M))VVDELVKVM(Oxidation (M))GRM(Oxidation (M))_"

    assert parse_modifications_from_maxquant_sequence(seq1, 2, xl, xl_mass) == {2: ('DSS', 138.06808), 9: ('Oxidation', 15.994915)}
    assert parse_modifications_from_maxquant_sequence(seq2, 2, xl, xl_mass) == {2: ('DSS', 138.06808), 9: ('Oxidation', 15.994915), 12: ('Oxidation', 15.994915)}
    assert parse_modifications_from_maxquant_sequence(seq3, 2, xl, xl_mass) == {2: ('DSS', 138.06808), 1: ('Oxidation', 15.994915), 10: ('Oxidation', 15.994915), 13: ('Oxidation', 15.994915)}

    with pytest.raises(
        RuntimeError,
        match="Could not parse sequence VVDEL. Is the sequence correctly formatted?"
    ):
        _r = parse_modifications_from_maxquant_sequence("VVDEL", 2, xl, xl_mass)

    with pytest.raises(
        RuntimeError,
        match="Modification at position 9 already exists!"
    ):
        _r = parse_modifications_from_maxquant_sequence(seq1, 9, xl, xl_mass)

    with pytest.raises(
        KeyError,
        match="Key Oxi not found in parameter 'modifications'. Are you missing a modification?"
    ):
        _r = parse_modifications_from_maxquant_sequence("_VVDELVKVM(Oxi (M))GR_", 2, xl, xl_mass)


def test2():
    from pyXLMS.parser import read_maxquant

    pr = read_maxquant(MAXQUANT1, crosslinker="DSS")
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
