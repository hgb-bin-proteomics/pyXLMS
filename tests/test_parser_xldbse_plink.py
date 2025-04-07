#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

PLINK2 = "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv"
PLINK3 = "data/plink3/Cas10_plus10_2025.04.07.filtered_cross-linked_spectra.csv"


def test1():
    from pyXLMS.parser import parse_spectrum_file_from_plink

    assert (
        parse_spectrum_file_from_plink(
            "XLpeplib_Beveridge_QEx-HFX_DSS_R1.20588.20588.3.0.dta"
        )
        == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    )


def test2():
    from pyXLMS.parser import parse_scan_nr_from_plink

    assert (
        parse_scan_nr_from_plink(
            "XLpeplib_Beveridge_QEx-HFX_DSS_R1.20588.20588.3.0.dta"
        )
        == 20588
    )


def test3():
    from pyXLMS.parser import read_plink
    from pyXLMS.transform import modifications_to_str as mts

    pr = read_plink(PLINK2)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "pLink"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 961

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "FDNLTKAER"
    assert mts(csm["alpha_modifications"]) == "(6:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 6
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [906]
    assert csm["alpha_proteins_peptide_positions"] == [901]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "YDENDKLIR"
    assert mts(csm["beta_modifications"]) == "(6:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [952]
    assert csm["beta_proteins_peptide_positions"] == [947]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(float("5.553153e-009"))
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 13098
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None
