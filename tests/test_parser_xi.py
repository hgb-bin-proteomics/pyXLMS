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
    from pyXLMS.parser_xi import detect_xi_filetype
    import pandas as pd

    xi = pd.read_csv(XISEARCH)
    csms = pd.read_csv(XIFDR_CSMS)
    xls = pd.read_csv(XIFDR_LINKS)
    err = pd.read_excel("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")

    assert detect_xi_filetype(xi) == "xisearch"
    assert detect_xi_filetype(csms) == "xifdr_csms"
    assert detect_xi_filetype(xls) == "xifdr_crosslinks"

    with pytest.raises(
        ValueError,
        match="Could not infer data source, are you sure you read a xi result file?",
    ):
        _r = detect_xi_filetype(err)


def test2():
    from pyXLMS.parser_xi import parse_modifications_from_xi_sequence as pseq

    assert pseq("KIECcmFDSVEISGVEDR") == {4: "cm"}
    assert pseq("KIECcmFDSVEMoxISGVEDR") == {4: "cm", 10: "ox"}
    assert pseq("KIECcmFDSVEISGVEDRMox") == {4: "cm", 17: "ox"}
    assert pseq("CcmKIECcmFDSVEISGVEDRMox") == {1: "cm", 5: "cm", 18: "ox"}


def test3():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    pr = p.read_xi(XISEARCH, verbose = 0)
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
    assert mts(csm["alpha_modifications"]) == "(3:[BS3|138.06807])"
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [866]
    assert csm["alpha_proteins_peptide_positions"] == [864]
    assert csm["alpha_score"] == pytest.approx(0.596154)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "SDKNR"
    assert mts(csm["beta_modifications"]) == "(3:[BS3|138.06807])"
    assert csm["beta_peptide_crosslink_position"] == 3
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [866]
    assert csm["beta_proteins_peptide_positions"] == [864]
    assert csm["beta_score"] == pytest.approx(0.596154)
    assert not csm["beta_decoy"]
    assert csm["crosslink-type"] == "intra"
    assert csm["score"] == pytest.approx(8.758549)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 2561
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None
