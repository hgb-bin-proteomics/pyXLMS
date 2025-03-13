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
    err = pd.read_excel("data/msannika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")

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

    assert pseq("KIECcmFDSVEISGVEDR") == {4: "Ccm"}
    assert pseq("KIECcmFDSVEMoxISGVEDR") == {4: "Ccm", 10: "Mox"}
    assert pseq("KIECcmFDSVEISGVEDRMox") == {4: "Ccm", 17: "Mox"}
    assert pseq("CcmKIECcmFDSVEISGVEDRMox") == {1: "Ccm", 5: "Ccm", 18: "Mox"}
