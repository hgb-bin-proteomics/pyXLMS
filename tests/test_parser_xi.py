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
    from pyXLMS.parser_xi import detect_xi_file
    import pandas as pd

    xi = pd.read_csv(XISEARCH)
    csms = pd.read_csv(XIFDR_CSMS)
    xls = pd.read_csv(XIFDR_LINKS)
    err = pd.read_excel("data/msannika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")

    assert detect_xi_file(xi) == "xisearch"
    assert detect_xi_file(csms) == "xifdr_csms"
    assert detect_xi_file(xls) == "xifdr_crosslinks"

    with pytest.raises(
        ValueError,
        match="Could not infer data source, are you sure you read a xi result file?",
    ):
        _r = detect_xi_file(err)
