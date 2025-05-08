#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

SCOUT_CSMS_9 = "data/scout/Cas9_Unfiltered_CSMs.csv"
SCOUT_CSMS_F_9 = "data/scout/Cas9_Filtered_CSMs.csv"
SCOUT_XL_9 = "data/scout/Cas9_Residue_Pairs.csv"
SCOUT_CSMS_10 = "data/scout/Cas10_Unfiltered_CSMs.csv"
SCOUT_CSMS_F_10 = "data/scout/Cas10_Filtered_CSMs.csv"
SCOUT_XL_10 = "data/scout/Cas10_Residue_Pairs.csv"
SCOUT_MM = "data/scout/Multi_Mod.csv"


def test1():
    from pyXLMS.parser import detect_scout_filetype
    import pandas as pd

    csms = pd.read_csv(SCOUT_CSMS_9)
    csms_f = pd.read_csv(SCOUT_CSMS_F_9)
    xls = pd.read_csv(SCOUT_XL_9)
    err = pd.read_excel("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")

    assert detect_scout_filetype(csms) == "scout_csms_unfiltered"
    assert detect_scout_filetype(csms_f) == "scout_csms_filtered"
    assert detect_scout_filetype(xls) == "scout_xl"

    with pytest.raises(
        ValueError,
        match="Could not infer data source, are you sure you read a Scout result file?",
    ):
        _r = detect_scout_filetype(err)
