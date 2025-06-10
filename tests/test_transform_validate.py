#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate
    pr = read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 826
    validated = validate(csms)
    assert len(validated) == 705


def test2():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate
    pr = read(["data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx"], engine="MS Annika", crosslinker="DSS")
    assert len(pr["crosslink-spectrum-matches"]) == 826
    assert len(pr["crosslinks"]) == 300
    validated = validate(pr)
    assert len(validated["crosslink-spectrum-matches"]) == 705
    assert len(validated["crosslinks"]) == 226


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate
    pr = read(["data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx"], engine="MS Annika", crosslinker="DSS")
    assert len(pr["crosslink-spectrum-matches"]) == 826
    assert len(pr["crosslinks"]) == 300
    validated = validate(pr, fdr=0.05)
    assert len(validated["crosslink-spectrum-matches"]) == 825
    assert len(validated["crosslinks"]) == 260
