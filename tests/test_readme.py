#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    import pyXLMS

    pr = pyXLMS.parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    _ = pyXLMS.transform.summary(pr)
    _ = pyXLMS.exporter.to_xifdr(
        pr["crosslink-spectrum-matches"], filename="msannika_CSMs_for_xiFDR.csv"
    )
    assert pr is not None
