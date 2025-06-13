#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.parser import read
    from pyXLMS.transform import summary

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    stats = summary(csms)
    assert stats["Number of CSMs"] == pytest.approx(826.0)
    assert stats["Number of unique CSMs"] == pytest.approx(826.0)
    assert stats["Number of intra CSMs"] == pytest.approx(803.0)
    assert stats["Number of inter CSMs"] == pytest.approx(23.0)
    assert stats["Number of target-target CSMs"] == pytest.approx(786.0)
    assert stats["Number of target-decoy CSMs"] == pytest.approx(39.0)
    assert stats["Number of decoy-decoy CSMs"] == pytest.approx(1.0)
    assert stats["Minimum CSM score"] == pytest.approx(1.11)
    assert stats["Maximum CSM score"] == pytest.approx(452.99)


def test2():
    from pyXLMS.parser import read
    from pyXLMS.transform import summary

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    stats = summary(pr)
    assert stats["Number of crosslinks"] == pytest.approx(300.0)
    assert stats["Number of unique crosslinks by peptide"] == pytest.approx(300.0)
    assert stats["Number of unique crosslinks by protein"] == pytest.approx(298.0)
    assert stats["Number of intra crosslinks"] == pytest.approx(279.0)
    assert stats["Number of inter crosslinks"] == pytest.approx(21.0)
    assert stats["Number of target-target crosslinks"] == pytest.approx(265.0)
    assert stats["Number of target-decoy crosslinks"] == pytest.approx(0.0)
    assert stats["Number of decoy-decoy crosslinks"] == pytest.approx(35.0)
    assert stats["Minimum crosslink score"] == pytest.approx(1.11)
    assert stats["Maximum crosslink score"] == pytest.approx(452.99)


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import summary

    pr = read(
        [
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        ],
        engine="MS Annika",
        crosslinker="DSS",
    )
    stats = summary(pr)
    assert stats["Number of CSMs"] == pytest.approx(826.0)
    assert stats["Number of unique CSMs"] == pytest.approx(826.0)
    assert stats["Number of intra CSMs"] == pytest.approx(803.0)
    assert stats["Number of inter CSMs"] == pytest.approx(23.0)
    assert stats["Number of target-target CSMs"] == pytest.approx(786.0)
    assert stats["Number of target-decoy CSMs"] == pytest.approx(39.0)
    assert stats["Number of decoy-decoy CSMs"] == pytest.approx(1.0)
    assert stats["Minimum CSM score"] == pytest.approx(1.11)
    assert stats["Maximum CSM score"] == pytest.approx(452.99)
    assert stats["Number of crosslinks"] == pytest.approx(300.0)
    assert stats["Number of unique crosslinks by peptide"] == pytest.approx(300.0)
    assert stats["Number of unique crosslinks by protein"] == pytest.approx(298.0)
    assert stats["Number of intra crosslinks"] == pytest.approx(279.0)
    assert stats["Number of inter crosslinks"] == pytest.approx(21.0)
    assert stats["Number of target-target crosslinks"] == pytest.approx(265.0)
    assert stats["Number of target-decoy crosslinks"] == pytest.approx(0.0)
    assert stats["Number of decoy-decoy crosslinks"] == pytest.approx(35.0)
    assert stats["Minimum crosslink score"] == pytest.approx(1.11)
    assert stats["Maximum crosslink score"] == pytest.approx(452.99)
