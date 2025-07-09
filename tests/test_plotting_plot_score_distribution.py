#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS import parser
    from pyXLMS import plotting
    pr = parser.read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")
    csms = pr["crosslink-spectrum-matches"]
    fig, ax = plotting.plot_score_distribution(csms)
    assert fig is not None
    assert ax is not None


def test2():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_score_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    with pytest.raises(
        ValueError,
        match=r"FDR must be given as a real number between 0 and 1, e\.g\. 0\.01 corresponds to 1\% FDR!",
    ):
        _plot = plot_score_distribution([])


def test4():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_score_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    with pytest.raises(
        ValueError,
        match=r"FDR must be given as a real number between 0 and 1, e\.g\. 0\.01 corresponds to 1\% FDR!",
    ):
        _plot = plot_score_distribution()
