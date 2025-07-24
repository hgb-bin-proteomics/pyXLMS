#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test1():
    from pyXLMS import parser
    from pyXLMS import plotting

    pr = parser.read_msannika(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
    )
    csms = pr["crosslink-spectrum-matches"]
    fig, ax = plotting.plot_crosslink_type_distribution(csms)
    assert fig is not None
    assert ax is not None


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test2():
    from pyXLMS import parser
    from pyXLMS import plotting

    pr = parser.read_msannika(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx"
    )
    crosslinks = pr["crosslinks"]
    fig, ax = plotting.plot_crosslink_type_distribution(crosslinks)
    assert fig is not None
    assert ax is not None


def test3():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_crosslink_type_distribution

    pr = parser.read_msannika(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
    )
    csms = pr["crosslink-spectrum-matches"]

    with pytest.raises(
        IndexError,
        match=r"At least two colors need to be given for the plot!",
    ):
        _plot = plot_crosslink_type_distribution(csms, colors=["#6d4bff"])


def test4():
    from pyXLMS.plotting import plot_crosslink_type_distribution

    with pytest.raises(
        ValueError,
        match=r"Can't plot crosslink type distribution if no crosslink-spectrum-matches or crosslinks are given!",
    ):
        _plot = plot_crosslink_type_distribution([])


def test5():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_crosslink_type_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    pr["crosslink-spectrum-matches"][0]["data_type"] = "peptide-spectrum-match"

    with pytest.raises(
        TypeError,
        match=r"Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!",
    ):
        _plot = plot_crosslink_type_distribution(pr["crosslink-spectrum-matches"])
