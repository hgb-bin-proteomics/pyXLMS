#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


@pytest.fixture(autouse=True)
def cleanup_figures():
    from matplotlib import use
    import matplotlib.pyplot as plt

    use("agg")

    yield

    plt.close(fig="all")


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test1():
    from pyXLMS import parser
    from pyXLMS import plotting

    pr = parser.read_msannika(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
    )
    csms = pr["crosslink-spectrum-matches"]
    fig, ax = plotting.plot_peptide_pair_distribution(csms)
    assert fig is not None
    assert ax is not None


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test2():
    from pyXLMS import parser
    from pyXLMS import plotting

    pr = parser.read_msannika(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
    )
    csms = pr["crosslink-spectrum-matches"]
    fig, ax = plotting.plot_peptide_pair_distribution(csms, top_n=2)
    assert fig is not None
    assert ax is not None


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test3():
    from pyXLMS import parser
    from pyXLMS import plotting

    pr = parser.read_msannika(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
    )
    csms = pr["crosslink-spectrum-matches"]
    fig, ax = plotting.plot_peptide_pair_distribution(csms, top_n=10000)
    assert fig is not None
    assert ax is not None


def test4():
    from pyXLMS.plotting import plot_peptide_pair_distribution

    with pytest.raises(
        ValueError,
        match=r"Can't plot peptide pair distribution if no crosslink-spectrum-matches are given!",
    ):
        _plot = plot_peptide_pair_distribution([])


def test5():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_peptide_pair_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    pr["crosslink-spectrum-matches"][0]["data_type"] = "peptide-spectrum-match"

    with pytest.raises(
        TypeError,
        match=r"Unsupported data type for input data! Parameter data has to be a list of crosslink-spectrum-match!",
    ):
        _plot = plot_peptide_pair_distribution(pr["crosslink-spectrum-matches"])
