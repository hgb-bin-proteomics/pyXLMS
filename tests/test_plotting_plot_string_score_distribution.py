#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
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


@pytest.mark.external
@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test1():
    from pyXLMS import parser
    from pyXLMS import plotting

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    fig, ax = plotting.plot_string_score_distribution(xls, organism="Homo sapiens")
    assert fig is not None
    assert ax is not None


@pytest.mark.external
@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test2():
    from pyXLMS import parser
    from pyXLMS import plotting

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    fig, ax = plotting.plot_string_score_distribution(
        xls, organism="Homo sapiens", plot_type="bar"
    )
    assert fig is not None
    assert ax is not None


@pytest.mark.external
@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test3():
    from pyXLMS import parser
    from pyXLMS import plotting

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    fig, ax = plotting.plot_string_score_distribution(
        xls, organism=9606, plot_type="hist"
    )
    assert fig is not None
    assert ax is not None


def test4():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_string_score_distribution

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    with pytest.raises(
        TypeError,
        match=r"Plot type has to be one of 'bar', or 'hist'!",
    ):
        _plot = plot_string_score_distribution(xls, plot_type="pie")


def test5():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_string_score_distribution

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    with pytest.raises(
        ValueError,
        match=r"Six colors are required for plot type 'bar'!",
    ):
        _plot = plot_string_score_distribution(
            xls, plot_type="bar", colors=["red", "green", "blue"]
        )


def test6():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_string_score_distribution

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    with pytest.raises(
        ValueError,
        match=r"Three colors are required for plot type 'hist'!",
    ):
        _plot = plot_string_score_distribution(
            xls, plot_type="hist", colors=["red", "green"]
        )


def test7():
    from pyXLMS import parser
    from pyXLMS.transform import filter_crosslink_type
    from pyXLMS.plotting import plot_string_score_distribution

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    intra = filter_crosslink_type(xls)["Intra"]
    with pytest.raises(
        ValueError,
        match=r"Can't plot STRING score distribution because data does not contain inter-links!",
    ):
        _plot = plot_string_score_distribution(intra)


def test8():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_string_score_distribution

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    with pytest.raises(
        ValueError,
        match=r"Input data does not have annotated STRING scores! In this case a valid organism has to be given!",
    ):
        _plot = plot_string_score_distribution(xls)


def test9():
    from pyXLMS.plotting import plot_string_score_distribution

    with pytest.raises(
        ValueError,
        match=r"Can't plot STRING score distribution if no crosslink-spectrum-matches or crosslinks are given!",
    ):
        _plot = plot_string_score_distribution([])


def test10():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_string_score_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    with pytest.raises(
        TypeError,
        match=r"Provided input is not a valid list of type CrosslinkSpectrumMatch or Crosslink!",
    ):
        _plot = plot_string_score_distribution([pr])
