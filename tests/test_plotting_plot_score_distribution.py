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
    fig, ax = plotting.plot_score_distribution(csms)
    assert fig is not None
    assert ax is not None


def test2():
    from pyXLMS.plotting import plot_score_distribution

    with pytest.raises(
        ValueError,
        match=r"Can't plot score distribution if no crosslink-spectrum-matches or crosslinks are given!",
    ):
        _plot = plot_score_distribution([])


def test3():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_score_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    with pytest.raises(
        TypeError,
        match=r"Provided input is not a valid list of type CrosslinkSpectrumMatch or Crosslink!",
    ):
        _plot = plot_score_distribution([pr])


def test4():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_score_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    csms = pr.csms()
    csms[0] = csms[0].copy_with_update({"score": None})

    with pytest.raises(
        ValueError,
        match="Attribute .* is missing in at least one element but is required!",
    ):
        _plot = plot_score_distribution(csms)


def test5():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_score_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    csms = pr.csms()
    csms[0] = csms[0].copy_with_update({"alpha_decoy": None})

    with pytest.raises(
        ValueError,
        match="Attribute .* is missing in at least one element but is required!",
    ):
        _plot = plot_score_distribution(csms)


def test6():
    from pyXLMS.parser import read
    from pyXLMS.plotting import plot_score_distribution

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    csms = pr.csms()
    csms[0] = csms[0].copy_with_update({"beta_decoy": None})

    with pytest.raises(
        ValueError,
        match="Attribute .* is missing in at least one element but is required!",
    ):
        _plot = plot_score_distribution(csms)
