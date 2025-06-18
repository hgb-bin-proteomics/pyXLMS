#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    targets = targets_only(result["crosslink-spectrum-matches"])
    assert len(targets) == 786


def test2():
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    targets = targets_only(result["crosslinks"])
    len(targets) == 265


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only

    result = read(
        [
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        ],
        engine="MS Annika",
        crosslinker="DSS",
    )
    result_targets = targets_only(result)
    assert len(result_targets["crosslink-spectrum-matches"]) == 786
    assert len(result_targets["crosslinks"]) == 265


def test4():
    from pyXLMS.transform import targets_only
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]
    with pytest.raises(
        RuntimeError,
        match=r"No target matches found! Are you sure your data is labelled\?",
    ):
        _r = targets_only(crosslinks)


def test5():
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    result["crosslink-spectrum-matches"] = []
    with pytest.raises(
        RuntimeError,
        match=r"No target crosslink-spectrum-matches found! Are you sure they are labelled\?",
    ):
        _r = targets_only(result)


def test6():
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    result["crosslinks"] = []
    with pytest.raises(
        RuntimeError,
        match=r"No target crosslinks found! Are you sure they are labelled\?",
    ):
        _r = targets_only(result)
