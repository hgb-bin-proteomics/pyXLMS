#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.pipelines import pipeline

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
        unique=True,
        validate={"fdr": 0.05, "formula": "(TD-DD)/TT"},
        targets_only=True,
    )
    assert pr["crosslink-spectrum-matches"] is not None
    assert len(pr["crosslink-spectrum-matches"]) == 786


def test2():
    from pyXLMS.pipelines import pipeline

    with pytest.raises(
        TypeError,
        match=r"Parameter unique has to be a dictionary of parameters for transform\.unique\(\), a boolean or None!",
    ):
        _pr = pipeline(
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
            engine="MS Annika",
            crosslinker="DSS",
            unique=1,
            validate=True,
            targets_only=True,
        )


def test3():
    from pyXLMS.pipelines import pipeline

    with pytest.raises(
        TypeError,
        match=r"Parameter validate has to be a dictionary of parameters for transform\.validate\(\), a boolean or None!",
    ):
        _pr = pipeline(
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
            engine="MS Annika",
            crosslinker="DSS",
            unique=True,
            validate=1,
            targets_only=True,
        )


def test4():
    from pyXLMS.pipelines import pipeline

    with pytest.raises(
        TypeError,
        match=r"Parameter targets\_only has to be a boolean or None!",
    ):
        _pr = pipeline(
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
            engine="MS Annika",
            crosslinker="DSS",
            unique=True,
            validate=True,
            targets_only=1,
        )
