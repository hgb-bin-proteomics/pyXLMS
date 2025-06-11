#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


F1 = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.mzid"
F2 = "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3.mzid"


def test1():
    from pyXLMS.parser import parse_scan_nr_from_mzid

    assert parse_scan_nr_from_mzid("scan=5321") == 5321


@pytest.mark.slow
def test2():
    from pyXLMS import parser as p

    pr = p.read_mzid(F1, verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "mzIdentML"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 786


@pytest.mark.slow
def test3():
    from pyXLMS import parser as p

    pr = p.read_mzid(F2, verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "mzIdentML"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 823
