#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.exporter import to_impxfdr
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    crosslinks = pr["crosslinks"]
    df = to_impxfdr(crosslinks, filename="crosslinks.xlsx")
    assert df.shape[0] == 225
    assert df.shape[1] == 11


def test2():
    from pyXLMS.exporter import to_impxfdr
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    crosslinks = pr["crosslinks"]
    df = to_impxfdr(crosslinks, filename=None)
    assert df.shape[0] == 225
    assert df.shape[1] == 11


def test3():
    from pyXLMS.exporter import to_impxfdr
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    df = to_impxfdr(csms, filename="csms.xlsx")
    assert df.shape[0] == 411
    assert df.shape[1] == 11


def test4():
    from pyXLMS.exporter import to_impxfdr
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    df = to_impxfdr(csms, filename=None, targets_only=False)
    assert df.shape[0] == 413
    assert df.shape[1] == 11


def test5():
    from pyXLMS.exporter import to_impxfdr
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]

    with pytest.raises(
        ValueError,
        match="Provided data does not contain any crosslinks or crosslink-spectrum-matches after filtering for targets only!",
    ):
        _df = to_impxfdr(crosslinks, filename=None)


def test6():
    from pyXLMS.exporter import to_impxfdr
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]

    with pytest.raises(
        RuntimeError,
        match="Can't export to IMP-X-FDR because not all necessary information is available!",
    ):
        _df = to_impxfdr(crosslinks, filename=None, targets_only=False)
