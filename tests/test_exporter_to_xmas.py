#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.exporter import to_xmas

    with pytest.raises(
        ValueError,
        match="Provided crosslinks contain no elements!",
    ):
        _df = to_xmas([], filename=None)


def test2():
    from pyXLMS.exporter import to_xmas
    from pyXLMS.data import create_csm_min
    from pyXLMS.data import create_crosslink_min

    csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    xl1 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    mixed = [xl1, csm1]
    with pytest.raises(
        TypeError, match="Not all elements in data have the same data type!"
    ):
        _df = to_xmas(mixed, filename=None)


def test3():
    from pyXLMS.exporter import to_xmas
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]
    df = to_xmas(crosslinks, filename="crosslinks_xmas.xlsx")
    assert df.shape[0] == 2
    assert df.shape[1] == 2


def test4():
    from pyXLMS.exporter import to_xmas
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]
    df = to_xmas(crosslinks, filename=None)
    assert df.shape[0] == 2
    assert df.shape[1] == 2


def test5():
    from pyXLMS.exporter import to_xmas
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    crosslinks = pr["crosslinks"]
    df = to_xmas(crosslinks, filename=None)
    assert df.shape[0] == 227
    assert df.shape[1] == 2


def test6():
    from pyXLMS.exporter import to_xmas
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    with pytest.raises(
        TypeError,
        match="Unsupported data type for input crosslinks! Parameter crosslinks has to be a list of crosslinks!",
    ):
        _df = to_xmas(csms, filename=None)
