#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.exporter import to_xlinkdb
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    crosslinks = pr["crosslinks"]
    df = to_xlinkdb(crosslinks, filename="crosslinksForXlinkDB")
    assert df.shape[0] == 227
    assert df.shape[1] == 7


def test2():
    from pyXLMS.exporter import to_xlinkdb
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    crosslinks = pr["crosslinks"]
    df = to_xlinkdb(crosslinks, filename=None)
    assert df.shape[0] == 227
    assert df.shape[1] == 7


def test3():
    from pyXLMS.exporter import to_xlinkdb
    from pyXLMS.parser import read

    pr = read(
        "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv",
        engine="xiSearch/xiFDR",
        crosslinker="DSS",
    )
    crosslinks = pr["crosslinks"]

    with pytest.raises(
        ValueError,
        match="Parameter filename must only contain alpha-numeric characters and no file extension!",
    ):
        _df = to_xlinkdb(crosslinks, filename="crosslinksForXlinkDB.tsv")


def test4():
    from pyXLMS.exporter import to_xlinkdb

    with pytest.raises(
        ValueError,
        match="Provided crosslinks contain no elements!",
    ):
        _df = to_xlinkdb([], filename="crosslinksForXlinkDB")


def test5():
    from pyXLMS.exporter import to_xlinkdb
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]

    with pytest.raises(
        RuntimeError,
        match="Can't export to XlinkDB because not all necessary information is available!",
    ):
        _df = to_xlinkdb(crosslinks, filename=None)
