#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

XINET_COLS = [
    "Protein1",
    "PepPos1",
    "PepSeq1",
    "LinkPos1",
    "Protein2",
    "PepPos2",
    "PepSeq2",
    "LinkPos2",
    "Id",
]


def test1():
    from pyXLMS.exporter import to_xinet
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only
    from pyXLMS.transform import filter_proteins

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    crosslinks = targets_only(pr)["crosslinks"]
    cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    df = to_xinet(cas9, filename="crosslinks_xiNET.csv")
    assert df.shape[0] == 253
    assert df.shape[1] == 10


def test2():
    from pyXLMS.exporter import to_xinet
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only
    from pyXLMS.transform import filter_proteins

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    crosslinks = targets_only(pr)["crosslinks"]
    cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    df = to_xinet(cas9, filename=None)
    assert df.shape[0] == 253
    assert df.shape[1] == 10


def test3():
    from pyXLMS.exporter import to_xinet
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only
    from pyXLMS.transform import filter_proteins

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    crosslinks = targets_only(pr)["crosslinks"]
    cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    df = to_xinet(cas9, filename=None)
    cols = df.columns.values.tolist()
    for col in XINET_COLS:
        assert col in cols
    assert "Score" in cols


def test4():
    from pyXLMS.exporter import to_xinet
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only
    from pyXLMS.transform import filter_proteins

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    crosslinks = targets_only(pr)["crosslinks"]
    cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
    for xl in cas9:
        xl["score"] = None
    df = to_xinet(cas9, filename=None)
    cols = df.columns.values.tolist()
    for col in XINET_COLS:
        assert col in cols
    assert "Score" not in cols


def test5():
    from pyXLMS.exporter import to_xinet

    with pytest.raises(
        ValueError,
        match="Provided crosslinks contain no elements!",
    ):
        _df = to_xinet([], filename="crosslinks_xiNET.csv")


def test6():
    from pyXLMS.exporter import to_xinet
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]

    with pytest.raises(
        RuntimeError,
        match="Can't export to xiNET because not all necessary information is available!",
    ):
        _df = to_xinet(crosslinks, filename=None)
