#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

XIFDR_COLS = [
    "run",
    "scan",
    "peptide1",
    "peptide2",
    "peptide link 1",
    "peptide link 2",
    "is decoy 1",
    "is decoy 2",
    "precursor charge",
    "accession1",
    "accession2",
    "peptide position 1",
    "peptide position 2",
    "score",
]


def test1():
    from pyXLMS.exporter import to_xifdr
    from pyXLMS.parser import read

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    df = to_xifdr(csms, filename="msannika_xiFDR.csv")
    assert df.shape[0] == 826
    assert df.shape[1] == 14
    cols = df.columns.values.tolist()
    for col in XIFDR_COLS:
        assert col in cols


def test2():
    from pyXLMS.exporter import to_xifdr
    from pyXLMS.parser import read

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    df = to_xifdr(csms, filename=None)
    assert df.shape[0] == 826
    assert df.shape[1] == 14
    cols = df.columns.values.tolist()
    for col in XIFDR_COLS:
        assert col in cols


def test3():
    from pyXLMS.exporter import to_xifdr

    with pytest.raises(
        ValueError,
        match="Provided crosslink-spectrum-matches contain no elements!",
    ):
        _df = to_xifdr([], filename=None)


def test4():
    from pyXLMS.exporter import to_xifdr
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]

    with pytest.raises(
        TypeError,
        match="Unsupported data type for input csms! Parameter csms has to be a list of crosslink-spectrum-matches!",
    ):
        _df = to_xifdr(crosslinks, filename=None)


def test5():
    from pyXLMS.exporter import to_xifdr
    from pyXLMS.data import create_csm_min

    csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    csm2 = create_csm_min("PEKPTIDE", 3, "PEPKTIDE", 4, "RUN_1", 2)
    csms = [csm1, csm2]

    with pytest.raises(
        RuntimeError,
        match="Can't export to xiFDR because not all necessary information is available!",
    ):
        _df = to_xifdr(csms, filename=None)


def test6():
    from pyXLMS.exporter import to_xifdr
    from pyXLMS.parser import read

    pr = read("data/scout/Cas9_Unfiltered_CSMs.csv", engine="Scout", crosslinker="DSSO")
    csms = pr["crosslink-spectrum-matches"]

    with pytest.raises(
        RuntimeError,
        match="Can't export to xiFDR because not all necessary information is available!",
    ):
        _df = to_xifdr(csms, filename=None)
