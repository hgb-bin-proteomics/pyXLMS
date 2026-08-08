#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import io

import pytest


def test1():
    # A xiFDR "Links" (crosslinks) export whose "PSMIDs" column holds numeric PSM
    # identifiers - as emitted by e.g. xiFDR 2.1.5.2 - carries no peptide sequences
    # to recover, so crosslink-level reading has to fail with a clear RuntimeError
    # that points to the CSM export instead of an opaque IndexError.
    from pyXLMS.parser import read_xi

    links = (
        "PSMIDs,Protein1,Protein2,fromSite,ToSite,Decoy1,Decoy2,Score\n"
        "123456,Cas9,Cas9,10,20,false,false,15.5\n"
    )
    with pytest.raises(RuntimeError, match="the 'PSMIDs' column is not in the expected"):
        _r = read_xi(io.StringIO(links), decoy_prefix="rev_", verbose=0)


def test2():
    # A xiFDR "Links" export in the expected "P1_<peptide> P2_<peptide> <pos1>
    # <pos2>" PSMIDs format still reads normally - the version guard must not reject
    # valid inputs.
    from pyXLMS import parser as p

    links = (
        "PSMIDs,Protein1,Protein2,fromSite,ToSite,Decoy1,Decoy2,Score\n"
        "P1_PEPTIDEK P2_KLESIER 4 2,Cas9,Cas9,10,20,false,false,15.5\n"
    )
    pr = p.read_xi(io.StringIO(links), decoy_prefix="rev_", verbose=0)

    assert pr["data_type"] == "parser_result"
    assert pr["crosslinks"] is not None
    assert pr["crosslink-spectrum-matches"] is None

    xls = pr["crosslinks"]
    assert len(xls) == 1
    xl = xls[0]
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["beta_proteins"] == ["Cas9"]
    # pyXLMS >=2.0 canonicalizes a symmetric crosslink's alpha/beta order, so check
    # the two endpoints as an unordered set; each peptide keeps its protein position.
    endpoints = {
        (xl["alpha_peptide"], xl["alpha_proteins_crosslink_positions"][0]),
        (xl["beta_peptide"], xl["beta_proteins_crosslink_positions"][0]),
    }
    assert endpoints == {("PEPTIDEK", 10), ("KLESIER", 20)}


def test3():
    # The scan column is lowercase "scan" in newer xiFDR versions (2.1.5.2, 2.2.1)
    # and uppercase "Scan" in older ones. __get_xifdr_scan accepts either spelling
    # so that a diagnostic message can never itself raise a KeyError while trying to
    # report a different error; missing both spellings is a clear KeyError. The
    # helper is internal and only reached on the modification-diagnostic path, so it
    # is exercised here directly.
    import pandas as pd

    from pyXLMS.parser import _parser_xldbse_xi as xi

    get_scan = getattr(xi, "__get_xifdr_scan")

    assert get_scan(pd.Series({"scan": 4321})) == 4321
    assert get_scan(pd.Series({"Scan": 1234})) == 1234
    # lowercase is preferred when both are present, matching __read_xifdr_csms
    assert get_scan(pd.Series({"scan": 1, "Scan": 2})) == 1

    with pytest.raises(KeyError, match="Neither 'scan' nor 'Scan'"):
        _r = get_scan(pd.Series({"ScanId": 9}))
