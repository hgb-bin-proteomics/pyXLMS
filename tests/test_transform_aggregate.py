#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        ["data/_test/aggregate/csms.txt", "data/_test/aggregate/xls.txt"],
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    assert len(pr["crosslinks"]) == 10
    unique_peptide = unique(pr, by="peptide")
    assert len(unique_peptide["crosslink-spectrum-matches"]) == 5
    assert len(unique_peptide["crosslinks"]) == 3


def test2():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        ["data/_test/aggregate/csms.txt", "data/_test/aggregate/xls.txt"],
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    assert len(pr["crosslinks"]) == 10
    unique_protein = unique(pr, by="protein")
    assert len(unique_protein["crosslink-spectrum-matches"]) == 5
    assert len(unique_protein["crosslinks"]) == 2


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read("data/_test/aggregate/csms.txt", engine="custom", crosslinker="DSS")
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_peptide = aggregate(pr["crosslink-spectrum-matches"], by="peptide")
    assert len(aggregate_peptide) == 3


def test4():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read("data/_test/aggregate/csms.txt", engine="custom", crosslinker="DSS")
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_protein = aggregate(pr["crosslink-spectrum-matches"], by="protein")
    assert len(aggregate_protein) == 2


def test5():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/csms.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    u = unique(pr["crosslink-spectrum-matches"])
    assert len(u) == 5
    assert [csm["alpha_peptide"] for csm in u] == [
        "KPEPTIDE",
        "KPEPTIDE",
        "PEKPTIDE",
        "PEKPTIDE",
        "PEPKTIDE",
    ]
    assert [csm["alpha_proteins"] for csm in u] == [
        ["PROTA"],
        ["PROTA"],
        ["PROTA"],
        ["PROTA"],
        ["PROTA", "PROTB"],
    ]


def test6():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/xls.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslinks"]) == 10
    u = unique(pr["crosslinks"])
    assert len(u) == 3
    assert [xl["alpha_peptide"] for xl in u] == ["KPEPTIDE", "PEKPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in u] == [
        ["PROTA"],
        ["PROTA"],
        ["PROTA", "PROTB"],
    ]


def test7():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/xls.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslinks"]) == 10
    u = unique(pr["crosslinks"], by="protein")
    assert len(u) == 2
    assert [xl["alpha_peptide"] for xl in u] == ["KPEPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in u] == [["PROTA"], ["PROTA", "PROTB"]]


def test8():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read("data/_test/aggregate/csms.txt", engine="custom", crosslinker="DSS")
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_peptide = aggregate(pr["crosslink-spectrum-matches"], by="peptide")
    assert len(aggregate_peptide) == 3
    assert [xl["alpha_peptide"] for xl in aggregate_peptide] == [
        "KPEPTIDE",
        "PEKPTIDE",
        "PEPKTIDE",
    ]
    assert [xl["alpha_proteins"] for xl in aggregate_peptide] == [
        ["PROTA"],
        ["PROTA"],
        ["PROTA", "PROTB"],
    ]


def test9():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read("data/_test/aggregate/csms.txt", engine="custom", crosslinker="DSS")
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_protein = aggregate(pr["crosslink-spectrum-matches"], by="protein")
    assert len(aggregate_protein) == 2
    assert [xl["alpha_peptide"] for xl in aggregate_protein] == ["KPEPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in aggregate_protein] == [
        ["PROTA"],
        ["PROTA", "PROTB"],
    ]


def test10():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/xls_min.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslinks"]) == 10
    u = unique(pr["crosslinks"])
    assert len(u) == 3
    assert [xl["alpha_peptide"] for xl in u] == ["KPEPTIDE", "PEKPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in u] == [
        ["PROTA"],
        ["PROTA"],
        ["PROTA", "PROTB"],
    ]


def test11():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/xls_min.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslinks"]) == 10
    err_str = (
        r"Grouping by protein crosslink position is only available if all crosslinks have defined protein crosslink positions!\n"
        r"This error might be fixable with 'transform\.reannotate_positions\(\)'\!"
    )
    with pytest.raises(ValueError, match=err_str):
        _u = unique(pr["crosslinks"], by="protein")
