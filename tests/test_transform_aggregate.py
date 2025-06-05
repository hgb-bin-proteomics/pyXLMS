#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


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
