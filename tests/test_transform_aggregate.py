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


def test12():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read(
        "data/_test/aggregate/csms_min.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_peptide = aggregate(pr["crosslink-spectrum-matches"])
    assert len(aggregate_peptide) == 3
    assert [xl["alpha_peptide"] for xl in aggregate_peptide] == [
        "KPEPTIDE",
        "PEKPTIDE",
        "PEPKTIDE",
    ]


def test13():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read(
        "data/_test/aggregate/csms_min.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = (
        r"Grouping by protein crosslink position is only available if all crosslink-spectrum-matches have defined protein crosslink positions!\n"
        r"This error might be fixable with 'transform\.reannotate_positions\(\)'\!"
    )
    with pytest.raises(ValueError, match=err_str):
        _aggregate_protein = aggregate(pr["crosslink-spectrum-matches"], by="protein")


def test14():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/csms_score.txt",
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
        ["PROTB", "PROTA"],
    ]
    assert u[0]["score"] == pytest.approx(2.5)
    assert u[1]["score"] == pytest.approx(4.5)
    assert u[2]["score"] == pytest.approx(6.5)
    assert u[3]["score"] == pytest.approx(8.5)
    assert u[4]["score"] == pytest.approx(10.5)


def test15():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    u = unique(pr["crosslink-spectrum-matches"], score="lower_better")
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
    assert u[0]["score"] == pytest.approx(1.5)
    assert u[1]["score"] == pytest.approx(3.5)
    assert u[2]["score"] == pytest.approx(5.5)
    assert u[3]["score"] == pytest.approx(7.5)
    assert u[4]["score"] == pytest.approx(9.5)


def test16():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/xls_score.txt",
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
        ["PROTB", "PROTA"],
    ]
    assert u[0]["score"] == pytest.approx(4.5)
    assert u[1]["score"] == pytest.approx(8.5)
    assert u[2]["score"] == pytest.approx(10.5)


def test17():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/xls_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslinks"]) == 10
    u = unique(pr["crosslinks"], score="lower_better")
    assert len(u) == 3
    assert [xl["alpha_peptide"] for xl in u] == ["KPEPTIDE", "PEKPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in u] == [
        ["PROTA"],
        ["PROTA"],
        ["PROTA", "PROTB"],
    ]
    assert u[0]["score"] == pytest.approx(1.5)
    assert u[1]["score"] == pytest.approx(5.5)
    assert u[2]["score"] == pytest.approx(9.5)


def test18():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/xls_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslinks"]) == 10
    u = unique(pr["crosslinks"], by="protein")
    assert len(u) == 2
    assert [xl["alpha_peptide"] for xl in u] == ["PEKPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in u] == [["PROTA"], ["PROTB", "PROTA"]]
    assert u[0]["score"] == pytest.approx(8.5)
    assert u[1]["score"] == pytest.approx(10.5)


def test19():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/xls_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslinks"]) == 10
    u = unique(pr["crosslinks"], by="protein", score="lower_better")
    assert len(u) == 2
    assert [xl["alpha_peptide"] for xl in u] == ["KPEPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in u] == [["PROTA"], ["PROTA", "PROTB"]]
    assert u[0]["score"] == pytest.approx(1.5)
    assert u[1]["score"] == pytest.approx(9.5)


def test20():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_peptide = aggregate(pr["crosslink-spectrum-matches"])
    assert len(aggregate_peptide) == 3
    assert [xl["alpha_peptide"] for xl in aggregate_peptide] == [
        "KPEPTIDE",
        "PEKPTIDE",
        "PEPKTIDE",
    ]
    assert [xl["alpha_proteins"] for xl in aggregate_peptide] == [
        ["PROTA"],
        ["PROTA"],
        ["PROTB", "PROTA"],
    ]
    assert aggregate_peptide[0]["score"] == pytest.approx(4.5)
    assert aggregate_peptide[1]["score"] == pytest.approx(8.5)
    assert aggregate_peptide[2]["score"] == pytest.approx(10.5)


def test21():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_peptide = aggregate(
        pr["crosslink-spectrum-matches"], score="lower_better"
    )
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
    assert aggregate_peptide[0]["score"] == pytest.approx(1.5)
    assert aggregate_peptide[1]["score"] == pytest.approx(5.5)
    assert aggregate_peptide[2]["score"] == pytest.approx(9.5)


def test22():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_protein = aggregate(pr["crosslink-spectrum-matches"], by="protein")
    assert len(aggregate_protein) == 2
    assert [xl["alpha_peptide"] for xl in aggregate_protein] == ["PEKPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in aggregate_protein] == [
        ["PROTA"],
        ["PROTB", "PROTA"],
    ]
    assert aggregate_protein[0]["score"] == pytest.approx(8.5)
    assert aggregate_protein[1]["score"] == pytest.approx(10.5)


def test23():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    aggregate_protein = aggregate(
        pr["crosslink-spectrum-matches"], by="protein", score="lower_better"
    )
    assert len(aggregate_protein) == 2
    assert [xl["alpha_peptide"] for xl in aggregate_protein] == ["KPEPTIDE", "PEPKTIDE"]
    assert [xl["alpha_proteins"] for xl in aggregate_protein] == [
        ["PROTA"],
        ["PROTA", "PROTB"],
    ]
    assert aggregate_protein[0]["score"] == pytest.approx(1.5)
    assert aggregate_protein[1]["score"] == pytest.approx(9.5)


def test24():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = (
        r"Parameter 'by' has to be one of 'peptide' or 'protein'! Option 'peptide' will group by peptide sequence and "
        r"peptide crosslink position while option 'protein' will group by protein identifier and protein crosslink position."
    )
    with pytest.raises(TypeError, match=err_str):
        _u = unique(pr["crosslink-spectrum-matches"], by="sequence")


def test25():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = (
        r"Parameter 'by' has to be one of 'peptide' or 'protein'! Option 'peptide' will group by peptide sequence and "
        r"peptide crosslink position while option 'protein' will group by protein identifier and protein crosslink position."
    )
    with pytest.raises(TypeError, match=err_str):
        _a = aggregate(pr["crosslink-spectrum-matches"], by="sequence")


def test26():
    from pyXLMS.parser import read
    from pyXLMS.transform import unique

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = (
        r"Parameter 'score' has to be one of 'higher_better' or 'lower_better'! If two identical crosslinks or crosslink-spectrum"
        r"-matches are found, the one with the higher score is kept if 'higher_better' is selected, and vice versa."
    )
    with pytest.raises(TypeError, match=err_str):
        _u = unique(pr["crosslink-spectrum-matches"], score="lower")


def test27():
    from pyXLMS.parser import read
    from pyXLMS.transform import aggregate

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = (
        r"Parameter 'score' has to be one of 'higher_better' or 'lower_better'! If two identical crosslinks or crosslink-spectrum"
        r"-matches are found, the one with the higher score is kept if 'higher_better' is selected, and vice versa."
    )
    with pytest.raises(TypeError, match=err_str):
        _a = aggregate(pr["crosslink-spectrum-matches"], score="lower")
