#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    # cross-linked: two peptides; proteins de-duplicated and sorted per site
    from pyXLMS.parser import parse_sites_from_plink

    sites = parse_sites_from_plink(
        "KALAAAGYDVEK(1)-AVAASKER(6)",
        "PROTA (64)-PROTB (52)/PROTA (64)-PROTC (55)/",
    )
    assert len(sites) == 2
    assert sites[0]["peptide"] == "KALAAAGYDVEK"
    assert sites[0]["peptide_position"] == 1
    assert sites[0]["proteins"] == ["PROTA"]
    assert sites[0]["protein_positions"] == [64]
    assert sites[1]["peptide"] == "AVAASKER"
    assert sites[1]["peptide_position"] == 6
    assert sites[1]["proteins"] == ["PROTB", "PROTC"]
    assert sites[1]["protein_positions"] == [52, 55]


def test2():
    # a non-cross-linked (mono-/loop-linked) sequence is rejected with a clear error
    from pyXLMS.parser import parse_sites_from_plink

    with pytest.raises(ValueError, match="Not a cross-linked pLink sequence"):
        _r = parse_sites_from_plink("APK(3)", "FUS (451)/")


def test3():
    from pyXLMS.parser import parse_sites_from_plink

    site_a, site_b = parse_sites_from_plink("PEPK(4)-KELS(1)", "PROT (10)-PROT (20)/")
    assert site_a["peptide"] == "PEPK"


def test4():
    from pyXLMS.parser import parse_sites_from_plink

    with pytest.raises(ValueError, match="Not a cross-linked pLink sequence"):
        _r = parse_sites_from_plink("APK(3)-PEPK(4)", "FUS (451)/")


def test5():
    from pyXLMS.parser import parse_sites_from_plink

    with pytest.raises(ValueError, match="Parsed more than two possible sites"):
        _r = parse_sites_from_plink("APK(3)-PEPK(4)-TIDE(4)", "FUS (451)-BAZ (31)/")


def test6():
    from pyXLMS.parser import parse_sites_from_plink

    with pytest.raises(ValueError, match="less or more than two possible sites"):
        _r = parse_sites_from_plink("APK(3)-PEPK(4)", "FUS (451)-BAZ (31)-BAR (32)/")


def test7():
    # isoform accessions contain a hyphen (e.g. "sp|Q14240-2|..."); the crosslink
    # is still two-sided. The site delimiter is the "-" that follows a position, so
    # the accession is not split mid-string and the ambiguous group is preserved.
    from pyXLMS.parser import parse_sites_from_plink

    site_a, site_b = parse_sites_from_plink(
        "ELAQQIQK(8)-AEVQK(5)",
        "sp|Q14240|IF4A2 (119)-sp|Q14240|IF4A2 (147)/"
        "sp|Q14240|IF4A2 (119)-sp|Q14240-2|IF4A2 (148)/"
        "sp|Q14240-2|IF4A2 (120)-sp|Q14240|IF4A2 (147)/"
        "sp|Q14240-2|IF4A2 (120)-sp|Q14240-2|IF4A2 (148)/",
    )
    assert site_a["peptide"] == "ELAQQIQK"
    assert site_a["proteins"] == ["sp|Q14240-2|IF4A2", "sp|Q14240|IF4A2"]
    assert site_a["protein_positions"] == [120, 119]
    assert site_b["peptide"] == "AEVQK"
    assert site_b["proteins"] == ["sp|Q14240-2|IF4A2", "sp|Q14240|IF4A2"]
    assert site_b["protein_positions"] == [148, 147]
