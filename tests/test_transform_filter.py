#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_target_decoy

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    target_and_decoys = filter_target_decoy(result["crosslink-spectrum-matches"])
    assert len(target_and_decoys["Target-Target"]) == 786
    assert len(target_and_decoys["Target-Decoy"]) == 39
    assert len(target_and_decoys["Decoy-Decoy"]) == 1


def test2():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_target_decoy

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    target_and_decoys = filter_target_decoy(result["crosslinks"])
    assert len(target_and_decoys["Target-Target"]) == 265
    assert len(target_and_decoys["Target-Decoy"]) == 0
    assert len(target_and_decoys["Decoy-Decoy"]) == 35


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_proteins

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    proteins_csms = filter_proteins(result["crosslink-spectrum-matches"], ["Cas9"])
    assert proteins_csms["Proteins"] == ["Cas9"]
    assert len(proteins_csms["Both"]) == 798
    assert len(proteins_csms["One"]) == 23


def test4():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_proteins

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    proteins_xls = filter_proteins(result["crosslinks"], ["Cas9"])
    assert proteins_xls["Proteins"] == ["Cas9"]
    assert len(proteins_xls["Both"]) == 274
    assert len(proteins_xls["One"]) == 21


def test5():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_crosslink_type

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    crosslink_type_filtered_csms = filter_crosslink_type(
        result["crosslink-spectrum-matches"]
    )
    assert len(crosslink_type_filtered_csms["Intra"]) == 803
    assert len(crosslink_type_filtered_csms["Inter"]) == 23


def test6():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_crosslink_type

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    crosslink_type_filtered_crosslinks = filter_crosslink_type(result["crosslinks"])
    assert len(crosslink_type_filtered_crosslinks["Intra"]) == 279
    assert len(crosslink_type_filtered_crosslinks["Inter"]) == 21


def test7():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_protein_distribution

    result = read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    proteins_csms = filter_protein_distribution(result["crosslink-spectrum-matches"])
    proteins_found = list(proteins_csms.keys())  # proteins found
    proteins = [
        "Cas9",
        "sp|MYG_HUMAN|",
        "sp|CAH1_HUMAN|",
        "sp|RETBP_HUMAN|",
        "sp|K1C15_SHEEP|",
    ]
    for p in proteins:
        assert p in proteins_found
    cas9 = len(proteins_csms["Cas9"])  # number of CSMs for protein Cas9
    assert cas9 == 728


def test8():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_peptide_pair_distribution

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    peptide_pairs = filter_peptide_pair_distribution(
        result["crosslink-spectrum-matches"]
    )
    peptide_pairs_found = list(peptide_pairs.keys())[:5]  # first 5 found peptide pairs
    peptide_pairs_should = [
        "GQKNSR-GQKNSR",
        "GQKNSR-GSQKDR",
        "SDKNR-SDKNR",
        "DKQSGK-DKQSGK",
        "DKQSGK-HSIKK",
    ]
    for p in peptide_pairs_should:
        assert p in peptide_pairs_found
    MTNFDKNLPNEK_SKLVSDFR = len(
        peptide_pairs["MTNFDKNLPNEK-SKLVSDFR"]
    )  # number of CSMs for peptide pair MTNFDKNLPNEK-SKLVSDFR
    assert MTNFDKNLPNEK_SKLVSDFR == 21


def test9():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_peptide_pair_distribution
    from pyXLMS.transform import aggregate

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    peptide_pairs = filter_peptide_pair_distribution(
        result["crosslink-spectrum-matches"]
    )
    assert len(peptide_pairs) == len(aggregate(result["crosslink-spectrum-matches"]))
