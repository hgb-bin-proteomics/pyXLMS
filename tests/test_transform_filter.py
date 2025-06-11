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
