#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


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
        result["crosslink-spectrum-matches"], prefix_decoys=False
    )
    peptide_pairs_found = list(peptide_pairs.keys())[:5]  # first 5 found peptide pairs
    peptide_pairs_should = [
        "GQKNSR:3-GQKNSR:3",
        "GQKNSR:3-GSQKDR:4",
        "SDKNR:3-SDKNR:3",
        "DKQSGK:2-DKQSGK:2",
        "DKQSGK:2-HSIKK:4",
    ]
    for p in peptide_pairs_should:
        assert p in peptide_pairs_found
    MTNFDKNLPNEK_SKLVSDFR = len(
        peptide_pairs["MTNFDKNLPNEK:6-SKLVSDFR:2"]
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
        result["crosslink-spectrum-matches"], prefix_decoys=False
    )
    assert len(peptide_pairs) == len(aggregate(result["crosslink-spectrum-matches"]))


def test10():
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
        "GQKNSR:3-GQKNSR:3",
        "GQKNSR:3-DECOY_GSQKDR:4",
        "SDKNR:3-SDKNR:3",
        "DKQSGK:2-DKQSGK:2",
        "DKQSGK:2-HSIKK:4",
    ]
    for p in peptide_pairs_should:
        assert p in peptide_pairs_found
    MTNFDKNLPNEK_SKLVSDFR = len(
        peptide_pairs["MTNFDKNLPNEK:6-SKLVSDFR:2"]
    )  # number of CSMs for peptide pair MTNFDKNLPNEK-SKLVSDFR
    assert MTNFDKNLPNEK_SKLVSDFR == 21


def test11():
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


def test12():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    residue_pairs = filter_residue_pair_distribution(
        result["crosslink-spectrum-matches"]
    )
    residue_pairs_found = list(residue_pairs.keys())[:5]
    residue_pairs_should = [
        "Cas9:779-Cas9:779",
        "Cas9:779-DECOY_Cas9:696",
        "Cas9:866-Cas9:866",
        "Cas9:677-Cas9:677",
        "Cas9:48-Cas9:677",
    ]
    for r in residue_pairs_should:
        assert r in residue_pairs_found
    assert len(residue_pairs["Cas9:1122-Cas9:884"]) == 22


def test13():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution
    from pyXLMS.transform import aggregate

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    residue_pairs = filter_residue_pair_distribution(
        result["crosslink-spectrum-matches"]
    )
    assert len(residue_pairs) == len(
        aggregate(result["crosslink-spectrum-matches"], by="protein")
    )


def test14():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    residue_pairs = filter_residue_pair_distribution(
        result["crosslink-spectrum-matches"], prefix_decoys=False
    )
    residue_pairs_found = list(residue_pairs.keys())[:5]
    residue_pairs_should = [
        "Cas9:779-Cas9:779",
        "Cas9:696-Cas9:779",
        "Cas9:866-Cas9:866",
        "Cas9:677-Cas9:677",
        "Cas9:48-Cas9:677",
    ]
    for r in residue_pairs_should:
        assert r in residue_pairs_found
    assert len(residue_pairs["Cas9:1122-Cas9:884"]) == 22


def test15():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution
    from pyXLMS.transform import aggregate

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    residue_pairs = filter_residue_pair_distribution(
        result["crosslink-spectrum-matches"], prefix_decoys=False
    )
    assert len(residue_pairs) == len(
        aggregate(result["crosslink-spectrum-matches"], by="protein")
    )


def test16():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = result.csms()
    csms[0] = csms[0].copy_with_update({"alpha_proteins": None})

    with pytest.raises(
        ValueError,
        match="Attribute .* is missing in at least one element but is required!",
    ):
        _ = filter_residue_pair_distribution(csms)


def test17():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = result.csms()
    csms[0] = csms[0].copy_with_update({"beta_proteins": None})

    with pytest.raises(
        ValueError,
        match="Attribute .* is missing in at least one element but is required!",
    ):
        _ = filter_residue_pair_distribution(csms)


def test18():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = result.csms()
    csms[0] = csms[0].copy_with_update({"alpha_proteins_crosslink_positions": None})

    with pytest.raises(
        ValueError,
        match="Attribute .* is missing in at least one element but is required!",
    ):
        _ = filter_residue_pair_distribution(csms)


def test19():
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = result.csms()
    csms[0] = csms[0].copy_with_update({"alpha_proteins_crosslink_positions": None})

    with pytest.raises(
        ValueError,
        match="Attribute .* is missing in at least one element but is required!",
    ):
        _ = filter_residue_pair_distribution(csms)


def test20():
    # this is technically test for transform.util.get_available_keys for @param always_revalidate
    from pyXLMS.parser import read
    from pyXLMS.transform import filter_residue_pair_distribution

    result = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = result.csms()
    csms[0] = csms[0].copy_with_update({"alpha_proteins_crosslink_positions": None})

    with pytest.raises(
        ValueError,
        match="Attribute .* is missing in at least one element but is required!",
    ):
        _ = filter_residue_pair_distribution(csms)
