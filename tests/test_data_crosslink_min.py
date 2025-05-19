#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS import data

    crosslink = data.create_crosslink_min(
        "PEPTIDE",
        1,
        "EDITPEP",
        3,
    )
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "partial"
    assert crosslink["alpha_peptide"] == "EDITPEP"
    assert crosslink["alpha_peptide_crosslink_position"] == 3
    assert crosslink["alpha_proteins"] is None
    assert crosslink["alpha_proteins_crosslink_positions"] is None
    assert crosslink["alpha_decoy"] is None
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 1
    assert crosslink["beta_proteins"] is None
    assert crosslink["beta_proteins_crosslink_positions"] is None
    assert crosslink["beta_decoy"] is None
    assert crosslink["crosslink_type"] == "inter"
    assert crosslink["score"] is None


def test2():
    from pyXLMS import data

    crosslink = data.create_crosslink_min("PEPTIDE", 1, "EDITPEP", 3, score=170.3)
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "partial"
    assert crosslink["alpha_peptide"] == "EDITPEP"
    assert crosslink["alpha_peptide_crosslink_position"] == 3
    assert crosslink["alpha_proteins"] is None
    assert crosslink["alpha_proteins_crosslink_positions"] is None
    assert crosslink["alpha_decoy"] is None
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 1
    assert crosslink["beta_proteins"] is None
    assert crosslink["beta_proteins_crosslink_positions"] is None
    assert crosslink["beta_decoy"] is None
    assert crosslink["crosslink_type"] == "inter"
    assert crosslink["score"] == pytest.approx(170.3)


def test3():
    from pyXLMS import data

    crosslink = data.create_crosslink_min(
        "PEPTIDE",
        1,
        "EDITPEP",
        3,
        proteins_a=["Protein"],
        xl_position_proteins_a=[15],
        decoy_a=True,
        proteins_b=["Protein"],
        xl_position_proteins_b=[37],
        decoy_b=False,
        score=170.3,
        additional_information={"charge": 3},
    )
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "full"
    assert crosslink["alpha_peptide"] == "EDITPEP"
    assert crosslink["alpha_peptide_crosslink_position"] == 3
    assert crosslink["alpha_proteins"] == ["Protein"]
    assert crosslink["alpha_proteins_crosslink_positions"] == [37]
    assert not crosslink["alpha_decoy"]
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 1
    assert crosslink["beta_proteins"] == ["Protein"]
    assert crosslink["beta_proteins_crosslink_positions"] == [15]
    assert crosslink["beta_decoy"]
    assert crosslink["crosslink_type"] == "intra"
    assert crosslink["score"] == pytest.approx(170.3)
    assert crosslink["additional_information"]["charge"] == 3
