#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS import data

    err = {"peptide": "PEPTIDE"}
    with pytest.raises(
        TypeError, match="Parameter csm is not a valid crosslink-spectrum-match!"
    ):
        _r = data.create_crosslink_from_csm(err)


def test2():
    from pyXLMS import data

    err = {"data_type": "peptide-spectrum-match", "peptide": "PEPTIDE", "scan": 1}
    with pytest.raises(
        TypeError, match="Parameter csm is not a valid crosslink-spectrum-match!"
    ):
        _r = data.create_crosslink_from_csm(err)


def test3():
    from pyXLMS import data

    csm = data.create_csm_min("PEPKTIDE", 4, "PEPTIKDE", 6, "RUN_1", 1)
    crosslink = data.create_crosslink_from_csm(csm)
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "partial"
    assert crosslink["alpha_peptide"] == "PEPKTIDE"
    assert crosslink["alpha_peptide_crosslink_position"] == 4
    assert crosslink["alpha_proteins"] is None
    assert crosslink["alpha_proteins_crosslink_positions"] is None
    assert crosslink["alpha_decoy"] is None
    assert crosslink["beta_peptide"] == "PEPTIKDE"
    assert crosslink["beta_peptide_crosslink_position"] == 6
    assert crosslink["beta_proteins"] is None
    assert crosslink["beta_proteins_crosslink_positions"] is None
    assert crosslink["beta_decoy"] is None
    assert crosslink["crosslink_type"] == "inter"
    assert crosslink["score"] is None
    assert crosslink["additional_information"] is None


def test4():
    from pyXLMS import data

    csm = data.create_csm_min("PEPKTIDE", 4, "PEPTIKDE", 6, "RUN_1", 1, score=170.3)
    crosslink = data.create_crosslink_from_csm(csm)
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "partial"
    assert crosslink["alpha_peptide"] == "PEPKTIDE"
    assert crosslink["alpha_peptide_crosslink_position"] == 4
    assert crosslink["alpha_proteins"] is None
    assert crosslink["alpha_proteins_crosslink_positions"] is None
    assert crosslink["alpha_decoy"] is None
    assert crosslink["beta_peptide"] == "PEPTIKDE"
    assert crosslink["beta_peptide_crosslink_position"] == 6
    assert crosslink["beta_proteins"] is None
    assert crosslink["beta_proteins_crosslink_positions"] is None
    assert crosslink["beta_decoy"] is None
    assert crosslink["crosslink_type"] == "inter"
    assert crosslink["score"] == pytest.approx(170.3)
    assert crosslink["additional_information"] is None
