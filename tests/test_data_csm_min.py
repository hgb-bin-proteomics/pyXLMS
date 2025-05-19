#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS import data

    csm = data.create_csm_min(
        "PEPTIDE",
        1,
        "EDITPEP",
        3,
        "RUN_1",
        1,
    )
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    # alpha
    assert csm["alpha_peptide"] == "EDITPEP"
    assert csm["alpha_modifications"] is None
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] is None
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    # beta
    assert csm["beta_peptide"] == "PEPTIDE"
    assert csm["beta_modifications"] is None
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] is None
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    # csm
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] is None
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] is None
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test2():
    from pyXLMS import data

    csm = data.create_csm_min("PEPTIDE", 1, "EDITPEP", 3, "RUN_1", 1, score=170.3)
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    # alpha
    assert csm["alpha_peptide"] == "EDITPEP"
    assert csm["alpha_modifications"] is None
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] is None
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    # beta
    assert csm["beta_peptide"] == "PEPTIDE"
    assert csm["beta_modifications"] is None
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] is None
    assert csm["beta_proteins_crosslink_positions"] is None
    assert csm["beta_proteins_peptide_positions"] is None
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    # csm
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(170.3)
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] is None
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test3():
    from pyXLMS import data

    csm = data.create_csm_min(
        "PEPTIDE",
        1,
        "EDITPEP",
        3,
        "RUN_1",
        1,
        modifications_a={1: (" Oxidation ", 15.994915)},
        proteins_a=["PROTEIN"],
        xl_position_proteins_a=[1],
        pep_position_proteins_a=[1],
        score_a=50.3,
        decoy_a=False,
        modifications_b={2: ("Oxidation", 15.994915)},
        proteins_b=["NIETORP", "PROTEIN"],
        xl_position_proteins_b=[7, 4],
        pep_position_proteins_b=[5, 2],
        score_b=170.3,
        decoy_b=False,
        score=170.3,
        charge=3,
        rt=23.4,
        im_cv=-50.0,
        additional_information={"m/z": 1337.1},
    )
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "full"
    # alpha
    assert csm["alpha_peptide"] == "EDITPEP"
    assert len(csm["alpha_modifications"]) == 1
    assert 2 in csm["alpha_modifications"]
    assert csm["alpha_modifications"][2][0] == "Oxidation"
    assert csm["alpha_modifications"][2][1] == pytest.approx(15.994915)
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert len(csm["alpha_proteins"]) == 2
    assert csm["alpha_proteins"][0] == "NIETORP"
    assert len(csm["alpha_proteins_crosslink_positions"]) == 2
    assert csm["alpha_proteins_crosslink_positions"][0] == 7
    assert len(csm["alpha_proteins_peptide_positions"]) == 2
    assert csm["alpha_proteins_peptide_positions"][0] == 5
    assert csm["alpha_score"] == pytest.approx(170.3)
    assert not csm["alpha_decoy"]
    # beta
    assert csm["beta_peptide"] == "PEPTIDE"
    assert len(csm["beta_modifications"]) == 1
    assert 1 in csm["beta_modifications"]
    assert csm["beta_modifications"][1][0] == "Oxidation"
    assert csm["beta_modifications"][1][1] == pytest.approx(15.994915)
    assert csm["beta_peptide_crosslink_position"] == 1
    assert len(csm["beta_proteins"]) == 1
    assert csm["beta_proteins"][0] == "PROTEIN"
    assert len(csm["beta_proteins_crosslink_positions"]) == 1
    assert csm["beta_proteins_crosslink_positions"][0] == 1
    assert len(csm["beta_proteins_peptide_positions"]) == 1
    assert csm["beta_proteins_peptide_positions"][0] == 1
    assert csm["beta_score"] == pytest.approx(50.3)
    assert not csm["beta_decoy"]
    # csm
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(170.3)
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(23.4)
    assert csm["ion_mobility"] == pytest.approx(-50.0)
    assert csm["additional_information"]["m/z"] == pytest.approx(1337.1)
