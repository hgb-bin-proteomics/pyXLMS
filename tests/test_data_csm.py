#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

def test1():
    from pyXLMS import data
    csm = data.create_csm("PEPTIDE", {1: (" Ox ", 16.0)}, 1, ["PROTEIN"], [1], [1], 50.3, False,
                          "EDITPEP", {2: ("Ox", 16.0)}, 3, ["NIETORP", "PROTEIN"], [5, 2], [5, 2], 170.3, False,
                          170.3, "RUN_1", 1, 3, 23.4, -50.0)
    assert csm["data_type"] == "crosslink-spectrum-match"
    # alpha
    assert csm["alpha_peptide"] == "EDITPEP"
    assert len(csm["alpha_modifications"]) == 1
    assert 2 in csm["alpha_modifications"]
    assert csm["alpha_modifications"][2][0] == "Ox"
    assert csm["alpha_modifications"][2][1] >= 15.95 and csm["alpha_modifications"][2][1] <= 16.05
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert len(csm["alpha_proteins"]) == 2
    assert csm["alpha_proteins"][0] == "NIETORP"
    assert len(csm["alpha_proteins_crosslink_positions"]) == 2
    assert csm["alpha_proteins_crosslink_positions"][0] == 5
    assert len(csm["alpha_proteins_peptide_positions"]) == 2
    assert csm["alpha_proteins_peptide_positions"][0] == 5
    assert csm["alpha_score"] >= 170.25 and csm["alpha_score"] <= 170.35
    assert csm["alpha_decoy"] == False
    # beta
    assert csm["beta_peptide"] == "PEPTIDE"
    assert len(csm["beta_modifications"]) == 1
    assert 1 in csm["beta_modifications"]
    assert csm["beta_modifications"][1][0] == "Ox"
    assert csm["beta_modifications"][1][1] >= 15.95 and csm["beta_modifications"][1][1] <= 16.05
    assert csm["beta_peptide_crosslink_position"] == 1
    assert len(csm["beta_proteins"]) == 1
    assert csm["beta_proteins"][0] == "PROTEIN"
    assert len(csm["beta_proteins_crosslink_positions"]) == 1
    assert csm["beta_proteins_crosslink_positions"][0] == 1
    assert len(csm["beta_proteins_peptide_positions"]) == 1
    assert csm["beta_proteins_peptide_positions"][0] == 1
    assert csm["beta_score"] >= 50.25 and csm["beta_score"] <= 50.35
    assert csm["beta_decoy"] == False
    # csm
    assert csm["score"] >= 170.25 and csm["score"] <= 170.35
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] == 3
    assert csm["retention_time"] >= 23.35 and csm["retention_time"] <= 23.45
    assert csm["ion_mobility"] >= -50.05 and csm["ion_mobility"] <= -49.95

def test2():
    from pyXLMS import data
    csm = data.create_csm("EDITPEP", {2: ("Ox", 16.0)}, 3, ["NIETORP", "PROTEIN"], [5, 2], [5, 2], 170.3, False,
                          "PEPTIDE", {1: (" Ox ", 16.0)}, 1, ["PROTEIN"], [1], [1], 50.3, False,
                          170.3, "RUN_1", 1, 3, 23.4, -50.0)
    assert csm["data_type"] == "crosslink-spectrum-match"
    # alpha
    assert csm["alpha_peptide"] == "EDITPEP"
    assert len(csm["alpha_modifications"]) == 1
    assert 2 in csm["alpha_modifications"]
    assert csm["alpha_modifications"][2][0] == "Ox"
    assert csm["alpha_modifications"][2][1] >= 15.95 and csm["alpha_modifications"][2][1] <= 16.05
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert len(csm["alpha_proteins"]) == 2
    assert csm["alpha_proteins"][0] == "NIETORP"
    assert len(csm["alpha_proteins_crosslink_positions"]) == 2
    assert csm["alpha_proteins_crosslink_positions"][0] == 5
    assert len(csm["alpha_proteins_peptide_positions"]) == 2
    assert csm["alpha_proteins_peptide_positions"][0] == 5
    assert csm["alpha_score"] >= 170.25 and csm["alpha_score"] <= 170.35
    assert csm["alpha_decoy"] == False
    # beta
    assert csm["beta_peptide"] == "PEPTIDE"
    assert len(csm["beta_modifications"]) == 1
    assert 1 in csm["beta_modifications"]
    assert csm["beta_modifications"][1][0] == "Ox"
    assert csm["beta_modifications"][1][1] >= 15.95 and csm["beta_modifications"][1][1] <= 16.05
    assert csm["beta_peptide_crosslink_position"] == 1
    assert len(csm["beta_proteins"]) == 1
    assert csm["beta_proteins"][0] == "PROTEIN"
    assert len(csm["beta_proteins_crosslink_positions"]) == 1
    assert csm["beta_proteins_crosslink_positions"][0] == 1
    assert len(csm["beta_proteins_peptide_positions"]) == 1
    assert csm["beta_proteins_peptide_positions"][0] == 1
    assert csm["beta_score"] >= 50.25 and csm["beta_score"] <= 50.35
    assert csm["beta_decoy"] == False
    # csm
    assert csm["score"] >= 170.25 and csm["score"] <= 170.35
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] == 3
    assert csm["retention_time"] >= 23.35 and csm["retention_time"] <= 23.45
    assert csm["ion_mobility"] >= -50.05 and csm["ion_mobility"] <= -49.95

def test3():
    from pyXLMS import data
    csm = data.create_csm("PEPTIDE  ", {1: ("    Ox ", 16.0)}, 3, ["   PROTEIN"], [3], [3], 50.3, False,
                          "   PEPTIDE", {2: ("Ox", 16.0)}, 1, ["PROTEIN  "], [1], [1], 170.3, False,
                          170.3, "RUN_1   ", 1, 3, 23.4, -50.0)
    assert csm["data_type"] == "crosslink-spectrum-match"
    # alpha
    assert csm["alpha_peptide"] == "PEPTIDE"
    assert len(csm["alpha_modifications"]) == 1
    assert 2 in csm["alpha_modifications"]
    assert csm["alpha_modifications"][2][0] == "Ox"
    assert csm["alpha_modifications"][2][1] >= 15.95 and csm["alpha_modifications"][2][1] <= 16.05
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert len(csm["alpha_proteins"]) == 1
    assert csm["alpha_proteins"][0] == "PROTEIN"
    assert len(csm["alpha_proteins_crosslink_positions"]) == 1
    assert csm["alpha_proteins_crosslink_positions"][0] == 1
    assert len(csm["alpha_proteins_peptide_positions"]) == 1
    assert csm["alpha_proteins_peptide_positions"][0] == 1
    assert csm["alpha_score"] >= 170.25 and csm["alpha_score"] <= 170.35
    assert csm["alpha_decoy"] == False
    # beta
    assert csm["beta_peptide"] == "PEPTIDE"
    assert len(csm["beta_modifications"]) == 1
    assert 1 in csm["beta_modifications"]
    assert csm["beta_modifications"][1][0] == "Ox"
    assert csm["beta_modifications"][1][1] >= 15.95 and csm["beta_modifications"][1][1] <= 16.05
    assert csm["beta_peptide_crosslink_position"] == 3
    assert len(csm["beta_proteins"]) == 1
    assert csm["beta_proteins"][0] == "PROTEIN"
    assert len(csm["beta_proteins_crosslink_positions"]) == 1
    assert csm["beta_proteins_crosslink_positions"][0] == 3
    assert len(csm["beta_proteins_peptide_positions"]) == 1
    assert csm["beta_proteins_peptide_positions"][0] == 3
    assert csm["beta_score"] >= 50.25 and csm["beta_score"] <= 50.35
    assert csm["beta_decoy"] == False
    # csm
    assert csm["score"] >= 170.25 and csm["score"] <= 170.35
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] == 3
    assert csm["retention_time"] >= 23.35 and csm["retention_time"] <= 23.45
    assert csm["ion_mobility"] >= -50.05 and csm["ion_mobility"] <= -49.95
