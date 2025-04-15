#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS import data

    csm = data.create_csm(
        "PEPTIDE",
        {1: (" Oxidation ", 15.994915)},
        1,
        ["PROTEIN"],
        [1],
        [1],
        50.3,
        False,
        "EDITPEP",
        {2: ("Oxidation", 15.994915)},
        3,
        ["NIETORP", "PROTEIN"],
        [7, 4],
        [5, 2],
        170.3,
        False,
        170.3,
        "RUN_1",
        1,
        3,
        23.4,
        -50.0,
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
    # assert csm["alpha_decoy"] == False
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
    # assert csm["beta_decoy"] == False
    assert not csm["beta_decoy"]
    # csm
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(170.3)
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(23.4)
    assert csm["ion_mobility"] == pytest.approx(-50.0)


def test2():
    from pyXLMS import data

    csm = data.create_csm(
        "PEPTIDE",
        {1: (" Oxidation ", 15.994915)},
        1,
        ["PROTEIN"],
        [1],
        [1],
        50.3,
        False,
        "EDITPEP",
        {2: ("Oxidation", 15.994915)},
        3,
        ["NIETORP"],
        [7],
        [5],
        170.3,
        False,
        170.3,
        "RUN_1",
        1,
        3,
        23.4,
        -50.0,
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
    assert len(csm["alpha_proteins"]) == 1
    assert csm["alpha_proteins"][0] == "NIETORP"
    assert len(csm["alpha_proteins_crosslink_positions"]) == 1
    assert csm["alpha_proteins_crosslink_positions"][0] == 7
    assert len(csm["alpha_proteins_peptide_positions"]) == 1
    assert csm["alpha_proteins_peptide_positions"][0] == 5
    assert csm["alpha_score"] == pytest.approx(170.3)
    # assert csm["alpha_decoy"] == False
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
    # assert csm["beta_decoy"] == False
    assert not csm["beta_decoy"]
    # csm
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(170.3)
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(23.4)
    assert csm["ion_mobility"] == pytest.approx(-50.0)


def test3():
    from pyXLMS import data

    csm = data.create_csm(
        "EDITPEP",
        {2: ("Oxidation", 15.994915)},
        3,
        ["NIETORP", "PROTEIN"],
        [7, 4],
        [5, 2],
        170.3,
        False,
        "PEPTIDE",
        {1: (" Oxidation ", 15.994915)},
        1,
        ["PROTEIN"],
        [1],
        [1],
        50.3,
        False,
        170.3,
        "RUN_1",
        1,
        3,
        23.4,
        -50.0,
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
    # assert csm["alpha_decoy"] == False
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
    # assert csm["beta_decoy"] == False
    assert not csm["beta_decoy"]
    # csm
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(170.3)
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(23.4)
    assert csm["ion_mobility"] == pytest.approx(-50.0)


def test4():
    from pyXLMS import data

    csm = data.create_csm(
        "PEPTIDE  ",
        {1: ("    Oxidation ", 15.994915)},
        3,
        ["   PROTEIN"],
        [5],
        [3],
        50.3,
        False,
        "   PEPTIDE",
        {2: ("Oxidation", 15.994915)},
        1,
        ["PROTEIN  "],
        [1],
        [1],
        170.3,
        False,
        170.3,
        "RUN_1   ",
        1,
        3,
        23.4,
        -50.0,
    )
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "full"
    # alpha
    assert csm["alpha_peptide"] == "PEPTIDE"
    assert len(csm["alpha_modifications"]) == 1
    assert 2 in csm["alpha_modifications"]
    assert csm["alpha_modifications"][2][0] == "Oxidation"
    assert csm["alpha_modifications"][2][1] == pytest.approx(15.994915)
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert len(csm["alpha_proteins"]) == 1
    assert csm["alpha_proteins"][0] == "PROTEIN"
    assert len(csm["alpha_proteins_crosslink_positions"]) == 1
    assert csm["alpha_proteins_crosslink_positions"][0] == 1
    assert len(csm["alpha_proteins_peptide_positions"]) == 1
    assert csm["alpha_proteins_peptide_positions"][0] == 1
    assert csm["alpha_score"] == pytest.approx(170.3)
    # assert csm["alpha_decoy"] == False
    assert not csm["alpha_decoy"]
    # beta
    assert csm["beta_peptide"] == "PEPTIDE"
    assert len(csm["beta_modifications"]) == 1
    assert 1 in csm["beta_modifications"]
    assert csm["beta_modifications"][1][0] == "Oxidation"
    assert csm["beta_modifications"][1][1] == pytest.approx(15.994915)
    assert csm["beta_peptide_crosslink_position"] == 3
    assert len(csm["beta_proteins"]) == 1
    assert csm["beta_proteins"][0] == "PROTEIN"
    assert len(csm["beta_proteins_crosslink_positions"]) == 1
    assert csm["beta_proteins_crosslink_positions"][0] == 5
    assert len(csm["beta_proteins_peptide_positions"]) == 1
    assert csm["beta_proteins_peptide_positions"][0] == 3
    assert csm["beta_score"] == pytest.approx(50.3)
    # assert csm["beta_decoy"] == False
    assert not csm["beta_decoy"]
    # csm
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(170.3)
    assert csm["spectrum_file"] == "RUN_1"
    assert csm["scan_nr"] == 1
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(23.4)
    assert csm["ion_mobility"] == pytest.approx(-50.0)


def test5():
    from pyXLMS import data

    csm = data.create_csm(
        "PEPTIDE",
        None,
        3,
        None,
        None,
        None,
        None,
        None,
        "PEPTIDE",
        None,
        1,
        None,
        None,
        None,
        None,
        None,
        None,
        "RUN_1",
        1,
        None,
        None,
        None,
    )
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    # alpha
    assert csm["alpha_peptide"] == "PEPTIDE"
    assert csm["alpha_modifications"] is None
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] is None
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    # beta
    assert csm["beta_peptide"] == "PEPTIDE"
    assert csm["beta_modifications"] is None
    assert csm["beta_peptide_crosslink_position"] == 3
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


def test6():
    from pyXLMS import data

    with pytest.raises(TypeError, match=f"modifications_a must be {dict}!"):
        _csm = data.create_csm(
            "PEPTIDE",
            (" Oxidation ", 15.994915),
            1,
            ["PROTEIN"],
            [1],
            [1],
            50.3,
            False,
            "EDITPEP",
            {2: ("Oxidation", 15.994915)},
            3,
            ["NIETORP", "PROTEIN"],
            [7, 4],
            [5, 2],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )


def test7():
    from pyXLMS import data

    with pytest.raises(
        TypeError, match=f"Dict values of modifications_b must be {tuple}!"
    ):
        _csm = data.create_csm(
            "PEPTIDE",
            {1: (" Oxidation ", 15.994915)},
            1,
            ["PROTEIN"],
            [1],
            [1],
            50.3,
            False,
            "EDITPEP",
            {2: {"Oxidation": 15.994915}},
            3,
            ["NIETORP", "PROTEIN"],
            [7, 4],
            [5, 2],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )


def test8():
    from pyXLMS import data

    with pytest.raises(TypeError, match=f"xl_position_peptide_b must be {int}!"):
        _csm = data.create_csm(
            "PEPTIDE",
            {1: (" Oxidation ", 15.994915)},
            1,
            ["PROTEIN"],
            [1],
            [1],
            50.3,
            False,
            "EDITPEP",
            {2: ("Oxidation", 15.994915)},
            "3",
            ["NIETORP", "PROTEIN"],
            [7, 4],
            [5, 2],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )


def test9():
    from pyXLMS import data

    with pytest.raises(
        TypeError, match=f"List values of xl_position_proteins_b must be {int}!"
    ):
        _csm = data.create_csm(
            "PEPTIDE",
            {1: (" Oxidation ", 15.994915)},
            1,
            ["PROTEIN"],
            [1],
            [1],
            50.3,
            False,
            "EDITPEP",
            {2: ("Oxidation", 15.994915)},
            3,
            ["NIETORP", "PROTEIN"],
            [7, "2"],
            [5, 2],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )


def test10():
    from pyXLMS import data

    with pytest.raises(
        ValueError,
        match="Crosslink position has to be given for every protein! Length of proteins_a and xl_position_proteins_a has to match!",
    ):
        _csm = data.create_csm(
            "PEPTIDE",
            {1: (" Oxidation ", 15.994915)},
            1,
            ["PROTEIN"],
            [1, 2],
            [1],
            50.3,
            False,
            "EDITPEP",
            {2: ("Oxidation", 15.994915)},
            3,
            ["NIETORP", "PROTEIN"],
            [7, 4],
            [5, 2],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )


def test11():
    from pyXLMS import data

    with pytest.raises(
        ValueError,
        match="Crosslink position has to be given for every protein! Length of proteins_b and xl_position_proteins_b has to match!",
    ):
        _csm = data.create_csm(
            "PEPTIDE",
            {1: (" Oxidation ", 15.994915)},
            1,
            ["PROTEIN"],
            [1],
            [1],
            50.3,
            False,
            "EDITPEP",
            {2: ("Oxidation", 15.994915)},
            3,
            ["NIETORP", "PROTEIN"],
            [7],
            [5, 2],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )


def test12():
    from pyXLMS import data

    with pytest.raises(
        ValueError,
        match="Peptide position has to be given for every protein! Length of proteins_a and pep_position_proteins_a has to match!",
    ):
        _csm = data.create_csm(
            "PEPTIDE",
            {1: (" Oxidation ", 15.994915)},
            1,
            ["PROTEIN"],
            [1],
            [],
            50.3,
            False,
            "EDITPEP",
            {2: ("Oxidation", 15.994915)},
            3,
            ["NIETORP", "PROTEIN"],
            [7, 4],
            [5, 2],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )


def test13():
    from pyXLMS import data

    with pytest.raises(
        ValueError,
        match="Peptide position has to be given for every protein! Length of proteins_b and pep_position_proteins_b has to match!",
    ):
        _csm = data.create_csm(
            "PEPTIDE",
            {1: (" Oxidation ", 15.994915)},
            1,
            ["PROTEIN"],
            [1],
            [1],
            50.3,
            False,
            "EDITPEP",
            {2: ("Oxidation", 15.994915)},
            3,
            ["NIETORP", "PROTEIN"],
            [7, 4],
            [5],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )


def test14():
    from pyXLMS import data

    with pytest.raises(
        ValueError,
        match="0-based value found! All positions must use 1-based indexing!",
    ):
        _csm = data.create_csm(
            "PEPTIDE",
            {1: (" Oxidation ", 15.994915)},
            1,
            ["PROTEIN"],
            [1],
            [1],
            50.3,
            False,
            "EDITPEP",
            {2: ("Oxidation", 15.994915)},
            3,
            ["NIETORP", "PROTEIN"],
            [7, 4],
            [5, 4],
            170.3,
            False,
            170.3,
            "RUN_1",
            1,
            3,
            23.4,
            -50.0,
        )
