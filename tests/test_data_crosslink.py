#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.data import check_input

    x = 1
    assert check_input(x, "x", int)


def test2():
    from pyXLMS import data

    x = 1
    assert data.check_input(x, "x", int)


def test3():
    from pyXLMS import data

    x = [1, 2, 3]
    assert data.check_input(x, "x", list)


def test4():
    from pyXLMS import data

    x = [1, 2, 3]
    assert data.check_input(x, "x", list, int)


def test5():
    from pyXLMS import data

    x = 1
    with pytest.raises(TypeError, match=f"x must be {str}!"):
        _i = data.check_input(x, "x", str)


def test6():
    from pyXLMS import data

    x = [1, 2, 3]
    with pytest.raises(TypeError, match=f"List values of x must be {str}!"):
        _i = data.check_input(x, "x", list, str)


def test7():
    from pyXLMS import data

    x = {"a": 1, "b": 2, "c": 3}
    with pytest.raises(TypeError, match=f"Dict values of x must be {str}!"):
        _i = data.check_input(x, "x", dict, str)


def test8():
    from pyXLMS import data

    x = 1
    assert data.check_input_multi(x, "x", [int, list])


def test9():
    from pyXLMS import data

    x = [1, 2, 3]
    assert data.check_input_multi(x, "x", [int, list])


def test10():
    from pyXLMS import data

    x = [1, 2, 3]
    assert data.check_input_multi(x, "x", [int, list], int)


def test11():
    from pyXLMS import data

    x = 1
    with pytest.raises(
        TypeError, match=f"x must be one of {','.join([str(str), str(list)])}!"
    ):
        _i = data.check_input_multi(x, "x", [str, list])


def test12():
    from pyXLMS import data

    x = [1, 2, 3]
    with pytest.raises(TypeError, match=f"List values of x must be {str}!"):
        _i = data.check_input_multi(x, "x", [str, list], str)


def test13():
    from pyXLMS import data

    x = {"a": 1, "b": 2, "c": 3}
    with pytest.raises(TypeError, match=f"Dict values of x must be {str}!"):
        _i = data.check_input_multi(x, "x", [str, dict], str)


def test14():
    from pyXLMS import data

    x = 1
    assert data.check_indexing(x)


def test15():
    from pyXLMS import data

    x = [1, 2]
    assert data.check_indexing(x)


def test16():
    from pyXLMS import data

    x = 0
    with pytest.raises(
        ValueError,
        match="0-based value found! All positions must use 1-based indexing!",
    ):
        _i = data.check_indexing(x)


def test17():
    from pyXLMS import data

    x = [1, 2, 0]
    with pytest.raises(
        ValueError,
        match="0-based value found! All positions must use 1-based indexing!",
    ):
        _i = data.check_indexing(x)


def test18():
    from pyXLMS import data

    crosslink = data.create_crosslink(
        "PEPTIDE",
        1,
        ["PROTEIN"],
        [1],
        False,
        "EDITPEP",
        3,
        ["NIETORP", "PROTEIN"],
        [5, 2],
        False,
        170.3,
    )
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "full"
    assert crosslink["alpha_peptide"] == "EDITPEP"
    assert crosslink["alpha_peptide_crosslink_position"] == 3
    assert len(crosslink["alpha_proteins"]) == 2
    assert crosslink["alpha_proteins"][0] == "NIETORP"
    assert len(crosslink["alpha_proteins_crosslink_positions"]) == 2
    assert crosslink["alpha_proteins_crosslink_positions"][0] == 5
    # assert crosslink["alpha_decoy"] == False
    assert not crosslink["alpha_decoy"]
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 1
    assert len(crosslink["beta_proteins"]) == 1
    assert crosslink["beta_proteins"][0] == "PROTEIN"
    assert len(crosslink["beta_proteins_crosslink_positions"]) == 1
    assert crosslink["beta_proteins_crosslink_positions"][0] == 1
    # assert crosslink["beta_decoy"] == False
    assert not crosslink["beta_decoy"]
    assert crosslink["crosslink_type"] == "intra"
    assert crosslink["score"] == pytest.approx(170.3)


def test19():
    from pyXLMS import data

    crosslink = data.create_crosslink(
        "EDITPEP",
        3,
        ["NIETORP", "PROTEIN"],
        [5, 2],
        False,
        "PEPTIDE",
        1,
        ["PROTEIN"],
        [1],
        False,
        170.3,
    )
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "full"
    assert crosslink["alpha_peptide"] == "EDITPEP"
    assert crosslink["alpha_peptide_crosslink_position"] == 3
    assert len(crosslink["alpha_proteins"]) == 2
    assert crosslink["alpha_proteins"][0] == "NIETORP"
    assert len(crosslink["alpha_proteins_crosslink_positions"]) == 2
    assert crosslink["alpha_proteins_crosslink_positions"][0] == 5
    # assert crosslink["alpha_decoy"] == False
    assert not crosslink["alpha_decoy"]
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 1
    assert len(crosslink["beta_proteins"]) == 1
    assert crosslink["beta_proteins"][0] == "PROTEIN"
    assert len(crosslink["beta_proteins_crosslink_positions"]) == 1
    assert crosslink["beta_proteins_crosslink_positions"][0] == 1
    # assert crosslink["beta_decoy"] == False
    assert not crosslink["beta_decoy"]
    assert crosslink["crosslink_type"] == "intra"
    assert crosslink["score"] == pytest.approx(170.3)


def test20():
    from pyXLMS import data

    crosslink = data.create_crosslink(
        "EDITPEP    ",
        3,
        ["    NIETORP", "PROTEIN  "],
        [5, 2],
        False,
        "     PEPTIDE",
        1,
        [" PROTEIN "],
        [1],
        False,
        170.3,
    )
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "full"
    assert crosslink["alpha_peptide"] == "EDITPEP"
    assert crosslink["alpha_peptide_crosslink_position"] == 3
    assert len(crosslink["alpha_proteins"]) == 2
    assert crosslink["alpha_proteins"][0] == "NIETORP"
    assert crosslink["alpha_proteins"][1] == "PROTEIN"
    assert len(crosslink["alpha_proteins_crosslink_positions"]) == 2
    assert crosslink["alpha_proteins_crosslink_positions"][0] == 5
    # assert crosslink["alpha_decoy"] == False
    assert not crosslink["alpha_decoy"]
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 1
    assert len(crosslink["beta_proteins"]) == 1
    assert crosslink["beta_proteins"][0] == "PROTEIN"
    assert len(crosslink["beta_proteins_crosslink_positions"]) == 1
    assert crosslink["beta_proteins_crosslink_positions"][0] == 1
    # assert crosslink["beta_decoy"] == False
    assert not crosslink["beta_decoy"]
    assert crosslink["crosslink_type"] == "intra"
    assert crosslink["score"] == pytest.approx(170.3)


def test21():
    from pyXLMS import data

    crosslink = data.create_crosslink(
        "PEPTIDE",
        3,
        ["PROTEIN"],
        [3],
        True,
        "PEPTIDE",
        1,
        ["PROTEIN"],
        [1],
        False,
        170.3,
    )
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "full"
    assert crosslink["alpha_peptide"] == "PEPTIDE"
    assert crosslink["alpha_peptide_crosslink_position"] == 1
    assert len(crosslink["alpha_proteins"]) == 1
    assert crosslink["alpha_proteins"][0] == "PROTEIN"
    assert len(crosslink["alpha_proteins_crosslink_positions"]) == 1
    assert crosslink["alpha_proteins_crosslink_positions"][0] == 1
    # assert crosslink["alpha_decoy"] == False
    assert not crosslink["alpha_decoy"]
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 3
    assert len(crosslink["beta_proteins"]) == 1
    assert crosslink["beta_proteins"][0] == "PROTEIN"
    assert len(crosslink["beta_proteins_crosslink_positions"]) == 1
    assert crosslink["beta_proteins_crosslink_positions"][0] == 3
    # assert crosslink["beta_decoy"] == True
    assert crosslink["beta_decoy"]
    assert crosslink["crosslink_type"] == "intra"
    assert crosslink["score"] == pytest.approx(170.3)


def test22():
    from pyXLMS import data

    crosslink = data.create_crosslink(
        "PEPTIDE",
        3,
        None,
        None,
        None,
        "PEPTIDE",
        1,
        None,
        None,
        None,
        None,
    )
    assert crosslink["data_type"] == "crosslink"
    assert crosslink["completeness"] == "partial"
    assert crosslink["alpha_peptide"] == "PEPTIDE"
    assert crosslink["alpha_peptide_crosslink_position"] == 1
    assert crosslink["alpha_proteins"] is None
    assert crosslink["alpha_proteins_crosslink_positions"] is None
    assert crosslink["alpha_decoy"] is None
    assert crosslink["beta_peptide"] == "PEPTIDE"
    assert crosslink["beta_peptide_crosslink_position"] == 3
    assert crosslink["beta_proteins"] is None
    assert crosslink["beta_proteins_crosslink_positions"] is None
    assert crosslink["beta_decoy"] is None
    assert crosslink["crosslink_type"] == "inter"
    assert crosslink["score"] is None


def test23():
    from pyXLMS import data

    with pytest.raises(TypeError, match=f"xl_position_peptide_a must be {int}!"):
        _crosslink = data.create_crosslink(
            "PEPTIDE",
            "3",
            ["PROTEIN"],
            [3],
            False,
            "PEPTIDE",
            1,
            ["PROTEIN"],
            [1],
            False,
            170.3,
        )


def test24():
    from pyXLMS import data

    with pytest.raises(
        TypeError, match=f"List values of xl_position_proteins_a must be {int}!"
    ):
        _crosslink = data.create_crosslink(
            "PEPTIDE",
            3,
            ["PROTEIN"],
            ["3"],
            False,
            "PEPTIDE",
            1,
            ["PROTEIN"],
            [1],
            True,
            170.3,
        )


def test25():
    from pyXLMS import data

    with pytest.raises(
        ValueError,
        match="Crosslink position has to be given for every protein! Length of proteins_a and xl_position_proteins_a has to match!",
    ):
        _crosslink = data.create_crosslink(
            "PEPTIDE",
            3,
            ["PROTEIN"],
            [3, 4],
            True,
            "PEPTIDE",
            1,
            ["PROTEIN"],
            [1],
            False,
            170.3,
        )


def test26():
    from pyXLMS import data

    with pytest.raises(
        ValueError,
        match="Crosslink position has to be given for every protein! Length of proteins_b and xl_position_proteins_b has to match!",
    ):
        _crosslink = data.create_crosslink(
            "PEPTIDE",
            3,
            ["PROTEIN"],
            [3],
            True,
            "PEPTIDE",
            1,
            ["PROTEIN"],
            [1, 2],
            True,
            170.3,
        )
