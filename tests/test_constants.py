#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.constants import AMINO_ACIDS

    assert "A" in AMINO_ACIDS
    assert "B" in AMINO_ACIDS


def test2():
    from pyXLMS.constants import AMINO_ACIDS_3TO1

    assert AMINO_ACIDS_3TO1["GLY"] == "G"


def test3():
    from pyXLMS.constants import AMINO_ACIDS_1TO3

    assert AMINO_ACIDS_1TO3["G"] == "GLY"


def test4():
    from pyXLMS.constants import CROSSLINKERS

    assert CROSSLINKERS["BS3"] == pytest.approx(138.06808)


def test5():
    from pyXLMS.constants import MODIFICATIONS

    assert MODIFICATIONS["Carbamidomethyl"] == pytest.approx(57.021464)
    assert MODIFICATIONS["BS3"] == pytest.approx(138.06808)


def test6():
    from pyXLMS.constants import XI_MODIFICATION_MAPPING

    assert XI_MODIFICATION_MAPPING["cm"][0] == "Carbamidomethyl"
    assert XI_MODIFICATION_MAPPING["cm"][1] == pytest.approx(57.021464)
    assert XI_MODIFICATION_MAPPING["ox"][0] == "Oxidation"
    assert XI_MODIFICATION_MAPPING["ox"][1] == pytest.approx(15.994915)


def test7():
    from pyXLMS.constants import SCOUT_MODIFICATION_MAPPING

    assert SCOUT_MODIFICATION_MAPPING["+57.021460"][0] == "Carbamidomethyl"
    assert SCOUT_MODIFICATION_MAPPING["+57.021460"][1] == pytest.approx(57.021464)
    assert SCOUT_MODIFICATION_MAPPING["Carbamidomethyl"][0] == "Carbamidomethyl"
    assert SCOUT_MODIFICATION_MAPPING["Carbamidomethyl"][1] == pytest.approx(57.021464)
    assert SCOUT_MODIFICATION_MAPPING["Oxidation of Methionine"][0] == "Oxidation"
    assert SCOUT_MODIFICATION_MAPPING["Oxidation of Methionine"][1] == pytest.approx(
        15.994915
    )


def test8():
    from pyXLMS.constants import MEROX_MODIFICATION_MAPPING

    assert MEROX_MODIFICATION_MAPPING["B"] == {
        "Amino Acid": "C",
        "Modification": ("Carbamidomethyl", 57.021464),
    }
    assert MEROX_MODIFICATION_MAPPING["m"] == {
        "Amino Acid": "M",
        "Modification": ("Oxidation", 15.994915),
    }
    assert MEROX_MODIFICATION_MAPPING == {
        "B": {"Amino Acid": "C", "Modification": ("Carbamidomethyl", 57.021464)},
        "m": {"Amino Acid": "M", "Modification": ("Oxidation", 15.994915)},
    }


def test9():
    from pyXLMS.constants import AMINO_ACIDS_REPLACEMENTS

    assert "X" in AMINO_ACIDS_REPLACEMENTS
