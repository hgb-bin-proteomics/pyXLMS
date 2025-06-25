#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.constants import AMINO_ACIDS

    assert "A" in AMINO_ACIDS
    assert "B" not in AMINO_ACIDS


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
