#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.parser import util as p

    assert p.format_sequence("PEP[K]TIDE") == "PEPKTIDE"
    assert p.format_sequence("PEPKdssoTIDE") == "PEPKTIDE"
    assert p.format_sequence("peptide", remove_lower=False) == "PEPTIDE"


def test2():
    from pyXLMS.parser import util as p

    assert p.get_bool_from_value(True)
    assert not p.get_bool_from_value(False)
    assert p.get_bool_from_value("True")
    assert not p.get_bool_from_value("False")
    assert p.get_bool_from_value("true")
    assert not p.get_bool_from_value("false")
    assert p.get_bool_from_value("t")
    assert not p.get_bool_from_value("f")
    assert p.get_bool_from_value(1)
    assert not p.get_bool_from_value(0)


def test3():
    from pyXLMS.parser import util as p

    value = 2
    with pytest.raises(
        ValueError, match=f"Cannot parse bool value from the given input {value}."
    ):
        _b = p.get_bool_from_value(value)


def test4():
    from pyXLMS.parser import util as p

    value = 2.0
    with pytest.raises(
        ValueError, match=f"Cannot parse bool value from the given input {value}."
    ):
        _b = p.get_bool_from_value(value)


def test5():
    from pyXLMS.parser import util as p

    with pytest.raises(
        ValueError, match=r"Cannot parse bool value from the given input \[\]."
    ):
        _b = p.get_bool_from_value([])
