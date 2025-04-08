#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from pyXLMS import transform

    modifications = {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)}
    modifications_str = "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    assert transform.modifications_to_str(modifications) == modifications_str


def test2():
    from pyXLMS import transform

    modifications = {1: ("Oxidation", 15.994915)}
    modifications_str = "(1:[Oxidation|15.994915])"
    assert transform.modifications_to_str(modifications) == modifications_str


def test3():
    from pyXLMS import transform

    modifications = dict()
    modifications_str = ""
    assert transform.modifications_to_str(modifications) == modifications_str


def test4():
    from pyXLMS import transform

    modifications = None
    assert transform.modifications_to_str(modifications) is None
