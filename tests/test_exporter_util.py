#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from pyXLMS.exporter.util import __get_filename

    assert __get_filename("filename", "xlsx") == "filename.xlsx"
    assert __get_filename("filename.xlsx", "xlsx") == "filename.xlsx"
    assert __get_filename("filename.csv", "xlsx") == "filename.csv.xlsx"
    assert __get_filename("filename.csv.xlsx", "xlsx") == "filename.csv.xlsx"
