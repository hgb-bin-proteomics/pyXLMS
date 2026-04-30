#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from pyXLMS.parser import read

    pr = read(
        "data/_test/annotate_string_scores/Nucleus_Rep1_Crosslinks.txt.xz",
        engine="MS Annika",
        crosslinker="DSBSO",
        format="txt",
        compression="xz",
        unsafe=True,
        verbose=0,
    )
    assert pr["crosslinks"] is not None
