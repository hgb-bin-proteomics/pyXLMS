#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import filter_proteins
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    cas9 = filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]
    export = to_alphalink2(
        cas9, fasta="data/_fasta/Cas9_plus10.fasta", filename_prefix="Cas9"
    )
    assert export["Exported files"][0] == "Cas9_AlphaLink.txt"
    assert export["Exported files"][1] == "Cas9_AlphaLink.fasta"
