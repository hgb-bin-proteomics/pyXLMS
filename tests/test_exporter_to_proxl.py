#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    xml = to_proxl(
        pr["crosslink-spectrum-matches"],
        fasta_filename="data/_fasta/Cas9_plus10.fasta",
        search_engine="MS Annika",
        search_engine_version="3.0.1",
        score="higher_better",
        crosslinker="DSS",
        filename="DSS_Cas9_ProXL.xml",
    )
    assert xml is not None
