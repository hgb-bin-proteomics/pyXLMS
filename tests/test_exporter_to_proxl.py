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


def test2():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl
    from lxml import etree

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    _xml = to_proxl(
        pr["crosslink-spectrum-matches"],
        fasta_filename="data/_fasta/Cas9_plus10.fasta",
        search_engine="MS Annika",
        search_engine_version="3.0.1",
        score="higher_better",
        crosslinker="DSS",
        filename="online.xml",
        schema_validation="online",
    )
    xmlschema_doc = etree.parse("data/_test/exporter/proxl/proxl-xml.xsd")
    xmlschema = etree.XMLSchema(xmlschema_doc)
    doc = etree.parse("online.xml")
    assert xmlschema.validate(doc)


def test3():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl
    from lxml import etree

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    _xml = to_proxl(
        pr["crosslink-spectrum-matches"],
        fasta_filename="data/_fasta/Cas9_plus10.fasta",
        search_engine="MS Annika",
        search_engine_version="3.0.1",
        score="higher_better",
        crosslinker="DSS",
        filename="offline.xml",
        schema_validation="offline",
    )
    xmlschema_doc = etree.parse("data/_test/exporter/proxl/proxl-xml.xsd")
    xmlschema = etree.XMLSchema(xmlschema_doc)
    doc = etree.parse("offline.xml")
    assert xmlschema.validate(doc)


def test4():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl
    from pyXLMS.transform import fasta_title_to_accession
    from lxml import etree

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
        fasta_filename_override="gTuSC-parsimonious-plusRev.fasta",
        fasta_title_to_accession=fasta_title_to_accession,
        schema_validation="offline",
    )
    xmlschema_doc = etree.parse("data/_test/exporter/proxl/proxl-xml.xsd")
    xmlschema = etree.XMLSchema(xmlschema_doc)
    doc = etree.parse("offline.xml")
    assert xmlschema.validate(doc)
    assert 'fasta_filename="gTuSC-parsimonious-plusRev.fasta"' in xml
    assert '"sp|RETBP_HUMAN|"' not in xml


def test5():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    with pytest.raises(
        TypeError,
        match="Parameter 'score' has to be one of 'higher_better' or 'lower_better'!",
    ):
        _xml = to_proxl(
            pr["crosslink-spectrum-matches"],
            fasta_filename="data/_fasta/Cas9_plus10.fasta",
            search_engine="MS Annika",
            search_engine_version="3.0.1",
            score="greater_better",
            crosslinker="DSS",
            schema_validation="offline",
        )
