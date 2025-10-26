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


def test6():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    with pytest.raises(
        TypeError,
        match="Parameter 'schema_validation' has to be one of 'online' or 'offline'!",
    ):
        _xml = to_proxl(
            pr["crosslink-spectrum-matches"],
            fasta_filename="data/_fasta/Cas9_plus10.fasta",
            search_engine="MS Annika",
            search_engine_version="3.0.1",
            score="higher_better",
            crosslinker="DSS",
            schema_validation="local",
        )


def test7():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    with pytest.raises(
        KeyError,
        match=(
            "Cannot infer crosslinker mass because crosslinker is not defined in "
            r"parameter 'modifications'\. Please specify crosslinker mass manually!"
        ),
    ):
        _xml = to_proxl(
            pr["crosslink-spectrum-matches"],
            fasta_filename="data/_fasta/Cas9_plus10.fasta",
            search_engine="MS Annika",
            search_engine_version="3.0.1",
            score="higher_better",
            crosslinker="SDA",
            schema_validation="offline",
        )


def test8():
    from pyXLMS.exporter import to_proxl

    with pytest.raises(
        ValueError,
        match="Provided crosslink-spectrum-matches contain no elements!",
    ):
        _xml = to_proxl(
            [],
            fasta_filename="data/_fasta/Cas9_plus10.fasta",
            search_engine="MS Annika",
            search_engine_version="3.0.1",
            score="higher_better",
            crosslinker="DSS",
            schema_validation="offline",
        )


def test9():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    with pytest.raises(
        TypeError,
        match="Unsupported data type for input csms! Parameter csms has to be a list of crosslink-spectrum-matches!",
    ):
        _xml = to_proxl(
            pr["crosslinks"],
            fasta_filename="data/_fasta/Cas9_plus10.fasta",
            search_engine="MS Annika",
            search_engine_version="3.0.1",
            score="higher_better",
            crosslinker="DSS",
            schema_validation="offline",
        )


def test10():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_proxl
    from pyXLMS.data import create_csm_min

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    with pytest.raises(
        RuntimeError,
        match="Can't export to ProXL because not all necessary information is available!",
    ):
        _xml = to_proxl(
            pr["crosslink-spectrum-matches"]
            + [create_csm_min("PEPK", 4, "KPEP", 1, "MSFile", 1)],
            fasta_filename="data/_fasta/Cas9_plus10.fasta",
            search_engine="MS Annika",
            search_engine_version="3.0.1",
            score="higher_better",
            crosslinker="DSS",
            schema_validation="offline",
        )


def test11():
    from pyXLMS import parser
    from pyXLMS import exporter
    from pyXLMS.transform import targets_only

    parser_result = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = parser_result["crosslink-spectrum-matches"]
    csms = targets_only(csms)
    xml = exporter.to_proxl(
        csms,
        fasta_filename="data/_fasta/Cas9_plus10.fasta",
        search_engine="MS Annika",
        search_engine_version="3.0.1",
        score="higher_better",
        crosslinker="DSS",
        filename="CSMS_exported_to_ProXL.xml",
        schema_validation="online",
    )
    assert xml is not None
