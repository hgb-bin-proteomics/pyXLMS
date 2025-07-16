#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


FASTA_SMALL = "data/_fasta/Cas9_plus10.fasta"
FASTA_LARGE = "data/_fasta/uniprotkb_proteome_UP000005640_AND_revi_2025_05_23.fasta"


def build_protein_db(fasta: str):
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    def fasta_title_to_accession(title: str):
        if "|" in title:
            return title.split("|")[1].strip()
        return title.strip()

    protein_db = dict()
    with open(fasta, "r", encoding="utf-8") as f:
        for item in SimpleFastaParser(f):
            protein_db[fasta_title_to_accession(item[0])] = item[1]
    return protein_db


def test1():
    from pyXLMS.transform.reannotate_positions import __get_proteins_and_positions

    peptide = "ADANLDK"
    with pytest.raises(RuntimeError, match=f"No match found for peptide {peptide}!"):
        _r = __get_proteins_and_positions(peptide, {})


def test2():
    from pyXLMS.transform.reannotate_positions import __get_proteins_and_positions

    peptide = "ADANLDK"
    with pytest.raises(RuntimeError, match=f"No match found for peptide {peptide}!"):
        _r = __get_proteins_and_positions(peptide, {"PROT1": "PEPKTIDE"})


def test3():
    from pyXLMS.transform.reannotate_positions import __get_proteins_and_positions

    protein_db = build_protein_db(FASTA_SMALL)

    peptide = "ADANLDK"
    prot_pos = __get_proteins_and_positions(peptide, protein_db)
    assert prot_pos[0] == ["Cas9"]
    assert prot_pos[1] == [1286]


def test4():
    from pyXLMS.transform.reannotate_positions import __get_proteins_and_positions

    protein_db = build_protein_db(FASTA_SMALL)

    peptide = "GNTDRHSIK"
    prot_pos = __get_proteins_and_positions(peptide, protein_db)
    assert prot_pos[0] == ["Cas9"]
    assert prot_pos[1] == [39]


def test5():
    from pyXLMS.transform.reannotate_positions import __get_proteins_and_positions

    protein_db = build_protein_db(FASTA_SMALL)

    peptide = "YYRDPCCCRPVSCQTVSRPVTFVPRCTR"
    prot_pos = __get_proteins_and_positions(peptide, protein_db)
    assert prot_pos[0] == ["KRA3_SHEEP"]
    assert prot_pos[1] == [22]


def test6():
    from pyXLMS.transform.reannotate_positions import __get_proteins_and_positions

    protein_db = build_protein_db(FASTA_LARGE)

    peptide = "KEESGF"
    prot_pos = __get_proteins_and_positions(peptide, protein_db)
    assert prot_pos[0] == ["A0A087X1C5", "P10635"]
    assert prot_pos[1] == [213, 213]


def test7():
    from pyXLMS.transform.reannotate_positions import __get_proteins_and_positions

    protein_db = build_protein_db(FASTA_LARGE)

    peptide = "PGTKLCHGDSELTSGLLAT"
    prot_pos = __get_proteins_and_positions(peptide, protein_db)
    assert prot_pos[0] == ["Q9Y6C7"]
    assert prot_pos[1] == [75]


def test8():
    from pyXLMS.transform.reannotate_positions import __get_proteins_and_positions

    protein_db = build_protein_db(FASTA_LARGE)

    peptide = "PGTKLCHGDSELTSGLLAT"
    # this is a performance test to see if it shows up in pytest --durations
    for i in range(1000):
        prot_pos = __get_proteins_and_positions(peptide, protein_db)
        assert prot_pos[0] == ["Q9Y6C7"]
        assert prot_pos[1] == [75]


def test9():
    from pyXLMS.transform import fasta_title_to_accession

    title = "sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1"
    assert fasta_title_to_accession(title) == "A0A087X1C5"
    assert fasta_title_to_accession("Cas9") == "Cas9"


def test10():
    from pyXLMS.transform import reannotate_positions

    with pytest.raises(
        TypeError,
        match="Can't annotate positions for dict. Dict has to be a valid 'parser_result'!",
    ):
        _r = reannotate_positions({"data_type": "error"}, FASTA_SMALL)


def test11():
    from pyXLMS.transform import reannotate_positions
    from pyXLMS.data import create_csm_min

    csms = [
        create_csm_min("ADANLDK", 7, "GNTDRHSIK", 9, "RUN_1", 1),
        create_csm_min("GNTDRHSIK", 9, "YYRDPCCCRPVSCQTVSRPVTFVPRCTR", 1, "RUN_1", 2),
    ]
    new_csms = reannotate_positions(csms, FASTA_SMALL)

    assert new_csms[0]["alpha_proteins"] == ["Cas9"]
    assert new_csms[0]["alpha_proteins_crosslink_positions"] == [1293]
    assert new_csms[0]["alpha_proteins_peptide_positions"] == [1287]
    assert new_csms[0]["beta_proteins"] == ["Cas9"]
    assert new_csms[0]["beta_proteins_crosslink_positions"] == [48]
    assert new_csms[0]["beta_proteins_peptide_positions"] == [40]
    assert new_csms[1]["alpha_proteins"] == ["Cas9"]
    assert new_csms[1]["alpha_proteins_crosslink_positions"] == [48]
    assert new_csms[1]["alpha_proteins_peptide_positions"] == [40]
    assert new_csms[1]["beta_proteins"] == ["KRA3_SHEEP"]
    assert new_csms[1]["beta_proteins_crosslink_positions"] == [23]
    assert new_csms[1]["beta_proteins_peptide_positions"] == [23]


def test12():
    from pyXLMS.transform import reannotate_positions
    from pyXLMS.data import create_crosslink_min

    xls = [
        create_crosslink_min("ADANLDK", 7, "GNTDRHSIK", 9),
        create_crosslink_min("GNTDRHSIK", 9, "YYRDPCCCRPVSCQTVSRPVTFVPRCTR", 1),
    ]
    new_xls = reannotate_positions(xls, FASTA_SMALL)

    assert new_xls[0]["alpha_proteins"] == ["Cas9"]
    assert new_xls[0]["alpha_proteins_crosslink_positions"] == [1293]
    assert new_xls[0]["beta_proteins"] == ["Cas9"]
    assert new_xls[0]["beta_proteins_crosslink_positions"] == [48]
    assert new_xls[1]["alpha_proteins"] == ["Cas9"]
    assert new_xls[1]["alpha_proteins_crosslink_positions"] == [48]
    assert new_xls[1]["beta_proteins"] == ["KRA3_SHEEP"]
    assert new_xls[1]["beta_proteins_crosslink_positions"] == [23]


def test13():
    from pyXLMS.transform import reannotate_positions
    from pyXLMS.data import create_csm_min
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.data import create_parser_result

    csms = [
        create_csm_min("ADANLDK", 7, "GNTDRHSIK", 9, "RUN_1", 1),
        create_csm_min("GNTDRHSIK", 9, "YYRDPCCCRPVSCQTVSRPVTFVPRCTR", 1, "RUN_1", 2),
    ]
    xls = [
        create_crosslink_min("ADANLDK", 7, "GNTDRHSIK", 9),
        create_crosslink_min("GNTDRHSIK", 9, "YYRDPCCCRPVSCQTVSRPVTFVPRCTR", 1),
    ]
    pr = create_parser_result("TEST", csms=csms, crosslinks=xls)
    new_pr = reannotate_positions(pr, FASTA_SMALL)

    assert new_pr["crosslink-spectrum-matches"][0]["alpha_proteins"] == ["Cas9"]
    assert new_pr["crosslink-spectrum-matches"][0][
        "alpha_proteins_crosslink_positions"
    ] == [1293]
    assert new_pr["crosslink-spectrum-matches"][0][
        "alpha_proteins_peptide_positions"
    ] == [1287]
    assert new_pr["crosslink-spectrum-matches"][0]["beta_proteins"] == ["Cas9"]
    assert new_pr["crosslink-spectrum-matches"][0][
        "beta_proteins_crosslink_positions"
    ] == [48]
    assert new_pr["crosslink-spectrum-matches"][0][
        "beta_proteins_peptide_positions"
    ] == [40]
    assert new_pr["crosslink-spectrum-matches"][1]["alpha_proteins"] == ["Cas9"]
    assert new_pr["crosslink-spectrum-matches"][1][
        "alpha_proteins_crosslink_positions"
    ] == [48]
    assert new_pr["crosslink-spectrum-matches"][1][
        "alpha_proteins_peptide_positions"
    ] == [40]
    assert new_pr["crosslink-spectrum-matches"][1]["beta_proteins"] == ["KRA3_SHEEP"]
    assert new_pr["crosslink-spectrum-matches"][1][
        "beta_proteins_crosslink_positions"
    ] == [23]
    assert new_pr["crosslink-spectrum-matches"][1][
        "beta_proteins_peptide_positions"
    ] == [23]
    assert new_pr["crosslinks"][0]["alpha_proteins"] == ["Cas9"]
    assert new_pr["crosslinks"][0]["alpha_proteins_crosslink_positions"] == [1293]
    assert new_pr["crosslinks"][0]["beta_proteins"] == ["Cas9"]
    assert new_pr["crosslinks"][0]["beta_proteins_crosslink_positions"] == [48]
    assert new_pr["crosslinks"][1]["alpha_proteins"] == ["Cas9"]
    assert new_pr["crosslinks"][1]["alpha_proteins_crosslink_positions"] == [48]
    assert new_pr["crosslinks"][1]["beta_proteins"] == ["KRA3_SHEEP"]
    assert new_pr["crosslinks"][1]["beta_proteins_crosslink_positions"] == [23]


def test14():
    from pyXLMS.transform import reannotate_positions
    from pyXLMS.data import create_csm_min

    csms = [create_csm_min("KEESGF", 1, "PGTKLCHGDSELTSGLLAT", 4, "RUN_1", 1)]
    new_csms = reannotate_positions(csms, FASTA_LARGE)

    assert new_csms[0]["alpha_proteins"] == ["A0A087X1C5", "P10635"]
    assert new_csms[0]["alpha_proteins_crosslink_positions"] == [214, 214]
    assert new_csms[0]["alpha_proteins_peptide_positions"] == [214, 214]
    assert new_csms[0]["beta_proteins"] == ["Q9Y6C7"]
    assert new_csms[0]["beta_proteins_crosslink_positions"] == [79]
    assert new_csms[0]["beta_proteins_peptide_positions"] == [76]


def test15():
    from pyXLMS.transform import reannotate_positions
    from pyXLMS.data import create_crosslink_min

    xls = [create_crosslink_min("KEESGF", 1, "PGTKLCHGDSELTSGLLAT", 4)]
    new_xls = reannotate_positions(xls, FASTA_LARGE)

    assert new_xls[0]["alpha_proteins"] == ["A0A087X1C5", "P10635"]
    assert new_xls[0]["alpha_proteins_crosslink_positions"] == [214, 214]
    assert new_xls[0]["beta_proteins"] == ["Q9Y6C7"]
    assert new_xls[0]["beta_proteins_crosslink_positions"] == [79]


def test16():
    from pyXLMS.transform import reannotate_positions
    from pyXLMS.data import create_csm_min
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.data import create_parser_result

    csms = [create_csm_min("KEESGF", 1, "PGTKLCHGDSELTSGLLAT", 4, "RUN_1", 1)]
    xls = [create_crosslink_min("KEESGF", 1, "PGTKLCHGDSELTSGLLAT", 4)]
    pr = create_parser_result("TEST", csms=csms, crosslinks=xls)
    new_pr = reannotate_positions(pr, FASTA_LARGE)

    assert new_pr["crosslink-spectrum-matches"][0]["alpha_proteins"] == [
        "A0A087X1C5",
        "P10635",
    ]
    assert new_pr["crosslink-spectrum-matches"][0][
        "alpha_proteins_crosslink_positions"
    ] == [214, 214]
    assert new_pr["crosslink-spectrum-matches"][0][
        "alpha_proteins_peptide_positions"
    ] == [214, 214]
    assert new_pr["crosslink-spectrum-matches"][0]["beta_proteins"] == ["Q9Y6C7"]
    assert new_pr["crosslink-spectrum-matches"][0][
        "beta_proteins_crosslink_positions"
    ] == [79]
    assert new_pr["crosslink-spectrum-matches"][0][
        "beta_proteins_peptide_positions"
    ] == [76]
    assert new_pr["crosslinks"][0]["alpha_proteins"] == ["A0A087X1C5", "P10635"]
    assert new_pr["crosslinks"][0]["alpha_proteins_crosslink_positions"] == [214, 214]
    assert new_pr["crosslinks"][0]["beta_proteins"] == ["Q9Y6C7"]
    assert new_pr["crosslinks"][0]["beta_proteins_crosslink_positions"] == [79]


def test17():
    from pyXLMS.transform import reannotate_positions
    from pyXLMS.data import create_csm_min
    from pyXLMS.data import create_parser_result

    # this is a performance test to see if it shows up in pytest --durations
    csms = [
        create_csm_min("KEESGF", 1, "PGTKLCHGDSELTSGLLAT", 4, "RUN_1", i + 1)
        for i in range(1000)
    ]
    pr = create_parser_result("TEST", csms=csms, crosslinks=None)
    new_pr = reannotate_positions(pr, FASTA_LARGE)

    assert new_pr["crosslink-spectrum-matches"][0]["alpha_proteins"] == [
        "A0A087X1C5",
        "P10635",
    ]
    assert new_pr["crosslink-spectrum-matches"][0][
        "alpha_proteins_crosslink_positions"
    ] == [214, 214]
    assert new_pr["crosslink-spectrum-matches"][0][
        "alpha_proteins_peptide_positions"
    ] == [214, 214]
    assert new_pr["crosslink-spectrum-matches"][0]["beta_proteins"] == ["Q9Y6C7"]
    assert new_pr["crosslink-spectrum-matches"][0][
        "beta_proteins_crosslink_positions"
    ] == [79]
    assert new_pr["crosslink-spectrum-matches"][0][
        "beta_proteins_peptide_positions"
    ] == [76]
    assert new_pr["crosslink-spectrum-matches"][-1]["alpha_proteins"] == [
        "A0A087X1C5",
        "P10635",
    ]
    assert new_pr["crosslink-spectrum-matches"][-1][
        "alpha_proteins_crosslink_positions"
    ] == [214, 214]
    assert new_pr["crosslink-spectrum-matches"][-1][
        "alpha_proteins_peptide_positions"
    ] == [214, 214]
    assert new_pr["crosslink-spectrum-matches"][-1]["beta_proteins"] == ["Q9Y6C7"]
    assert new_pr["crosslink-spectrum-matches"][-1][
        "beta_proteins_crosslink_positions"
    ] == [79]
    assert new_pr["crosslink-spectrum-matches"][-1][
        "beta_proteins_peptide_positions"
    ] == [76]
    assert new_pr["crosslinks"] is None
