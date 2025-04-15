#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


XL = "data/pyxlms/xl.txt"
CSM = "data/pyxlms/csm.txt"
XL_MIN = "data/pyxlms/xl_min.txt"
CSM_MIN = "data/pyxlms/csm_min.txt"
XL_NULL = "data/pyxlms/xl_null.txt"
CSM_NULL = "data/pyxlms/csm_null.txt"
XL_FORMAT = "data/pyxlms/xl_format.txt"
CSM_FORMAT = "data/pyxlms/csm_format.txt"
ANNIKA_XL = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt"
ANNIKA_CSM = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"


def test1():
    from pyXLMS.parser import pyxlms_modification_str_parser as mp

    assert mp("(1:[DSS|138.06808])") == {1: ("DSS", 138.06808)}
    assert mp("(1:[DSS|138.06808]);(7:[Oxidation|15.994915])") == {
        1: ("DSS", 138.06808),
        7: ("Oxidation", 15.994915),
    }


def test2():
    from pyXLMS import parser as p

    parser_result = p.read_custom(XL_MIN)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 2

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "partial"
    assert first_crosslink["alpha_peptide"] == "KPEPTIDE"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] is None
    assert first_crosslink["alpha_proteins_crosslink_positions"] is None
    assert first_crosslink["alpha_decoy"] is None
    assert first_crosslink["beta_peptide"] == "PEPKTIDE"
    assert first_crosslink["beta_peptide_crosslink_position"] == 4
    assert first_crosslink["beta_proteins"] is None
    assert first_crosslink["beta_proteins_crosslink_positions"] is None
    assert first_crosslink["beta_decoy"] is None
    assert first_crosslink["crosslink_type"] == "inter"
    assert first_crosslink["score"] is None

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "partial"
    assert last_crosslink["alpha_peptide"] == "EKTIDE"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 2
    assert last_crosslink["alpha_proteins"] is None
    assert last_crosslink["alpha_proteins_crosslink_positions"] is None
    assert last_crosslink["alpha_decoy"] is None
    assert last_crosslink["beta_peptide"] == "PEKPIDE"
    assert last_crosslink["beta_peptide_crosslink_position"] == 3
    assert last_crosslink["beta_proteins"] is None
    assert last_crosslink["beta_proteins_crosslink_positions"] is None
    assert last_crosslink["beta_decoy"] is None
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] is None


def test3():
    from pyXLMS import parser as p

    parser_result = p.read_custom(XL)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 2

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KPEPTIDE"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [11]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "PEPKTIDE"
    assert first_crosslink["beta_peptide_crosslink_position"] == 4
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [15]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(100.3)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "EKTIDE"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 2
    assert last_crosslink["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [11, 13]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "PEKPIDE"
    assert last_crosslink["beta_peptide_crosslink_position"] == 3
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [3]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(3.14159)


def test4():
    from pyXLMS import parser as p

    parser_result = p.read_custom(XL_NULL)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 2

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "partial"
    assert first_crosslink["alpha_peptide"] == "KPEPTIDE"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [11]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "PEPKTIDE"
    assert first_crosslink["beta_peptide_crosslink_position"] == 4
    assert first_crosslink["beta_proteins"] is None
    assert first_crosslink["beta_proteins_crosslink_positions"] is None
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "inter"
    assert first_crosslink["score"] == pytest.approx(100.3)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "EKTIDE"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 2
    assert last_crosslink["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [11, 13]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "PEKPIDE"
    assert last_crosslink["beta_peptide_crosslink_position"] == 3
    assert last_crosslink["beta_proteins"] is None
    assert last_crosslink["beta_proteins_crosslink_positions"] is None
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(3.14159)


def test5():
    from pyXLMS import parser as p

    parser_result = p.read_custom(XL, column_mapping={"Sequence A": "Alpha Peptide"})
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 2

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KPEPTIDE"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [11]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "PEPKTIDE"
    assert first_crosslink["beta_peptide_crosslink_position"] == 4
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [15]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(100.3)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "EKTIDE"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 2
    assert last_crosslink["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [11, 13]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "PEKPIDE"
    assert last_crosslink["beta_peptide_crosslink_position"] == 3
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [3]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(3.14159)
