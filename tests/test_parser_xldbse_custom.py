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
XL_REV1 = "data/pyxlms/xl_rev1.txt"
CSM_REV1 = "data/pyxlms/csm_rev1.txt"
XL_REV2 = "data/pyxlms/xl_rev2.txt"
CSM_REV2 = "data/pyxlms/csm_rev2.txt"
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
    assert last_crosslink["completeness"] == "partial"
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


def test6():
    from pyXLMS import parser as p

    parser_result = p.read_custom(XL_REV1)
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


def test7():
    from pyXLMS import parser as p

    parser_result = p.read_custom(XL_REV2)
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
    assert first_crosslink["alpha_decoy"]
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


def test8():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_custom(CSM_MIN)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 2

    first_csm = csms[0]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KPEPTIDE"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] is None
    assert first_csm["alpha_proteins_crosslink_positions"] is None
    assert first_csm["alpha_proteins_peptide_positions"] is None
    assert first_csm["alpha_score"] is None
    assert first_csm["alpha_decoy"] is None
    assert first_csm["beta_peptide"] == "PEPKTIDE"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 4
    assert first_csm["beta_proteins"] is None
    assert first_csm["beta_proteins_crosslink_positions"] is None
    assert first_csm["beta_proteins_peptide_positions"] is None
    assert first_csm["beta_score"] is None
    assert first_csm["beta_decoy"] is None
    assert first_csm["crosslink_type"] == "inter"
    assert first_csm["score"] is None
    assert first_csm["spectrum_file"] == "S1.raw"
    assert first_csm["scan_nr"] == 1
    assert first_csm["charge"] is None
    assert first_csm["retention_time"] is None
    assert first_csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "EKTIDE"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] is None
    assert last_csm["alpha_proteins_crosslink_positions"] is None
    assert last_csm["alpha_proteins_peptide_positions"] is None
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "PEKPIDE"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 3
    assert last_csm["beta_proteins"] is None
    assert last_csm["beta_proteins_crosslink_positions"] is None
    assert last_csm["beta_proteins_peptide_positions"] is None
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] is None
    assert last_csm["spectrum_file"] == "S1.raw"
    assert last_csm["scan_nr"] == 2
    assert last_csm["charge"] is None
    assert last_csm["retention_time"] is None
    assert last_csm["ion_mobility"] is None


def test9():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_custom(CSM)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 2

    first_csm = csms[0]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "KPEPTIDE"
    assert mts(first_csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [13]
    assert first_csm["alpha_proteins_peptide_positions"] == [13]
    assert first_csm["alpha_score"] == pytest.approx(87.53)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "PEPKTIDE"
    assert mts(first_csm["beta_modifications"]) == "(4:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 4
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [17]
    assert first_csm["beta_proteins_peptide_positions"] == [14]
    assert first_csm["beta_score"] == pytest.approx(100.3)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(87.53)
    assert first_csm["spectrum_file"] == "S1.raw"
    assert first_csm["scan_nr"] == 1
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(14.3)
    assert first_csm["ion_mobility"] == pytest.approx(50.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "EKTIDEM"
    assert (
        mts(last_csm["alpha_modifications"])
        == "(2:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    )
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [33, 21]
    assert last_csm["alpha_proteins_peptide_positions"] == [32, 20]
    assert last_csm["alpha_score"] == pytest.approx(5.3)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "PEKPIDE"
    assert mts(last_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 3
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [28]
    assert last_csm["beta_proteins_peptide_positions"] == [26]
    assert last_csm["beta_score"] == pytest.approx(34.89)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(5.4)
    assert last_csm["spectrum_file"] == "S1.raw"
    assert last_csm["scan_nr"] == 2
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(37.332)
    assert last_csm["ion_mobility"] == pytest.approx(-70.0)


def test10():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_custom(CSM_NULL)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 2

    first_csm = csms[0]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KPEPTIDE"
    assert mts(first_csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [13]
    assert first_csm["alpha_proteins_peptide_positions"] == [13]
    assert first_csm["alpha_score"] == pytest.approx(87.53)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "PEPKTIDE"
    assert mts(first_csm["beta_modifications"]) == "(4:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 4
    assert first_csm["beta_proteins"] is None
    assert first_csm["beta_proteins_crosslink_positions"] is None
    assert first_csm["beta_proteins_peptide_positions"] is None
    assert first_csm["beta_score"] == pytest.approx(100.3)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "inter"
    assert first_csm["score"] == pytest.approx(87.53)
    assert first_csm["spectrum_file"] == "S1.raw"
    assert first_csm["scan_nr"] == 1
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(14.3)
    assert first_csm["ion_mobility"] == pytest.approx(50.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "EKTIDEM"
    assert (
        mts(last_csm["alpha_modifications"])
        == "(2:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    )
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [33, 21]
    assert last_csm["alpha_proteins_peptide_positions"] == [32, 20]
    assert last_csm["alpha_score"] == pytest.approx(5.3)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "PEKPIDE"
    assert mts(last_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 3
    assert last_csm["beta_proteins"] is None
    assert last_csm["beta_proteins_crosslink_positions"] is None
    assert last_csm["beta_proteins_peptide_positions"] is None
    assert last_csm["beta_score"] == pytest.approx(34.89)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(5.4)
    assert last_csm["spectrum_file"] == "S1.raw"
    assert last_csm["scan_nr"] == 2
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(37.332)
    assert last_csm["ion_mobility"] == pytest.approx(-70.0)


def test11():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_custom(
        CSM_FORMAT, column_mapping={"Sequence A": "Alpha Peptide"}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 2

    first_csm = csms[0]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "KPEPTIDE"
    assert mts(first_csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [13]
    assert first_csm["alpha_proteins_peptide_positions"] == [13]
    assert first_csm["alpha_score"] == pytest.approx(87.53)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "PEPKTIDE"
    assert mts(first_csm["beta_modifications"]) == "(4:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 4
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [17]
    assert first_csm["beta_proteins_peptide_positions"] == [14]
    assert first_csm["beta_score"] == pytest.approx(100.3)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(87.53)
    assert first_csm["spectrum_file"] == "S1.raw"
    assert first_csm["scan_nr"] == 1
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(14.3)
    assert first_csm["ion_mobility"] == pytest.approx(50.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "EKTIDEM"
    assert (
        mts(last_csm["alpha_modifications"])
        == "(2:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    )
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [33, 21]
    assert last_csm["alpha_proteins_peptide_positions"] == [32, 20]
    assert last_csm["alpha_score"] == pytest.approx(5.3)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "PEKPIDE"
    assert mts(last_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 3
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [28]
    assert last_csm["beta_proteins_peptide_positions"] == [26]
    assert last_csm["beta_score"] == pytest.approx(34.89)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(5.4)
    assert last_csm["spectrum_file"] == "S1.raw"
    assert last_csm["scan_nr"] == 2
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(37.332)
    assert last_csm["ion_mobility"] == pytest.approx(-70.0)


def test12():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_custom(CSM_REV1)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 2

    first_csm = csms[0]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "KPEPTIDE"
    assert mts(first_csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [13]
    assert first_csm["alpha_proteins_peptide_positions"] == [13]
    assert first_csm["alpha_score"] == pytest.approx(87.53)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "PEPKTIDE"
    assert mts(first_csm["beta_modifications"]) == "(4:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 4
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [17]
    assert first_csm["beta_proteins_peptide_positions"] == [14]
    assert first_csm["beta_score"] == pytest.approx(100.3)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(87.53)
    assert first_csm["spectrum_file"] == "S1.raw"
    assert first_csm["scan_nr"] == 1
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(14.3)
    assert first_csm["ion_mobility"] == pytest.approx(50.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "EKTIDEM"
    assert (
        mts(last_csm["alpha_modifications"])
        == "(2:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    )
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [33, 21]
    assert last_csm["alpha_proteins_peptide_positions"] == [32, 20]
    assert last_csm["alpha_score"] == pytest.approx(5.3)
    assert not last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "PEKPIDE"
    assert mts(last_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 3
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [28]
    assert last_csm["beta_proteins_peptide_positions"] == [26]
    assert last_csm["beta_score"] == pytest.approx(34.89)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(5.4)
    assert last_csm["spectrum_file"] == "S1.raw"
    assert last_csm["scan_nr"] == 2
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(37.332)
    assert last_csm["ion_mobility"] == pytest.approx(-70.0)


def test13():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_custom(CSM_REV2)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 2

    first_csm = csms[0]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "KPEPTIDE"
    assert mts(first_csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [13]
    assert first_csm["alpha_proteins_peptide_positions"] == [13]
    assert first_csm["alpha_score"] == pytest.approx(87.53)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "PEPKTIDE"
    assert mts(first_csm["beta_modifications"]) == "(4:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 4
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [17]
    assert first_csm["beta_proteins_peptide_positions"] == [14]
    assert first_csm["beta_score"] == pytest.approx(100.3)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(87.53)
    assert first_csm["spectrum_file"] == "S1.raw"
    assert first_csm["scan_nr"] == 1
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(14.3)
    assert first_csm["ion_mobility"] == pytest.approx(50.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "EKTIDEM"
    assert (
        mts(last_csm["alpha_modifications"])
        == "(2:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    )
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [33, 21]
    assert last_csm["alpha_proteins_peptide_positions"] == [32, 20]
    assert last_csm["alpha_score"] == pytest.approx(5.3)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "PEKPIDE"
    assert mts(last_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 3
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [28]
    assert last_csm["beta_proteins_peptide_positions"] == [26]
    assert last_csm["beta_score"] == pytest.approx(34.89)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(5.4)
    assert last_csm["spectrum_file"] == "S1.raw"
    assert last_csm["scan_nr"] == 2
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(37.332)
    assert last_csm["ion_mobility"] == pytest.approx(-70.0)


def test14():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_custom([CSM, XL])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
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

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 2

    first_csm = csms[0]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "KPEPTIDE"
    assert mts(first_csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [13]
    assert first_csm["alpha_proteins_peptide_positions"] == [13]
    assert first_csm["alpha_score"] == pytest.approx(87.53)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "PEPKTIDE"
    assert mts(first_csm["beta_modifications"]) == "(4:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 4
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [17]
    assert first_csm["beta_proteins_peptide_positions"] == [14]
    assert first_csm["beta_score"] == pytest.approx(100.3)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(87.53)
    assert first_csm["spectrum_file"] == "S1.raw"
    assert first_csm["scan_nr"] == 1
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(14.3)
    assert first_csm["ion_mobility"] == pytest.approx(50.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "EKTIDEM"
    assert (
        mts(last_csm["alpha_modifications"])
        == "(2:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    )
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] == ["Cas10", "Cas11"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [33, 21]
    assert last_csm["alpha_proteins_peptide_positions"] == [32, 20]
    assert last_csm["alpha_score"] == pytest.approx(5.3)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "PEKPIDE"
    assert mts(last_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 3
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [28]
    assert last_csm["beta_proteins_peptide_positions"] == [26]
    assert last_csm["beta_score"] == pytest.approx(34.89)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(5.4)
    assert last_csm["spectrum_file"] == "S1.raw"
    assert last_csm["scan_nr"] == 2
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(37.332)
    assert last_csm["ion_mobility"] == pytest.approx(-70.0)


def test15():
    from pyXLMS import parser as p
    from pyXLMS.transform import to_dataframe
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika([ANNIKA_XL, ANNIKA_CSM])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300
    to_dataframe(crosslinks).to_csv("XL.csv", index=False)

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826
    to_dataframe(csms).to_csv("CSM.csv", index=False)

    parser_result = p.read_custom(["XL.csv", "CSM.csv"])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "Custom"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "full"
    assert csm["alpha_peptide"] == "KIECFDSVEISGVEDR"
    assert (
        mts(csm["alpha_modifications"])
        == "(1:[DSS|138.06808]);(4:[Carbamidomethyl|57.021464])"
    )
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [575]
    assert csm["alpha_proteins_peptide_positions"] == [575]
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)
