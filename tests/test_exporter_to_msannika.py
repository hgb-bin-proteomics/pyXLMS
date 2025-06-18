#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

XIFDR_CSMS = "data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv"
XIFDR_LINKS = "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv"


def test1():
    from pyXLMS.exporter import get_msannika_crosslink_sequence

    assert get_msannika_crosslink_sequence("PEPKTIDE", 4) == "PEP[K]TIDE"
    assert get_msannika_crosslink_sequence("KPEPTIDE", 1) == "[K]PEPTIDE"
    assert get_msannika_crosslink_sequence("PEPTIDEK", 8) == "PEPTIDE[K]"
    with pytest.raises(
        ValueError,
        match=r"Crosslink position outside of range! Must be in range \[1, 8\]\.",
    ):
        _seq = get_msannika_crosslink_sequence("PEPTIDEK", 9)


def test2():
    from pyXLMS.exporter import to_msannika
    from pyXLMS.data import create_crosslink_min

    with pytest.raises(
        TypeError,
        match="Parameter 'format' has to be one of 'csv', 'tsv', or 'xlsx'!",
    ):
        xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
        xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
        crosslinks = [xl1, xl2]
        _df = to_msannika(crosslinks, "crosslinks.parquet", format="parquet")


def test3():
    from pyXLMS.exporter import to_msannika

    with pytest.raises(
        ValueError,
        match="Provided data does not contain any crosslinks or crosslink-spectrum-matches!",
    ):
        _df = to_msannika([])


def test4():
    from pyXLMS.exporter import to_msannika
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]
    df = to_msannika(crosslinks)
    assert df.shape[0] == 2
    assert df.shape[1] == 11


def test5():
    from pyXLMS.exporter import to_msannika
    from pyXLMS.data import create_crosslink_min

    xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    crosslinks = [xl1, xl2]
    df = to_msannika(crosslinks, filename="crosslinks.csv", format="csv")
    assert df.shape[0] == 2
    assert df.shape[1] == 11


def test6():
    from pyXLMS.exporter import to_msannika
    from pyXLMS.data import create_csm_min

    csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    csm2 = create_csm_min("PEKPTIDE", 3, "PEPKTIDE", 4, "RUN_1", 2)
    csms = [csm1, csm2]
    df = to_msannika(csms)
    assert df.shape[0] == 2
    assert df.shape[1] == 20


def test7():
    from pyXLMS.exporter import to_msannika
    from pyXLMS.data import create_csm_min

    csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    csm2 = create_csm_min("PEKPTIDE", 3, "PEPKTIDE", 4, "RUN_1", 2)
    csms = [csm1, csm2]
    df = to_msannika(csms, filename="csms.csv", format="csv")
    assert df.shape[0] == 2
    assert df.shape[1] == 20


def test8():
    from pyXLMS import parser as p
    from pyXLMS import exporter as e

    xi = p.read_xi(XIFDR_CSMS, parse_modifications=False, modifications={}, verbose=0)
    _ = e.to_msannika(
        xi["crosslink-spectrum-matches"], filename="xi_csms_to_annika.csv", format="csv"
    )
    with pytest.raises(KeyError, match="Modifications A"):
        _pr = p.read_msannika("xi_csms_to_annika.csv", sep=",")


def test9():
    from pyXLMS import parser as p
    from pyXLMS import exporter as e
    from pyXLMS.transform import modifications_to_str as mts

    xi = p.read_xi(XIFDR_CSMS, verbose=0)
    _ = e.to_msannika(
        xi["crosslink-spectrum-matches"], filename="xi_csms_to_annika.csv", format="csv"
    )
    pr = p.read_msannika("xi_csms_to_annika.csv", sep=",", parse_modifications=False)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MS Annika"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 413

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KIECFDSVEISGVEDR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [575]
    assert csm["alpha_proteins_peptide_positions"] == [575]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KIECFDSVEISGVEDR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [575]
    assert csm["beta_proteins_peptide_positions"] == [575]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(27.268)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 19140
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[332]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NPIDFLEAKGYK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 9
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [1189]
    assert csm["alpha_proteins_peptide_positions"] == [1181]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "QWYKNK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 4
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [488]
    assert csm["beta_proteins_peptide_positions"] == [485]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(10.071)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 21757
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "QITKHVAQILDSR"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 4
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [933]
    assert csm["alpha_proteins_peptide_positions"] == [930]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "QQLPEKYK"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [350]
    assert csm["beta_proteins_peptide_positions"] == [345]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(6.808)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.mgf"
    assert csm["scan_nr"] == 12190
    assert csm["charge"] == 5
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test10():
    from pyXLMS import parser as p
    from pyXLMS import exporter as e

    xi = p.read_xi(XIFDR_LINKS, parse_modifications=False, modifications={}, verbose=0)
    _ = e.to_msannika(xi["crosslinks"], filename="xi_xls_to_annika.xlsx", format="xlsx")
    pr = p.read_msannika("xi_xls_to_annika.xlsx", parse_modifications=False)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MS Annika"
    assert pr["crosslink-spectrum-matches"] is None
    assert pr["crosslinks"] is not None

    xls = pr["crosslinks"]
    assert len(xls) == 227

    xl = xls[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "VVDELVKVMGR"
    assert xl["alpha_peptide_crosslink_position"] == 7
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [753]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "VVDELVKVMGR"
    assert xl["beta_peptide_crosslink_position"] == 7
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [753]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(40.679)

    xl = xls[30]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "FDNLTKAER"
    assert xl["alpha_peptide_crosslink_position"] == 6
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [952]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "YDENDKLIR"
    assert xl["beta_peptide_crosslink_position"] == 6
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [906]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(25.616)

    xl = xls[219]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "NPIDFLEAKGYK"
    assert xl["alpha_peptide_crosslink_position"] == 9
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [488]
    assert xl["alpha_decoy"]
    assert xl["beta_peptide"] == "QWYKNK"
    assert xl["beta_peptide_crosslink_position"] == 4
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [1189]
    assert xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(10.071)

    xl = xls[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "LSKSR"
    assert xl["alpha_peptide_crosslink_position"] == 3
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [222]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "LSKSR"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [222]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(9.619)
