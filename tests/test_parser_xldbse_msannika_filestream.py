#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


XL_TSV = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt"
XL_XLSX = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx"
CSMS_TSV = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt"
CSMS_XLSX = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"


def test1():
    from pyXLMS import parser as p

    with open(XL_TSV, "r", encoding="utf-8") as f:
        parser_result = p.read_msannika(f, format="tsv")
        f.close()

    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

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


def test2():
    from pyXLMS import parser as p

    with open(XL_XLSX, "rb") as f:
        parser_result = p.read_msannika(f, format="xlsx")
        f.close()

    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

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


def test3():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    with open(CSMS_TSV, "r", encoding="utf-8") as f:
        parser_result = p.read_msannika(f, format="tsv")
        f.close()

    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[825]

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


def test4():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    with open(CSMS_XLSX, "rb") as f:
        parser_result = p.read_msannika(f, format="xlsx")
        f.close()

    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[825]

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
