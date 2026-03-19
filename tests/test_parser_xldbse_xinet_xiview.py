#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


MSANNIKA = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt"
XIVIEW = "data/xiview/DDX39B_LCSDA_shared_links_open_clamped.csv"


def test1():
    from pyXLMS import parser as p
    from pyXLMS import exporter as e

    parser_result = p.read_msannika(
        MSANNIKA, parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    e.to_xinet(crosslinks, filename="export_toxinet.csv")
    parser_result = p.read_xinet("export_toxinet.csv")
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "xiNET/xiVIEW"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "partial"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert first_crosslink["alpha_decoy"] is None
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert first_crosslink["beta_decoy"] is None
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "partial"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"] is None
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"] is None
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)


def test2():
    from pyXLMS import parser as p
    from pyXLMS import exporter as e

    parser_result = p.read_msannika(
        MSANNIKA, parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    e.to_xiview(crosslinks, filename="export_toxiview.csv", minimal=False)
    parser_result = p.read_xiview("export_toxiview.csv")
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "xiNET/xiVIEW"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "partial"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert first_crosslink["alpha_decoy"] is None
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert first_crosslink["beta_decoy"] is None
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "partial"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"] is None
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"] is None
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)


def test3():
    from pyXLMS import parser as p
    from pyXLMS import exporter as e

    parser_result = p.read_msannika(
        MSANNIKA, parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    e.to_xiview(crosslinks, filename="export_toxiview_min.csv", minimal=True)

    with pytest.raises(
        KeyError, match="Could not get a suitable column for the peptide sequence!"
    ):
        _parser_result = p.read_xiview("export_toxiview_min.csv")


def test4():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xiview(XIVIEW)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "xiNET/xiVIEW"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 124

    first_csm = csms[0]
    last_csm = csms[123]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "MTPVGTASNVKAQAAKEAQHAQLVAVAEDK"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["DECOY_decoy:P11940"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [1]
    assert first_csm["alpha_proteins_peptide_positions"] == [1]
    assert first_csm["alpha_score"] is None
    assert first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "TVAPTAAAAAAARPPGMTQTSTNAR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 1
    assert first_csm["beta_proteins"] == ["DECOY_decoy:P11940"]
    assert first_csm["beta_proteins_crosslink_positions"] == [131]
    assert first_csm["beta_proteins_peptide_positions"] == [131]
    assert first_csm["beta_score"] is None
    assert first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(20.280395)
    assert first_csm["spectrum_file"] == ""
    assert first_csm["scan_nr"] == 20
    assert first_csm["charge"] == 6
    assert first_csm["retention_time"] is None
    assert first_csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "GSGGGSSSGR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 2
    assert last_csm["alpha_proteins"] == ["DECOY_decoy:cont_Q3TTY5"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [51]
    assert last_csm["alpha_proteins_peptide_positions"] == [50]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "KDVDSCYMDK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["cont_Q3TTY5"]
    assert last_csm["beta_proteins_crosslink_positions"] == [308]
    assert last_csm["beta_proteins_peptide_positions"] == [308]
    assert last_csm["beta_score"] is None
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(6.280929)
    assert last_csm["spectrum_file"] == ""
    assert last_csm["scan_nr"] == 364
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] is None
    assert last_csm["ion_mobility"] is None


def test5():
    from pyXLMS import parser as p

    with pytest.raises(
        KeyError,
        match="Could not get a suitable column or value for the spectrum file name!",
    ):
        _parser_result = p.read_xiview(XIVIEW, verbose=2)


def test6():
    from pyXLMS import parser as p

    with pytest.raises(TypeError, match="Verbose level has to be one of 0, 1, or 2!"):
        _parser_result = p.read_xiview(XIVIEW, verbose=3)
