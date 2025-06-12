#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


XL_TSV = "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_Crosslinks.txt"
XL_XLSX = "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_Crosslinks.xlsx"
CSMS_TSV = "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_CSMs.txt"
CSMS_XLSX = "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_CSMs.xlsx"


def test1():
    from pyXLMS import parser as p

    parser_result = p.read_xlinkx(XL_TSV, parse_modifications=False, modifications={})
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 236

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)


def test2():
    from pyXLMS import parser as p

    parser_result = p.read_xlinkx(XL_XLSX, parse_modifications=False, modifications={})
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 236

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)


def test3():
    from pyXLMS import parser as p

    parser_result = p.read_xlinkx(
        [XL_TSV, XL_XLSX], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 472

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)


def test4():
    from pyXLMS import parser as p

    parser_result = p.read_xlinkx(
        [XL_TSV, XL_TSV], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 472

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)


def test5():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(CSMS_TSV, parse_modifications=False, modifications={})
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None


def test6():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(
        CSMS_XLSX, parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None


def test7():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(
        [CSMS_TSV, CSMS_XLSX], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823 * 2

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None


def test8():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(
        [CSMS_TSV, CSMS_TSV], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823 * 2

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None


def test9():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(
        [XL_XLSX, CSMS_XLSX], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 236

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None


def test10():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(
        [XL_XLSX, CSMS_TSV], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 236

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None


def test11():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(
        [XL_TSV, CSMS_XLSX], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 236

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None


def test12():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(
        [XL_TSV, CSMS_TSV], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 236

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None


def test13():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_xlinkx(
        [XL_XLSX, CSMS_XLSX, CSMS_TSV], parse_modifications=False, modifications={}
    )
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "XlinkX"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 236

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[-1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 1
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [49]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "LVDSTDKADLR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 7
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [152]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(445.6)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "LESEFVYGDYKVYDVR"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["Cas9"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [1018]
    assert not last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "VVDELVKVMGR"
    assert last_crosslink["beta_peptide_crosslink_position"] == 7
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [753]
    assert not last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "intra"
    assert last_crosslink["score"] == pytest.approx(11.49)

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 823 * 2

    first_csm = csms[0]
    csm = csms[21]
    last_csm = csms[-1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "partial"
    assert first_csm["alpha_peptide"] == "KNLIGALLFDSGETAEATR"
    assert mts(first_csm["alpha_modifications"]) is None
    assert first_csm["alpha_peptide_crosslink_position"] == 1
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [49]
    assert first_csm["alpha_proteins_peptide_positions"] == [49]
    assert first_csm["alpha_score"] is None
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "LVDSTDKADLR"
    assert mts(first_csm["beta_modifications"]) is None
    assert first_csm["beta_peptide_crosslink_position"] == 7
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [152]
    assert first_csm["beta_proteins_peptide_positions"] == [146]
    assert first_csm["beta_score"] is None
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(445.6)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert first_csm["scan_nr"] == 43526
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(88.612 * 60.0)
    assert first_csm["ion_mobility"] is None

    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "NLGKVGSKCCK"
    assert mts(csm["alpha_modifications"]) is None
    assert csm["alpha_peptide_crosslink_position"] == 8
    assert csm["alpha_proteins"] == ["spALBU_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [460]
    assert csm["alpha_proteins_peptide_positions"] == [453]
    assert csm["alpha_score"] is None
    assert csm["alpha_decoy"] is None
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) is None
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] is None
    assert csm["beta_decoy"] is None
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(297.32)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert csm["scan_nr"] == 46493
    assert csm["charge"] == 3
    assert csm["retention_time"] == pytest.approx(94.285 * 60.0)
    assert csm["ion_mobility"] is None

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "partial"
    assert last_csm["alpha_peptide"] == "ESILPKR"
    assert mts(last_csm["alpha_modifications"]) is None
    assert last_csm["alpha_peptide_crosslink_position"] == 6
    assert last_csm["alpha_proteins"] == ["Cas9"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [1117]
    assert last_csm["alpha_proteins_peptide_positions"] == [1112]
    assert last_csm["alpha_score"] is None
    assert last_csm["alpha_decoy"] is None
    assert last_csm["beta_peptide"] == "KMIAK"
    assert mts(last_csm["beta_modifications"]) is None
    assert last_csm["beta_peptide_crosslink_position"] == 1
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1024]
    assert last_csm["beta_proteins_peptide_positions"] == [1024]
    assert last_csm["beta_score"] is None
    assert last_csm["beta_decoy"] is None
    assert last_csm["crosslink_type"] == "intra"
    assert last_csm["score"] == pytest.approx(7.45)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw"
    assert last_csm["scan_nr"] == 13466
    assert last_csm["charge"] == 4
    assert last_csm["retention_time"] == pytest.approx(36.2167 * 60.0)
    assert last_csm["ion_mobility"] is None
