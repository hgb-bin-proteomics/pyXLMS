#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

PLINK2_CSMS = "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv"
PLINK3_CSMS = "data/plink3/Cas10_plus10_2025.04.07.filtered_cross-linked_spectra.csv"
PLINK2_XLS = "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_peptides.csv"
PLINK3_XLS = "data/plink3/Cas10_plus10_2025.04.07.filtered_cross-linked_peptides.csv"


def test1():
    from pyXLMS.parser import detect_plink_filetype

    assert (
        detect_plink_filetype(
            "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_peptides.csv"
        )
        == "crosslinks"
    )


def test2():
    from pyXLMS.parser import detect_plink_filetype

    assert (
        detect_plink_filetype(
            "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv"
        )
        == "crosslink-spectrum-matches"
    )


def test3():
    from pyXLMS.parser import detect_plink_filetype

    assert detect_plink_filetype(PLINK2_CSMS) == "crosslink-spectrum-matches"
    assert detect_plink_filetype(PLINK3_CSMS) == "crosslink-spectrum-matches"
    assert detect_plink_filetype(PLINK2_XLS) == "crosslinks"
    assert detect_plink_filetype(PLINK3_XLS) == "crosslinks"
    with open(PLINK2_CSMS, "r", encoding="utf-8") as f:
        assert detect_plink_filetype(f) == "crosslink-spectrum-matches"
    with open(PLINK3_CSMS, "r", encoding="utf-8") as f:
        assert detect_plink_filetype(f) == "crosslink-spectrum-matches"
    with open(PLINK2_XLS, "r", encoding="utf-8") as f:
        assert detect_plink_filetype(f) == "crosslinks"
    with open(PLINK3_XLS, "r", encoding="utf-8") as f:
        assert detect_plink_filetype(f) == "crosslinks"


def test4():
    from pyXLMS.parser import detect_plink_filetype

    with pytest.raises(
        ValueError,
        match="The provided file seems not to be one of the support pLink input files!",
    ):
        _ = detect_plink_filetype(
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt"
        )


def test5():
    from pyXLMS.parser.parser_xldbse_plink import (
        __read_plink_cross_linked_peptides_file,
    )

    plink2 = __read_plink_cross_linked_peptides_file(PLINK2_XLS)
    plink3 = __read_plink_cross_linked_peptides_file(PLINK3_XLS)
    assert plink2 is not None
    assert plink3 is not None


def test6():
    from pyXLMS.parser.parser_xldbse_plink import (
        __read_plink_cross_linked_peptides_file,
    )

    with open(PLINK2_XLS, "r", encoding="utf-8") as f:
        plink2 = __read_plink_cross_linked_peptides_file(f)
        f.close()
    with open(PLINK3_XLS, "r", encoding="utf-8") as f:
        plink3 = __read_plink_cross_linked_peptides_file(f)
        f.close()
    assert plink2 is not None
    assert plink3 is not None


def test7():
    from pyXLMS.parser.parser_xldbse_plink import (
        __read_plink_cross_linked_peptides_file,
    )

    with pytest.raises(
        ValueError,
        match="The provided file seems not to be a pLink cross-linked peptides file!",
    ):
        _ = __read_plink_cross_linked_peptides_file(PLINK2_CSMS)


def test8():
    from pyXLMS.parser import read_plink

    crosslinks = read_plink(
        "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_peptides.csv"
    )
    assert len(crosslinks["crosslinks"]) > 0


def test9():
    from pyXLMS.parser import read_plink

    pr = read_plink(PLINK2_XLS)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "pLink"
    assert pr["crosslink-spectrum-matches"] is None
    assert pr["crosslinks"] is not None

    xls = pr["crosslinks"]
    assert len(xls) == 252

    xl = xls[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "AGFIKR"
    assert xl["alpha_peptide_crosslink_position"] == 5
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [922]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "AGFIKR"
    assert xl["beta_peptide_crosslink_position"] == 5
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [922]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None

    xl = xls[1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "AGFIKR"
    assert xl["alpha_peptide_crosslink_position"] == 5
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [922]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "GQKNSR"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [779]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None

    xl = xls[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "LKSVK"
    assert xl["alpha_peptide_crosslink_position"] == 2
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [1162]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "YPKLESEFVYGDYK"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [1007]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None


def test10():
    from pyXLMS.parser import read_plink

    with open(PLINK2_XLS, "r", encoding="utf-8") as f:
        pr = read_plink(f)
        f.close()

    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "pLink"
    assert pr["crosslink-spectrum-matches"] is None
    assert pr["crosslinks"] is not None

    xls = pr["crosslinks"]
    assert len(xls) == 252

    xl = xls[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "AGFIKR"
    assert xl["alpha_peptide_crosslink_position"] == 5
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [922]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "AGFIKR"
    assert xl["beta_peptide_crosslink_position"] == 5
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [922]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None

    xl = xls[1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "AGFIKR"
    assert xl["alpha_peptide_crosslink_position"] == 5
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [922]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "GQKNSR"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [779]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None

    xl = xls[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "LKSVK"
    assert xl["alpha_peptide_crosslink_position"] == 2
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [1162]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "YPKLESEFVYGDYK"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [1007]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None


def test11():
    from pyXLMS.parser import read_plink

    pr = read_plink(PLINK3_XLS)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "pLink"
    assert pr["crosslink-spectrum-matches"] is None
    assert pr["crosslinks"] is not None

    xls = pr["crosslinks"]
    assert len(xls) == 2

    xl = xls[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "LFKGHPETLEK"
    assert xl["alpha_peptide_crosslink_position"] == 3
    assert xl["alpha_proteins"] == ["sp|MYG_HUMAN|"]
    assert xl["alpha_proteins_crosslink_positions"] == [34]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "LFKGHPETLEK"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["sp|MYG_HUMAN|"]
    assert xl["beta_proteins_crosslink_positions"] == [34]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None

    xl = xls[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "LKYENEMALR"
    assert xl["alpha_peptide_crosslink_position"] == 2
    assert xl["alpha_proteins"] == ["sp|K1C15_SHEEP|"]
    assert xl["alpha_proteins_crosslink_positions"] == [192]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "LKYENEMALR"
    assert xl["beta_peptide_crosslink_position"] == 2
    assert xl["beta_proteins"] == ["sp|K1C15_SHEEP|"]
    assert xl["beta_proteins_crosslink_positions"] == [192]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None


def test12():
    from pyXLMS.parser import read_plink

    with open(PLINK3_XLS, "r", encoding="utf-8") as f:
        pr = read_plink(f)
        f.close()

    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "pLink"
    assert pr["crosslink-spectrum-matches"] is None
    assert pr["crosslinks"] is not None

    xls = pr["crosslinks"]
    assert len(xls) == 2

    xl = xls[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "LFKGHPETLEK"
    assert xl["alpha_peptide_crosslink_position"] == 3
    assert xl["alpha_proteins"] == ["sp|MYG_HUMAN|"]
    assert xl["alpha_proteins_crosslink_positions"] == [34]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "LFKGHPETLEK"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["sp|MYG_HUMAN|"]
    assert xl["beta_proteins_crosslink_positions"] == [34]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None

    xl = xls[-1]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "partial"
    assert xl["alpha_peptide"] == "LKYENEMALR"
    assert xl["alpha_peptide_crosslink_position"] == 2
    assert xl["alpha_proteins"] == ["sp|K1C15_SHEEP|"]
    assert xl["alpha_proteins_crosslink_positions"] == [192]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "LKYENEMALR"
    assert xl["beta_peptide_crosslink_position"] == 2
    assert xl["beta_proteins"] == ["sp|K1C15_SHEEP|"]
    assert xl["beta_proteins_crosslink_positions"] == [192]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] is None


def test13():
    from pyXLMS.parser import read_plink
    from pyXLMS.transform import modifications_to_str as mts

    with open(PLINK2_CSMS, "r", encoding="utf-8") as f:
        pr = read_plink(f)
        f.close()

    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "pLink"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 961

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "FDNLTKAER"
    assert mts(csm["alpha_modifications"]) == "(6:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 6
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [906]
    assert csm["alpha_proteins_peptide_positions"] == [901]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "YDENDKLIR"
    assert mts(csm["beta_modifications"]) == "(6:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [952]
    assert csm["beta_proteins_peptide_positions"] == [947]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(float("5.553153e-009"))
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 13098
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[118]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KIECFDSVEISGVEDR"
    assert (
        mts(csm["alpha_modifications"])
        == "(1:[DSS|138.06808]);(4:[Carbamidomethyl|57.021464])"
    )
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [575]
    assert csm["alpha_proteins_peptide_positions"] == [575]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KIECFDSVEISGVEDR"
    assert (
        mts(csm["beta_modifications"])
        == "(1:[DSS|138.06808]);(4:[Carbamidomethyl|57.021464])"
    )
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [575]
    assert csm["beta_proteins_peptide_positions"] == [575]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(float("8.010787e-004"))
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 19140
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[951]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GKSDNVPSEEVVK"
    assert mts(csm["alpha_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [870]
    assert csm["alpha_proteins_peptide_positions"] == [869]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "VKYVTEGMR"
    assert (
        mts(csm["beta_modifications"])
        == "(2:[DSS|138.06808]);(8:[Oxidation|15.994915])"
    )
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [532]
    assert csm["beta_proteins_peptide_positions"] == [531]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(float("4.356499e-001"))
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 11060
    assert csm["charge"] == 4
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "LKYENEMALR"
    assert mts(csm["alpha_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["sp|K1C15_SHEEP|"]
    assert csm["alpha_proteins_crosslink_positions"] == [192]
    assert csm["alpha_proteins_peptide_positions"] == [191]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TKYNALK"
    assert mts(csm["beta_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["sp|CTRB_BOVIN|"]
    assert csm["beta_proteins_crosslink_positions"] == [145]
    assert csm["beta_proteins_peptide_positions"] == [144]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "inter"
    assert csm["score"] == pytest.approx(float("4.868460e-001"))
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 11341
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None


def test14():
    from pyXLMS.parser import read_plink
    from pyXLMS.transform import modifications_to_str as mts

    with open(PLINK3_CSMS, "r", encoding="utf-8") as f:
        pr = read_plink(f)
        f.close()

    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "pLink"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 3

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "LFKGHPETLEK"
    assert mts(csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 3
    assert csm["alpha_proteins"] == ["sp|K1C15_SHEEP|", "sp|MYG_HUMAN|"]
    assert csm["alpha_proteins_crosslink_positions"] == [192, 34]
    assert csm["alpha_proteins_peptide_positions"] == [190, 32]
    assert csm["alpha_score"] is None
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "LKYENEMALR"
    assert (
        mts(csm["beta_modifications"])
        == "(2:[DSS|138.06808]);(7:[Oxidation|15.994915])"
    )
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["sp|K1C15_SHEEP|", "sp|MYG_HUMAN|"]
    assert csm["beta_proteins_crosslink_positions"] == [192, 34]
    assert csm["beta_proteins_peptide_positions"] == [191, 33]
    assert csm["beta_score"] is None
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(float("0.332623"))
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 14421
    assert csm["charge"] == 3
    assert csm["retention_time"] is None
    assert csm["ion_mobility"] is None
