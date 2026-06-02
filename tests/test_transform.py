#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS import transform

    modifications = {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)}
    modifications_str = "(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])"
    assert transform.modifications_to_str(modifications) == modifications_str


def test2():
    from pyXLMS import transform

    modifications = {1: ("Oxidation", 15.994915)}
    modifications_str = "(1:[Oxidation|15.994915])"
    assert transform.modifications_to_str(modifications) == modifications_str


def test3():
    from pyXLMS import transform

    modifications = dict()
    modifications_str = ""
    assert transform.modifications_to_str(modifications) == modifications_str


def test4():
    from pyXLMS import transform

    modifications = None
    assert transform.modifications_to_str(modifications) is None


def test5():
    from pyXLMS import transform
    from pyXLMS import data

    data_list = [
        data.create_crosslink_min("PEPK", 4, "PKEP", 2),
        data.create_crosslink_min("KPEP", 1, "PEKP", 3),
    ]
    assert transform.assert_data_type_same(data_list)

    data_list = [
        data.create_crosslink_min("PEPK", 4, "PKEP", 2),
        data.create_csm_min("KPEP", 1, "PEKP", 3, "RUN_1", 1),
    ]
    assert not transform.assert_data_type_same(data_list)


def test6():
    from pyXLMS import transform
    from pyXLMS import data

    data_list = [
        data.create_crosslink_min("PEPK", 4, "PKEP", 2),
        data.create_crosslink_min("KPEP", 1, "PEKP", 3),
    ]
    k = transform.get_available_keys(data_list)

    assert k["data_type"]
    assert k["completeness"]
    assert k["alpha_peptide"]
    assert k["alpha_peptide_crosslink_position"]
    assert not k["alpha_proteins"]
    assert not k["alpha_proteins_crosslink_positions"]
    assert not k["alpha_decoy"]
    assert k["beta_peptide"]
    assert k["beta_peptide_crosslink_position"]
    assert not k["beta_proteins"]
    assert not k["beta_proteins_crosslink_positions"]
    assert not k["beta_decoy"]
    assert k["crosslink_type"]
    assert not k["score"]
    assert not k["additional_information"]


def test7():
    from pyXLMS import transform
    from pyXLMS import data

    data_list = [
        data.create_crosslink_min("PEPK", 4, "PKEP", 2, score=170.3),
        data.create_crosslink_min("KPEP", 1, "PEKP", 3, score=13.5),
    ]
    k = transform.get_available_keys(data_list)

    assert k["data_type"]
    assert k["completeness"]
    assert k["alpha_peptide"]
    assert k["alpha_peptide_crosslink_position"]
    assert not k["alpha_proteins"]
    assert not k["alpha_proteins_crosslink_positions"]
    assert not k["alpha_decoy"]
    assert k["beta_peptide"]
    assert k["beta_peptide_crosslink_position"]
    assert not k["beta_proteins"]
    assert not k["beta_proteins_crosslink_positions"]
    assert not k["beta_decoy"]
    assert k["crosslink_type"]
    assert k["score"]
    assert not k["additional_information"]


def test8():
    from pyXLMS import transform
    from pyXLMS import data

    data_list = [
        data.create_crosslink(
            "PEPK",
            4,
            ["PROT1"],
            [4],
            False,
            "PKEP",
            2,
            ["PROT2"],
            [2],
            False,
            170.3,
            {},
        ),
        data.create_crosslink(
            "KPEP", 1, ["PROT3"], [1], False, "PEKP", 3, ["PROT4"], [3], True, 13.5, {}
        ),
    ]
    k = transform.get_available_keys(data_list)

    assert k["data_type"]
    assert k["completeness"]
    assert k["alpha_peptide"]
    assert k["alpha_peptide_crosslink_position"]
    assert k["alpha_proteins"]
    assert k["alpha_proteins_crosslink_positions"]
    assert k["alpha_decoy"]
    assert k["beta_peptide"]
    assert k["beta_peptide_crosslink_position"]
    assert k["beta_proteins"]
    assert k["beta_proteins_crosslink_positions"]
    assert k["beta_decoy"]
    assert k["crosslink_type"]
    assert k["score"]
    assert k["additional_information"]


def test9():
    from pyXLMS import transform
    from pyXLMS import data

    data_list = [
        data.create_csm_min("PEPK", 4, "PKEP", 2, "RUN_1", 1),
        data.create_csm_min("KPEP", 1, "PEKP", 3, "RUN_1", 2),
    ]
    k = transform.get_available_keys(data_list)

    assert k["data_type"]
    assert k["completeness"]
    assert k["alpha_peptide"]
    assert not k["alpha_modifications"]
    assert k["alpha_peptide_crosslink_position"]
    assert not k["alpha_proteins"]
    assert not k["alpha_proteins_crosslink_positions"]
    assert not k["alpha_proteins_peptide_positions"]
    assert not k["alpha_score"]
    assert not k["alpha_decoy"]
    assert k["beta_peptide"]
    assert not k["beta_modifications"]
    assert k["beta_peptide_crosslink_position"]
    assert not k["beta_proteins"]
    assert not k["beta_proteins_crosslink_positions"]
    assert not k["beta_proteins_peptide_positions"]
    assert not k["beta_score"]
    assert not k["beta_decoy"]
    assert k["crosslink_type"]
    assert not k["score"]
    assert k["spectrum_file"]
    assert k["scan_nr"]
    assert not k["charge"]
    assert not k["retention_time"]
    assert not k["ion_mobility"]
    assert not k["additional_information"]


def test10():
    from pyXLMS import transform
    from pyXLMS import data

    data_list = [
        data.create_csm_min("PEPK", 4, "PKEP", 2, "RUN_1", 1, score=170.3),
        data.create_csm_min("KPEP", 1, "PEKP", 3, "RUN_1", 2, score=13.5),
    ]
    k = transform.get_available_keys(data_list)

    assert k["data_type"]
    assert k["completeness"]
    assert k["alpha_peptide"]
    assert not k["alpha_modifications"]
    assert k["alpha_peptide_crosslink_position"]
    assert not k["alpha_proteins"]
    assert not k["alpha_proteins_crosslink_positions"]
    assert not k["alpha_proteins_peptide_positions"]
    assert not k["alpha_score"]
    assert not k["alpha_decoy"]
    assert k["beta_peptide"]
    assert not k["beta_modifications"]
    assert k["beta_peptide_crosslink_position"]
    assert not k["beta_proteins"]
    assert not k["beta_proteins_crosslink_positions"]
    assert not k["beta_proteins_peptide_positions"]
    assert not k["beta_score"]
    assert not k["beta_decoy"]
    assert k["crosslink_type"]
    assert k["score"]
    assert k["spectrum_file"]
    assert k["scan_nr"]
    assert not k["charge"]
    assert not k["retention_time"]
    assert not k["ion_mobility"]
    assert not k["additional_information"]


def test11():
    from pyXLMS import transform
    from pyXLMS import data

    data_list = [
        data.create_csm(
            "PEPK",
            {},
            4,
            ["PROT1"],
            [4],
            [1],
            170.3,
            False,
            "PKEP",
            {},
            2,
            ["PROT2"],
            [2],
            [1],
            150.4,
            False,
            150.4,
            "RUN_1",
            1,
            4,
            213.123,
            -50.0,
            {},
        ),
        data.create_csm(
            "KPEP",
            {},
            1,
            ["PROT3"],
            [1],
            [1],
            171.3,
            False,
            "PEKP",
            {},
            3,
            ["PROT4"],
            [3],
            [1],
            13.5,
            True,
            13.5,
            "RUN_1",
            2,
            3,
            213.124,
            -70.0,
            {},
        ),
    ]
    k = transform.get_available_keys(data_list)

    assert k["data_type"]
    assert k["completeness"]
    assert k["alpha_peptide"]
    assert k["alpha_modifications"]
    assert k["alpha_peptide_crosslink_position"]
    assert k["alpha_proteins"]
    assert k["alpha_proteins_crosslink_positions"]
    assert k["alpha_proteins_peptide_positions"]
    assert k["alpha_score"]
    assert k["alpha_decoy"]
    assert k["beta_peptide"]
    assert k["beta_modifications"]
    assert k["beta_peptide_crosslink_position"]
    assert k["beta_proteins"]
    assert k["beta_proteins_crosslink_positions"]
    assert k["beta_proteins_peptide_positions"]
    assert k["beta_score"]
    assert k["beta_decoy"]
    assert k["crosslink_type"]
    assert k["score"]
    assert k["spectrum_file"]
    assert k["scan_nr"]
    assert k["charge"]
    assert k["retention_time"]
    assert k["ion_mobility"]
    assert k["additional_information"]


def test12():
    from pyXLMS import parser
    from pyXLMS import transform

    pr = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    assert transform.display(pr, return_str=True) == (
        """Data Type:                            parser_result
Completeness:                         full
Identifying Search Engine:            MS Annika
Number of Crosslink-Spectrum-Matches: 826
Number of Crosslinks:                 300"""
    )

    csms = pr["crosslink-spectrum-matches"]
    assert transform.display(csms[0], return_str=True) == (
        """Data Type:                          crosslink-spectrum-match
Completeness:                       full
Alpha Peptide:                      GQKNSR
Alpha Modifications:                {3: ('DSS', 138.06808)}
Alpha Peptide Crosslink Position:   3
Alpha Proteins:                     ['Cas9']
Alpha Proteins Crosslink Positions: [779]
Alpha Proteins Peptide Positions:   [777]
Alpha Peptide Score:                119.82548987540834
Alpha Decoy:                        False
Beta Peptide:                       GQKNSR
Beta Modifications:                 {3: ('DSS', 138.06808)}
Beta Peptide Crosslink Position:    3
Beta Proteins:                      ['Cas9']
Beta Proteins Crosslink Positions:  [779]
Beta Proteins Peptide Positions:    [777]
Beta Peptide Score:                 119.82547820493929
Beta Decoy:                         False
Crosslink Type:                     intra
CSM Score:                          119.82547820493929
Spectrum File:                      XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw
Scan Number:                        2257
Precursor Charge:                   3
Retention Time:                     733.1895599999999
Ion Mobility/FAIMS CV:              0.0"""
    )

    xls = pr["crosslinks"]
    assert transform.display(xls[0], return_str=True) == (
        """Data Type:                          crosslink
Completeness:                       full
Alpha Peptide:                      GQKNSR
Alpha Peptide Crosslink Position:   3
Alpha Proteins:                     ['Cas9']
Alpha Proteins Crosslink Positions: [779]
Alpha Decoy:                        False
Beta Peptide:                       GQKNSR
Beta Peptide Crosslink Position:    3
Beta Proteins:                      ['Cas9']
Beta Proteins Crosslink Positions:  [779]
Beta Decoy:                         False
Crosslink Type:                     intra
Crosslink Score:                    119.82547820493929"""
    )


def test13():
    from pyXLMS import parser
    from pyXLMS import transform

    pr = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    assert isinstance(
        transform.display(pr, show_additional_information=True, return_str=True), str
    )

    csms = pr["crosslink-spectrum-matches"]
    assert "Additional Information" in transform.display(
        csms[0], show_additional_information=True, return_str=True
    )
    xls = pr["crosslinks"]
    assert "Additional Information" in transform.display(
        xls[0], show_additional_information=True, return_str=True
    )


def test14():
    from pyXLMS import parser
    from pyXLMS import transform

    pr = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    assert transform.display(pr) is None

    csms = pr["crosslink-spectrum-matches"]
    assert transform.display(csms[0]) is None
    xls = pr["crosslinks"]
    assert transform.display(xls[0]) is None


def test15():
    from pyXLMS.transform import check_available_keys
    from pyXLMS import data

    data_list = [
        data.create_crosslink_min("PEPK", 4, "PKEP", 2),
        data.create_crosslink_min("KPEP", 1, "PEKP", 3),
    ]
    assert check_available_keys(["alpha_peptide"], data_list)

    with pytest.raises(
        ValueError,
        match=r"Attribute 'score' is missing in at least one element but is required!",
    ):
        _ = check_available_keys(["score"], data_list)
