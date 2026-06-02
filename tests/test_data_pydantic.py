#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from pyXLMS.data import CrosslinkSpectrumMatch as CSM

    csm = CSM(
        alpha_peptide="PEKP",
        alpha_peptide_crosslink_position=3,
        beta_peptide="TKIDE",
        beta_peptide_crosslink_position=2,
        spectrum_file="dsso.mzML",
        scan_nr=1,
    )
    assert csm is not None


def test2():
    from pyXLMS.data import CrosslinkSpectrumMatch as CSM

    csm = CSM(
        alpha_peptide="PEKP",
        alpha_peptide_crosslink_position=3,
        beta_peptide="TKIDE",
        beta_peptide_crosslink_position=2,
        spectrum_file="dsso.mzML",
        scan_nr=1,
    )
    csm_copy = csm.copy_with_update(update={"scan_nr": 2})
    assert csm_copy.scan_nr == 2


def test3():
    from pyXLMS import parser

    pr = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    assert csms[0].display() is None


def test4():
    from pyXLMS.data import create_csm_min

    csm = create_csm_min("PEPKTIDE", 4, "KPEPTIDE", 1, "RUN_1", 1)
    assert csm.to_proforma() == "KPEPTIDE//PEPKTIDE"


def test5():
    from pyXLMS.data import create_csm_min

    csm = create_csm_min("PEPKTIDE", 4, "KPEPTIDE", 1, "RUN_1", 1)
    assert (
        csm.to_proforma(crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PEPTIDE//PEPK[Xlink:DSSO]TIDE"
    )


def test6():
    from pyXLMS.data import create_csm_min

    csm = create_csm_min(
        "PEPKTIDE",
        4,
        "KPMEPTIDE",
        1,
        "RUN_1",
        1,
        modifications_b={3: ("Oxidation", 15.994915)},
    )
    assert (
        csm.to_proforma(crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE"
    )


def test7():
    from pyXLMS.data import create_csm_min

    csm = create_csm_min(
        "PEPKTIDE",
        4,
        "KPMEPTIDE",
        1,
        "RUN_1",
        1,
        modifications_b={3: ("Oxidation", 15.994915)},
        charge=3,
    )
    assert (
        csm.to_proforma(crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE/3"
    )


def test8():
    from pyXLMS.data import create_csm_min

    csm = create_csm_min(
        "PEPKTIDE",
        4,
        "KPMEPTIDE",
        1,
        "RUN_1",
        1,
        modifications_a={4: ("DSSO", 158.00376)},
        modifications_b={1: ("DSSO", 158.00376), 3: ("Oxidation", 15.994915)},
        charge=3,
    )
    assert (
        csm.to_proforma() == "K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3"
    )


def test9():
    from pyXLMS.data import create_csm_min

    csm = create_csm_min(
        "PEPKTIDE",
        4,
        "KPMEPTIDE",
        1,
        "RUN_1",
        1,
        modifications_a={4: ("DSSO", 158.00376)},
        modifications_b={1: ("DSSO", 158.00376), 3: ("Oxidation", 15.994915)},
        charge=3,
    )
    assert (
        csm.to_proforma(crosslinker="Xlink:DSSO")
        == "K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3"
    )


def test10():
    from pyXLMS.data import Crosslink

    xl = Crosslink(
        alpha_peptide="PEKP",
        alpha_peptide_crosslink_position=3,
        beta_peptide="TKIDE",
        beta_peptide_crosslink_position=2,
    )
    assert xl is not None


def test11():
    from pyXLMS.data import Crosslink

    xl = Crosslink(
        alpha_peptide="PEKP",
        alpha_peptide_crosslink_position=3,
        alpha_proteins=["PROT"],
        beta_peptide="PEKP",
        beta_peptide_crosslink_position=3,
        beta_proteins=["PROT"],
    )
    xl_copy = xl.copy_with_update(
        update={"additional_information": {"homomeric": True}}
    )
    assert xl_copy["additional_information"]["homomeric"]


def test12():
    from pyXLMS import parser

    pr = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    xls = pr["crosslinks"]
    assert xls[0].display() is None


def test13():
    from pyXLMS.data import create_crosslink_min

    xl = create_crosslink_min("PEPKTIDE", 4, "KPEPTIDE", 1)
    assert xl.to_proforma() == "KPEPTIDE//PEPKTIDE"


def test14():
    from pyXLMS.data import create_crosslink_min

    xl = create_crosslink_min("PEPKTIDE", 4, "KPEPTIDE", 1)
    assert (
        xl.to_proforma(crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PEPTIDE//PEPK[Xlink:DSSO]TIDE"
    )


def test15():
    from pyXLMS.data import Crosslink
    from pyXLMS.data import ParserResult

    xl = Crosslink(
        alpha_peptide="PEKP",
        alpha_peptide_crosslink_position=3,
        beta_peptide="TKIDE",
        beta_peptide_crosslink_position=2,
    )
    pr = ParserResult(search_engine="My Search Engine", crosslinks=[xl])
    assert pr is not None


def test16():
    from pyXLMS.data import Crosslink
    from pyXLMS.data import ParserResult

    pr = ParserResult(search_engine="My Search Engine")
    xl = Crosslink(
        alpha_peptide="PEKP",
        alpha_peptide_crosslink_position=3,
        beta_peptide="TKIDE",
        beta_peptide_crosslink_position=2,
    )
    pr_copy = pr.copy_with_update(update={"crosslinks": [xl]})
    assert len(pr_copy.xls()) == 1


def test17():
    from pyXLMS import parser

    pr = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    assert pr.display() is None
