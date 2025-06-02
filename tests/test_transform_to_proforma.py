#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.transform import to_proforma

    psm = {"peptide": "PEPKTIDE"}

    with pytest.raises(
        TypeError,
        match=r"Unsupported data type for input data! Parameter data has to be a \(list of\) crosslink or crosslink-spectrum-match!",
    ):
        _r = to_proforma(psm)


def test2():
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.transform import to_proforma

    xl = create_crosslink_min("PEPKTIDE", 4, "KPEPTIDE", 1)
    assert to_proforma(xl) == "KPEPTIDE//PEPKTIDE"


def test3():
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.transform import to_proforma

    xl = create_crosslink_min("PEPKTIDE", 4, "KPEPTIDE", 1)
    assert (
        to_proforma(xl, crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PEPTIDE//PEPK[Xlink:DSSO]TIDE"
    )


def test4():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

    csm = create_csm_min("PEPKTIDE", 4, "KPEPTIDE", 1, "RUN_1", 1)
    assert to_proforma(csm) == "KPEPTIDE//PEPKTIDE"


def test5():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

    csm = create_csm_min("PEPKTIDE", 4, "KPEPTIDE", 1, "RUN_1", 1)
    assert (
        to_proforma(csm, crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PEPTIDE//PEPK[Xlink:DSSO]TIDE"
    )


def test6():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

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
        to_proforma(csm, crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE"
    )


def test7():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

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
        to_proforma(csm, crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE/3"
    )


def test8():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

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
        to_proforma(csm) == "K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3"
    )


def test9():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

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
        to_proforma(csm, crosslinker="Xlink:DSSO")
        == "K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3"
    )


def test10():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

    csm = create_csm_min(
        "PEPKTIDE",
        4,
        "MKPMEPTIDE",
        2,
        "RUN_1",
        1,
        modifications_a={4: ("DSSO", 158.00376)},
        modifications_b={
            2: ("DSSO", 158.00376),
            4: ("Oxidation", 15.994915),
            1: ("Oxidation", 15.994915),
        },
        charge=3,
    )
    assert (
        to_proforma(csm)
        == "M[+15.994915]K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3"
    )


def test11():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

    csm_1 = create_csm_min(
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
    csm_2 = create_csm_min(
        "PEPKTIDE",
        4,
        "MKPMEPTIDE",
        2,
        "RUN_1",
        2,
        modifications_a={4: ("DSSO", 158.00376)},
        modifications_b={
            2: ("DSSO", 158.00376),
            4: ("Oxidation", 15.994915),
            1: ("Oxidation", 15.994915),
        },
        charge=3,
    )
    proforma = to_proforma([csm_1, csm_2])
    assert proforma[0] == "K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3"
    assert (
        proforma[1]
        == "M[+15.994915]K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3"
    )


def test12():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma
    from pyXLMS.transform import modifications_to_str as mts

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
        to_proforma(csm, crosslinker="Xlink:DSSO")
        == "K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE/3"
    )
    assert mts(csm["alpha_modifications"]) == "(3:[Oxidation|15.994915])"


def test13():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

    csm = create_csm_min(
        "PEPKTIDE",
        4,
        "MKPMEPTIDE",
        1,
        "RUN_1",
        1,
        modifications_a={4: ("DSSO", 158.00376)},
        modifications_b={
            0: ("DSSO", 158.00376),
            4: ("Oxidation", 15.994915),
            1: ("Oxidation", 15.994915),
        },
        charge=3,
    )
    assert (
        to_proforma(csm)
        == "[+158.00376]-M[+15.994915]KPM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3"
    )


def test14():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

    csm = create_csm_min(
        "PEPKTIDE",
        4,
        "MKPMEPTIDE",
        10,
        "RUN_1",
        1,
        modifications_a={4: ("DSSO", 158.00376)},
        modifications_b={
            11: ("DSSO", 158.00376),
            4: ("Oxidation", 15.994915),
            1: ("Oxidation", 15.994915),
        },
        charge=3,
    )
    assert (
        to_proforma(csm)
        == "M[+15.994915]KPM[+15.994915]EPTIDE-[+158.00376]//PEPK[+158.00376]TIDE/3"
    )


def test15():
    from pyXLMS.data import create_csm_min
    from pyXLMS.transform import to_proforma

    csm = create_csm_min(
        "PEPKTIDE",
        4,
        "MKPMEPTIDE",
        10,
        "RUN_1",
        1,
        modifications_a={4: ("DSSO", 158.00376)},
        modifications_b={
            11: ("DSSO", 158.00376),
            4: ("Oxidation", 15.994915),
            1: ("Oxidation", float("nan")),
        },
        charge=3,
    )
    assert (
        to_proforma(csm)
        == "MKPM[+15.994915]EPTIDE-[+158.00376]//PEPK[+158.00376]TIDE/3"
    )
