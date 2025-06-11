#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest
from typing import List, Dict, Any


def get_fdr_strict(data: List[Dict[str, Any]]) -> float:
    D = 0
    T = 0
    for item in data:
        if not item["alpha_decoy"] and not item["beta_decoy"]:
            T += 1
        else:
            D += 1
    return D / T


def get_fdr_relaxed(data: List[Dict[str, Any]]) -> float:
    D = 0
    DT = 0
    T = 0
    for item in data:
        if not item["alpha_decoy"] and not item["beta_decoy"]:
            T += 1
        elif item["alpha_decoy"] and item["beta_decoy"]:
            D += 1
        else:
            DT += 1
    if (DT - D) < 0.0:
        raise RuntimeError("Negative FDR!")
    return (DT - D) / T


def test1():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 826
    validated = validate(csms)
    assert len(validated) == 705


def test2():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        [
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        ],
        engine="MS Annika",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 826
    assert len(pr["crosslinks"]) == 300
    validated = validate(pr)
    assert len(validated["crosslink-spectrum-matches"]) == 705
    assert len(validated["crosslinks"]) == 226


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        [
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        ],
        engine="MS Annika",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 826
    assert len(pr["crosslinks"]) == 300
    validated = validate(pr, fdr=0.05)
    assert len(validated["crosslink-spectrum-matches"]) == 825
    assert len(validated["crosslinks"]) == 260


def test4():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    with pytest.raises(
        ValueError,
        match=r"FDR must be given as a real number between 0 and 1, e\.g\. 0\.01 corresponds to 1\% FDR!",
    ):
        _validated = validate(pr, fdr=1.0)


def test5():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    err_str = (
        r"Parameter 'formula' has to be one of 'D\/T', '\(TD\+DD\)\/TT' or '\(TD\-DD\)\/TT'! Where D and DD is the number of decoys, T and TT the number of targets, "
        r"and TD the number of target-decoys!"
    )
    with pytest.raises(
        TypeError,
        match=err_str,
    ):
        _validated = validate(pr, formula="T/D")


def test6():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    err_str = (
        r"Parameter 'score' has to be one of 'higher_better' or 'lower_better'! If two identical crosslinks or crosslink-spectrum"
        r"-matches are found, the one with the higher score is kept if 'higher_better' is selected, and vice versa."
    )
    with pytest.raises(
        TypeError,
        match=err_str,
    ):
        _validated = validate(pr, score="lower")


def test7():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        "data/pyxlms/csm_min.txt",
        engine="custom",
        crosslinker="DSS",
    )

    err_str = (
        r"Can't validate data if 'score' or target\/decoy labels are missing! Selecting 'ignore_missing_labels \= True' will ignore crosslinks and crosslink-spectrum-matches "
        r"that don't have a valid target\/decoy label and filter them out!"
    )
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _validated = validate(pr)


def test8():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        "data/_test/validate/csms.txt",
        engine="custom",
        crosslinker="DSS",
    )

    err_str = r"Can't estimate FDR with formula '\(TD\-DD\)\/TT' when there are not TD matches! Please select the default formula instead!"
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _validated = validate(pr)
