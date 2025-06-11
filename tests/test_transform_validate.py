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
        _validated = validate(pr, formula="(TD-DD)/TT")


def test9():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate

    pr = read(
        "data/_test/validate/csms.txt",
        engine="custom",
        crosslinker="DSS",
    )

    err_str = r"None of the data passes the desired FDR threshold! This is usually due to decoys with very good scores."
    with pytest.warns(
        RuntimeWarning,
        match=err_str,
    ):
        validated = validate(pr)
        assert validated["crosslink-spectrum-matches"] == []


def test10():
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
    validated = validate(
        pr["crosslink-spectrum-matches"], fdr=0.01, formula="(TD+DD)/TT"
    )
    assert get_fdr_strict(validated) < 0.01


def test11():
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
    validated = validate(
        pr["crosslink-spectrum-matches"], fdr=0.01, formula="(TD-DD)/TT"
    )
    assert get_fdr_relaxed(validated) < 0.01


def test12():
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
    for csm in pr["crosslink-spectrum-matches"]:
        csm["score"] = -csm["score"]
    validated = validate(
        pr["crosslink-spectrum-matches"],
        fdr=0.01,
        formula="(TD+DD)/TT",
        score="lower_better",
    )
    assert get_fdr_strict(validated) < 0.01


def test13():
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
    for csm in pr["crosslink-spectrum-matches"]:
        csm["score"] = -csm["score"]
    validated = validate(
        pr["crosslink-spectrum-matches"],
        fdr=0.01,
        formula="(TD-DD)/TT",
        score="lower_better",
    )
    assert get_fdr_relaxed(validated) < 0.01


@pytest.mark.slow
def test14():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate
    from pyXLMS.transform import filter_crosslink_type

    pr = read(
        "data/_test/validate/csms.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
        unsafe=True,
        verbose=0,
    )
    for csm in pr["crosslink-spectrum-matches"]:
        csm["score"] = -csm["score"]
    validated = validate(
        pr["crosslink-spectrum-matches"],
        fdr=0.01,
        formula="(TD+DD)/TT",
        score="lower_better",
        separate_intra_inter=True,
    )
    separate = filter_crosslink_type(validated)
    assert get_fdr_strict(separate["Intra"]) < 0.01
    assert get_fdr_strict(separate["Inter"]) < 0.01


@pytest.mark.slow
def test15():
    from pyXLMS.parser import read
    from pyXLMS.transform import validate
    from pyXLMS.transform import filter_crosslink_type

    pr = read(
        "data/_test/validate/csms.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
        unsafe=True,
        verbose=0,
    )
    for csm in pr["crosslink-spectrum-matches"]:
        csm["score"] = -csm["score"]
    validated = validate(
        pr["crosslink-spectrum-matches"],
        fdr=0.01,
        formula="(TD-DD)/TT",
        score="lower_better",
        separate_intra_inter=True,
    )
    separate = filter_crosslink_type(validated)
    assert get_fdr_relaxed(separate["Intra"]) < 0.01
    assert get_fdr_relaxed(separate["Inter"]) < 0.01


def test16():
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
    for i, csm in enumerate(pr["crosslink-spectrum-matches"]):
        csm["score"] = -csm["score"]
        if i < 5:
            csm["alpha_decoy"] = None
        elif i < 10:
            csm["beta_decoy"] = None
    validated = validate(
        pr["crosslink-spectrum-matches"],
        fdr=0.01,
        formula="(TD+DD)/TT",
        score="lower_better",
        ignore_missing_labels=True,
    )
    assert get_fdr_strict(validated) < 0.01


def test17():
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
    for i, csm in enumerate(pr["crosslink-spectrum-matches"]):
        csm["score"] = -csm["score"]
        if i < 5:
            csm["alpha_decoy"] = None
        elif i < 10:
            csm["beta_decoy"] = None
    validated = validate(
        pr["crosslink-spectrum-matches"],
        fdr=0.01,
        formula="(TD-DD)/TT",
        score="lower_better",
        ignore_missing_labels=True,
    )
    assert get_fdr_relaxed(validated) < 0.01


def test18():
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
    for i, csm in enumerate(pr["crosslink-spectrum-matches"]):
        csm["score"] = -csm["score"]
        if i < 5:
            csm["completeness"] = "partial"
            csm["alpha_decoy"] = None
        elif i < 10:
            csm["completeness"] = "partial"
            csm["beta_decoy"] = None
    err_str = (
        r"Can't validate data if 'score' or target\/decoy labels are missing! Selecting 'ignore_missing_labels \= True' will ignore crosslinks and crosslink-spectrum-matches "
        r"that don't have a valid target\/decoy label and filter them out!"
    )
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _validated = validate(
            pr["crosslink-spectrum-matches"],
            fdr=0.01,
            formula="(TD+DD)/TT",
            score="lower_better",
            ignore_missing_labels=False,
        )


def test19():
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
    for i, csm in enumerate(pr["crosslink-spectrum-matches"]):
        csm["score"] = -csm["score"]
        if i < 5:
            csm["completeness"] = "partial"
            csm["alpha_decoy"] = None
        elif i < 10:
            csm["completeness"] = "partial"
            csm["beta_decoy"] = None
    err_str = (
        r"Can't validate data if 'score' or target\/decoy labels are missing! Selecting 'ignore_missing_labels \= True' will ignore crosslinks and crosslink-spectrum-matches "
        r"that don't have a valid target\/decoy label and filter them out!"
    )
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _validated = validate(
            pr["crosslink-spectrum-matches"],
            fdr=0.01,
            formula="(TD-DD)/TT",
            score="lower_better",
            ignore_missing_labels=False,
        )
