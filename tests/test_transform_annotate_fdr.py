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
    from pyXLMS.transform import annotate_fdr

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    csms = annotate_fdr(csms)
    validated_csms = [
        csm
        for csm in csms
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert len(validated_csms) == 705


def test2():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

    pr = read(
        [
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
            "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
        ],
        engine="MS Annika",
        crosslinker="DSS",
    )
    pr = annotate_fdr(pr)
    validated_csms = [
        csm
        for csm in pr["crosslink-spectrum-matches"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert len(validated_csms) == 705
    validated_xls = [
        xl
        for xl in pr["crosslinks"]
        if xl["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert len(validated_xls) == 226


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 826
    csms = annotate_fdr(csms)
    assert len(csms) == 826
    for csm in csms:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    validated = [
        csm
        for csm in csms
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert len(validated) == 705


def test4():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
    pr = annotate_fdr(pr)
    assert len(pr["crosslink-spectrum-matches"]) == 826
    assert len(pr["crosslinks"]) == 300
    for csm in pr["crosslink-spectrum-matches"]:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    for xl in pr["crosslinks"]:
        assert "pyXLMS_annotated_FDR" in xl["additional_information"]
    validated_csms = [
        csm
        for csm in pr["crosslink-spectrum-matches"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    validated_xls = [
        xl
        for xl in pr["crosslinks"]
        if xl["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert len(validated_csms) == 705
    assert len(validated_xls) == 226


def test5():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
    pr = annotate_fdr(pr)
    assert len(pr["crosslink-spectrum-matches"]) == 826
    assert len(pr["crosslinks"]) == 300
    for csm in pr["crosslink-spectrum-matches"]:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    for xl in pr["crosslinks"]:
        assert "pyXLMS_annotated_FDR" in xl["additional_information"]
    validated_csms = [
        csm
        for csm in pr["crosslink-spectrum-matches"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.05
    ]
    validated_xls = [
        xl
        for xl in pr["crosslinks"]
        if xl["additional_information"]["pyXLMS_annotated_FDR"] <= 0.05
    ]
    assert len(validated_csms) == 825
    assert len(validated_xls) == 260


def test6():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
        _annotated = annotate_fdr(pr, formula="T/D")


def test7():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
    )

    err_str = r"Parameter 'score' has to be one of 'higher_better' or 'lower_better'!"
    with pytest.raises(
        TypeError,
        match=err_str,
    ):
        _annotated = annotate_fdr(pr, score="lower")


def test8():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

    pr = read(
        "data/pyxlms/csm_min.txt",
        engine="custom",
        crosslinker="DSS",
    )

    err_str = (
        r"Can't annotate data if 'score' or target\/decoy labels are missing! Selecting 'ignore_missing_labels \= True' will ignore crosslinks and crosslink-spectrum-matches "
        r"that don't have a valid target\/decoy label and filter them out!"
    )
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _annotated = annotate_fdr(pr)


def test9():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

    pr = read(
        "data/_test/validate/csms.txt",
        engine="custom",
        crosslinker="DSS",
    )

    err_str = r"Can't annotate FDR with formula '\(TD\-DD\)\/TT' when there are no TD matches! Please select the default formula instead!"
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _annotated = annotate_fdr(pr, formula="(TD-DD)/TT")


def test10():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
    pr = annotate_fdr(pr, formula="(TD+DD)/TT")
    assert len(pr["crosslink-spectrum-matches"]) == 826
    assert len(pr["crosslinks"]) == 300
    for csm in pr["crosslink-spectrum-matches"]:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    for xl in pr["crosslinks"]:
        assert "pyXLMS_annotated_FDR" in xl["additional_information"]
    validated_csms = [
        csm
        for csm in pr["crosslink-spectrum-matches"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    validated_xls = [
        xl
        for xl in pr["crosslinks"]
        if xl["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert get_fdr_strict(validated_csms) < 0.01
    assert get_fdr_strict(validated_xls) < 0.01


def test11():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
    pr["crosslinks"] = None
    pr = annotate_fdr(pr, formula="(TD-DD)/TT")
    assert len(pr["crosslink-spectrum-matches"]) == 826
    assert pr["crosslinks"] is None
    for csm in pr["crosslink-spectrum-matches"]:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    validated_csms = [
        csm
        for csm in pr["crosslink-spectrum-matches"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert get_fdr_relaxed(validated_csms) < 0.01


def test12():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
    csms = annotate_fdr(
        pr["crosslink-spectrum-matches"],
        formula="(TD+DD)/TT",
        score="lower_better",
    )
    assert len(csms) == 826
    for csm in csms:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    validated_csms = [
        csm
        for csm in csms
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert get_fdr_strict(validated_csms) < 0.01


def test13():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
    csms = annotate_fdr(
        pr["crosslink-spectrum-matches"],
        formula="(TD-DD)/TT",
        score="lower_better",
    )
    assert len(csms) == 826
    for csm in csms:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    validated_csms = [
        csm
        for csm in csms
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert get_fdr_relaxed(validated_csms) < 0.01


@pytest.mark.slow
def test14():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr
    from pyXLMS.transform import filter_crosslink_type

    pr = read(
        "data/_test/validate/csms_25000.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
        unsafe=True,
        verbose=0,
    )
    for csm in pr["crosslink-spectrum-matches"]:
        csm["score"] = -csm["score"]
    csms = annotate_fdr(
        pr["crosslink-spectrum-matches"],
        formula="(TD+DD)/TT",
        score="lower_better",
        separate_intra_inter=True,
    )
    for csm in csms:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    separate = filter_crosslink_type(csms)
    intra = [
        csm
        for csm in separate["Intra"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    inter = [
        csm
        for csm in separate["Inter"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert get_fdr_strict(intra) < 0.01
    assert get_fdr_strict(inter) < 0.01


@pytest.mark.slow
def test15():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr
    from pyXLMS.transform import filter_crosslink_type

    pr = read(
        "data/_test/validate/csms_25000.xlsx",
        engine="MS Annika",
        crosslinker="DSS",
        unsafe=True,
        verbose=0,
    )
    for csm in pr["crosslink-spectrum-matches"]:
        csm["score"] = -csm["score"]
    csms = annotate_fdr(
        pr["crosslink-spectrum-matches"],
        formula="(TD-DD)/TT",
        score="lower_better",
        separate_intra_inter=True,
    )
    for csm in csms:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    separate = filter_crosslink_type(csms)
    intra = [
        csm
        for csm in separate["Intra"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    inter = [
        csm
        for csm in separate["Inter"]
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert get_fdr_relaxed(intra) < 0.01
    assert get_fdr_relaxed(inter) < 0.01


def test16():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
    csms = annotate_fdr(
        pr["crosslink-spectrum-matches"],
        formula="(TD+DD)/TT",
        score="lower_better",
        ignore_missing_labels=True,
    )
    for csm in csms:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    validated_csms = [
        csm
        for csm in csms
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert get_fdr_strict(validated_csms) < 0.01


def test17():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
    csms = annotate_fdr(
        pr["crosslink-spectrum-matches"],
        formula="(TD-DD)/TT",
        score="lower_better",
        ignore_missing_labels=True,
    )
    for csm in csms:
        assert "pyXLMS_annotated_FDR" in csm["additional_information"]
    validated_csms = [
        csm
        for csm in csms
        if csm["additional_information"]["pyXLMS_annotated_FDR"] <= 0.01
    ]
    assert get_fdr_relaxed(validated_csms) < 0.01


def test18():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
        r"Can't annotate data if 'score' or target\/decoy labels are missing! Selecting 'ignore_missing_labels \= True' will ignore crosslinks and crosslink-spectrum-matches "
        r"that don't have a valid target\/decoy label and filter them out!"
    )
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _annotated = annotate_fdr(
            pr["crosslink-spectrum-matches"],
            formula="(TD+DD)/TT",
            score="lower_better",
            ignore_missing_labels=False,
        )


def test19():
    from pyXLMS.parser import read
    from pyXLMS.transform import annotate_fdr

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
        r"Can't annotate data if 'score' or target\/decoy labels are missing! Selecting 'ignore_missing_labels \= True' will ignore crosslinks and crosslink-spectrum-matches "
        r"that don't have a valid target\/decoy label and filter them out!"
    )
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _annotated = annotate_fdr(
            pr["crosslink-spectrum-matches"],
            formula="(TD-DD)/TT",
            score="lower_better",
            ignore_missing_labels=False,
        )
