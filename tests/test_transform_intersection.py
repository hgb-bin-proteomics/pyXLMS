#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import aggregate
    from pyXLMS.transform import intersection

    msannika = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    msannika = aggregate(msannika["crosslink-spectrum-matches"])
    assert len(msannika) == 235
    maxquant = pipeline(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    maxquant = aggregate(maxquant["crosslink-spectrum-matches"])
    assert len(maxquant) == 226
    crosslinks_intersection = intersection(msannika, maxquant)
    assert len(crosslinks_intersection) == 206


def test2():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import aggregate
    from pyXLMS.transform import intersection

    msannika = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    msannika = aggregate(msannika["crosslink-spectrum-matches"])
    assert len(msannika) == 235
    maxquant = pipeline(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    maxquant = aggregate(maxquant["crosslink-spectrum-matches"])
    assert len(maxquant) == 226
    crosslinks_intersection = intersection(msannika, maxquant)
    len(crosslinks_intersection) == 206
    plink = pipeline(
        "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv",
        engine="pLink",
        crosslinker="DSS",
    )
    plink = aggregate(plink["crosslink-spectrum-matches"])
    assert len(plink) == 252
    crosslinks_intersection = intersection(crosslinks_intersection, plink)
    assert len(crosslinks_intersection) == 203


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import intersection

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = (
        r"Parameter 'use' has to be one 'better_score', 'data_a', or 'data_b'! Option 'better_score' will return the intersection "
        r"with the better score \(if scores are available\)\. Option 'data_a' will return the intersection using elements of 'data_a'\. "
        r"Option 'data_b' will return the intersection using elements of 'data_b'\."
    )
    with pytest.raises(TypeError, match=err_str):
        _u = intersection(
            pr["crosslink-spectrum-matches"],
            pr["crosslink-spectrum-matches"],
            use="Set A",
        )


def test4():
    from pyXLMS.parser import read
    from pyXLMS.transform import intersection

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = (
        r"Parameter 'by' has to be one of 'peptide' or 'protein'! Option 'peptide' will group by peptide sequence and "
        r"peptide crosslink position while option 'protein' will group by protein identifier and protein crosslink position."
    )
    with pytest.raises(TypeError, match=err_str):
        _u = intersection(
            pr["crosslink-spectrum-matches"],
            pr["crosslink-spectrum-matches"],
            by="sequence",
        )


def test5():
    from pyXLMS.parser import read
    from pyXLMS.transform import intersection

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = (
        r"Parameter 'score' has to be one of 'higher_better' or 'lower_better'! If two identical crosslinks or crosslink-spectrum"
        r"-matches are found, the one with the higher score is kept if 'higher_better' is selected, and vice versa."
    )
    with pytest.raises(TypeError, match=err_str):
        _u = intersection(
            pr["crosslink-spectrum-matches"],
            pr["crosslink-spectrum-matches"],
            score="lower",
        )


def test6():
    from pyXLMS.parser import read
    from pyXLMS.transform import intersection

    pr = read(
        "data/_test/aggregate/csms_score.txt",
        engine="custom",
        crosslinker="DSS",
    )
    assert len(pr["crosslink-spectrum-matches"]) == 10
    err_str = r"Verbose level has to be one of 0, 1, or 2!"
    with pytest.raises(TypeError, match=err_str):
        _u = intersection(
            pr["crosslink-spectrum-matches"],
            pr["crosslink-spectrum-matches"],
            verbose=3,
        )


def test7():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import aggregate
    from pyXLMS.transform import intersection

    msannika = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    msannika = aggregate(msannika["crosslink-spectrum-matches"])
    maxquant = pipeline(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    maxquant = maxquant["crosslink-spectrum-matches"]
    err_str = r"Parameters 'data_a' and 'data_b' have to be of the same data type!"
    with pytest.raises(TypeError, match=err_str):
        _u = intersection(msannika, maxquant)


def test8():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import aggregate
    from pyXLMS.transform import intersection

    msannika = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    msannika = aggregate(msannika["crosslink-spectrum-matches"])
    maxquant = pipeline(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    maxquant = aggregate(maxquant["crosslink-spectrum-matches"])
    maxquant[-5]["completeness"] = "partial"
    maxquant[-5]["score"] = None
    err_str = (
        r"Can't intersect based on score because not all data have associated scores!"
    )
    with pytest.raises(ValueError, match=err_str):
        _u = intersection(msannika, maxquant)


def test9():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import intersection

    msannika = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    msannika = msannika["crosslink-spectrum-matches"]
    maxquant = pipeline(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    maxquant = maxquant["crosslink-spectrum-matches"]
    u = intersection(msannika, maxquant, verbose=0)
    assert u is not None


def test10():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import intersection

    msannika = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    msannika = msannika["crosslink-spectrum-matches"]
    maxquant = pipeline(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    maxquant = maxquant["crosslink-spectrum-matches"]
    warn_str = r"Creating intersection of crosslink-spectrum-matches\. Be sure that this makes sense for your data!"
    with pytest.warns(RuntimeWarning, match=warn_str):
        _u = intersection(msannika, maxquant)


def test11():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import intersection

    msannika = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    msannika = msannika["crosslink-spectrum-matches"]
    maxquant = pipeline(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    maxquant = maxquant["crosslink-spectrum-matches"]
    err_str = (
        r"Can't create intersection of crosslink-spectrum-matches for verbose level 2!"
    )
    with pytest.raises(RuntimeError, match=err_str):
        _u = intersection(msannika, maxquant, verbose=2)
