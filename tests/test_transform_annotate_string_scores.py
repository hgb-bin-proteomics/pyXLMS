#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


@pytest.mark.external
def test1():
    from pyXLMS.transform import get_string_network

    network = get_string_network(
        ["CDC42", "CDK1", "KIF23", "PLK1", "RAC2", "RACGAP1", "RHOA", "RHOB"], 9606
    )
    assert len(network) == 28


def test2():
    from pyXLMS.transform._annotate_string_scores import STRING_ORGANISMS

    assert STRING_ORGANISMS["Homo sapiens"] == 9606
    assert STRING_ORGANISMS["Mus musculus"] == 10090
    assert STRING_ORGANISMS["Arabidopsis thaliana"] == 3702
    assert STRING_ORGANISMS["Saccharomyces cerevisiae"] == 4932
    assert STRING_ORGANISMS["Drosophila melanogaster"] == 7227
    assert STRING_ORGANISMS["Danio rerio"] == 7955
    assert STRING_ORGANISMS["Caenorhabditis elegans"] == 6239
    assert STRING_ORGANISMS["Escherichia coli str. K-12 substr. MG1655"] == 511145
    assert STRING_ORGANISMS["Pseudomonas aeruginosa PAO1"] == 208964


@pytest.mark.external
def test3():
    from pyXLMS.transform import get_string_ids

    assert get_string_ids(["p53", "BRCA1", "cdk2", "Q99835"], organism=9606) == {
        "p53": "9606.ENSP00000269305",
        "BRCA1": "9606.ENSP00000418960",
        "cdk2": "9606.ENSP00000266970",
        "Q99835": "9606.ENSP00000249373",
    }


@pytest.mark.external
def test4():
    from pyXLMS.transform import get_string_ids

    assert get_string_ids(
        ["p53", "BRCA1", "cdk2", "Q99835"], organism="Homo sapiens"
    ) == {
        "p53": "9606.ENSP00000269305",
        "BRCA1": "9606.ENSP00000418960",
        "cdk2": "9606.ENSP00000266970",
        "Q99835": "9606.ENSP00000249373",
    }


def test5():
    from pyXLMS.transform import get_string_ids

    with pytest.raises(
        KeyError,
        match=r"Could not resolve organism s. pyogenes, please specify taxon identifier manually!",
    ):
        _ = get_string_ids(["p53", "BRCA1"], organism="s. pyogenes")


@pytest.mark.external
def test6():
    from pyXLMS.transform import get_string_network

    assert get_string_network(
        ["9606.ENSP00000269305", "9606.ENSP00000418960"], organism=9606
    ) == {
        "9606.ENSP00000269305_9606.ENSP00000418960": {
            "A": "9606.ENSP00000269305",
            "B": "9606.ENSP00000418960",
            "score": 0.999,
            "nscore": 0.0,
            "fscore": 0.0,
            "pscore": 0.0,
            "ascore": 0.067,
            "escore": 0.895,
            "dscore": 0.5,
            "tscore": 0.999,
        }
    }


@pytest.mark.external
def test7():
    from pyXLMS.transform import get_string_network

    assert get_string_network(
        ["9606.ENSP00000269305", "9606.ENSP00000418960"], organism="Homo sapiens"
    ) == {
        "9606.ENSP00000269305_9606.ENSP00000418960": {
            "A": "9606.ENSP00000269305",
            "B": "9606.ENSP00000418960",
            "score": 0.999,
            "nscore": 0.0,
            "fscore": 0.0,
            "pscore": 0.0,
            "ascore": 0.067,
            "escore": 0.895,
            "dscore": 0.5,
            "tscore": 0.999,
        }
    }


def test8():
    from pyXLMS.transform import get_string_network

    with pytest.raises(
        KeyError,
        match=r"Could not resolve organism s. pyogenes, please specify taxon identifier manually!",
    ):
        _ = get_string_network(["p53", "BRCA1"], organism="s. pyogenes")


@pytest.mark.external
def test9():
    from pyXLMS import parser
    from pyXLMS.transform import filter_crosslink_type
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    xls = annotate_string_scores(xls, organism="Homo sapiens")
    inter = filter_crosslink_type(xls)["Inter"]
    example = inter[4]  # example link with STRING score
    assert example["additional_information"][
        "pyXLMS_annotated_STRING_interactions"
    ] == [
        {
            "A": "9606.ENSP00000441875",
            "B": "9606.ENSP00000479488",
            "score": 0.999,
            "nscore": 0.0,
            "fscore": 0.0,
            "pscore": 0.068,
            "ascore": 0.923,
            "escore": 0.973,
            "dscore": 0.9,
            "tscore": 0.988,
        }
    ]
    assert example["additional_information"]["pyXLMS_annotated_STRING_score"] == 0.999


@pytest.mark.external
def test10():
    from pyXLMS import parser
    from pyXLMS.transform import filter_crosslink_type
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    xls = annotate_string_scores(xls, organism=9606)
    inter = filter_crosslink_type(xls)["Inter"]
    example = inter[4]  # example link with STRING score
    assert example["additional_information"][
        "pyXLMS_annotated_STRING_interactions"
    ] == [
        {
            "A": "9606.ENSP00000441875",
            "B": "9606.ENSP00000479488",
            "score": 0.999,
            "nscore": 0.0,
            "fscore": 0.0,
            "pscore": 0.068,
            "ascore": 0.923,
            "escore": 0.973,
            "dscore": 0.9,
            "tscore": 0.988,
        }
    ]
    assert example["additional_information"]["pyXLMS_annotated_STRING_score"] == 0.999


@pytest.mark.external
def test11():
    from pyXLMS import parser
    from pyXLMS.transform import filter_crosslink_type
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    xls = annotate_string_scores(xls, organism="Homo sapiens")
    inter = filter_crosslink_type(xls)["Inter"]
    for item in inter:
        assert "pyXLMS_annotated_STRING_interactions" in item["additional_information"]
        assert "pyXLMS_annotated_STRING_score" in item["additional_information"]


@pytest.mark.external
def test12():
    from pyXLMS import parser
    from pyXLMS.transform import filter_crosslink_type
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    for xl in xls:
        xl.additional_information = None
    xls = annotate_string_scores(xls, organism="Homo sapiens")
    inter = filter_crosslink_type(xls)["Inter"]
    for item in inter:
        assert "pyXLMS_annotated_STRING_interactions" in item["additional_information"]
        assert "pyXLMS_annotated_STRING_score" in item["additional_information"]


def test13():
    from pyXLMS.transform import annotate_string_scores

    with pytest.raises(
        KeyError,
        match=r"Could not resolve organism s. pyogenes, please specify taxon identifier manually!",
    ):
        _ = annotate_string_scores(["p53", "BRCA1"], organism="s. pyogenes")


@pytest.mark.external
def test14():
    from pyXLMS import parser
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    xls[0] = xls[0].copy_with_update({"alpha_proteins": None})
    with pytest.warns(
        RuntimeWarning,
        match=r"Some of your crosslink-spectrum-matches/crosslinks do not have associated proteins. Their STRING scores will be nan!",
    ):
        _ = annotate_string_scores(xls, organism="Homo sapiens")


def test15():
    from pyXLMS import parser
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    xls[0] = xls[0].copy_with_update({"alpha_proteins": None})
    with pytest.raises(
        RuntimeError,
        match=r"Some of your crosslink-spectrum-matches/crosslinks do not have associated proteins!",
    ):
        _ = annotate_string_scores(xls, organism="Homo sapiens", verbose=2)


def test16():
    from pyXLMS import parser
    from pyXLMS.transform import filter_crosslink_type
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    xls = pr["crosslinks"]
    intra = filter_crosslink_type(xls)["Intra"]
    with pytest.raises(
        ValueError,
        match=r"Can't annotate STRING scores for input data because it does not contain inter-links!",
    ):
        _ = annotate_string_scores(intra, organism="Homo sapiens")


@pytest.mark.external
def test17():
    from pyXLMS import parser
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_msannika(
        "data/_test/annotate_string_scores/Nucleus_Rep1_Crosslinks.txt",
        unsafe=True,
        verbose=0,
    )
    xls = pr["crosslinks"]
    with pytest.raises(
        RuntimeError,
        match=r"More than 2000 proteins/STRING IDs specified: 17295. Please reduce the number of proteins for a successful request!",
    ):
        _ = annotate_string_scores(xls, organism="Homo sapiens", verbose=2)


@pytest.mark.external
def test18():
    from pyXLMS import parser
    from pyXLMS.transform import filter_crosslink_type
    from pyXLMS.transform import annotate_string_scores

    pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_CSMs.parquet")
    csms = pr["crosslink-spectrum-matches"]
    csms = annotate_string_scores(csms, organism="Homo sapiens")
    inter = filter_crosslink_type(csms)["Inter"]
    for item in inter:
        assert "pyXLMS_annotated_STRING_interactions" in item["additional_information"]
        assert "pyXLMS_annotated_STRING_score" in item["additional_information"]
