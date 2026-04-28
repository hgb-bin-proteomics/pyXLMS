#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from pyXLMS.transform import get_string_network

    network = get_string_network(
        ["CDC42", "CDK1", "KIF23", "PLK1", "RAC2", "RACGAP1", "RHOA", "RHOB"], 9606
    )
    assert len(network) == 28


def test2():
    from pyXLMS.transform.annotate_string_scores import STRING_ORGANISMS

    assert STRING_ORGANISMS["Homo sapiens"] == 9606
    assert STRING_ORGANISMS["Mus musculus"] == 10090
    assert STRING_ORGANISMS["Arabidopsis thaliana"] == 3702
    assert STRING_ORGANISMS["Saccharomyces cerevisiae"] == 4932
    assert STRING_ORGANISMS["Drosophila melanogaster"] == 7227
    assert STRING_ORGANISMS["Danio rerio"] == 7955
    assert STRING_ORGANISMS["Caenorhabditis elegans"] == 6239
    assert STRING_ORGANISMS["Escherichia coli str. K-12 substr. MG1655"] == 511145
    assert STRING_ORGANISMS["Pseudomonas aeruginosa PAO1"] == 208964
