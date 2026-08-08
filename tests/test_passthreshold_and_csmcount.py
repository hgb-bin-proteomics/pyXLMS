#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def test1():
    # Every Scout crosslink carries its integer CSM count. pyXLMS >=2.0 records are
    # frozen and forbid extra attributes, so this optional per-record metric is
    # attached through the (mutable) additional_information field rather than as a
    # new top-level key.
    from pyXLMS.parser import read_scout

    crosslinks = read_scout("data/scout/Cas9_Residue_Pairs.csv", crosslinker="DSSO")[
        "crosslinks"
    ]
    assert len(crosslinks) > 0
    for xl in crosslinks:
        assert isinstance(xl.additional_information["csm_count"], int)
        assert xl.additional_information["csm_count"] >= 1
    assert crosslinks[0].additional_information["csm_count"] == 26


def test2():
    # mzIdentML passThreshold is carried through for every CSM via
    # additional_information (same frozen-record reason as above).
    from pyXLMS.parser import read_mzid

    # this Scout mzid has spectrum IDs without "scan=", so read_mzid emits the same
    # RuntimeWarning the other mzid tests expect.
    with pytest.warns(RuntimeWarning):
        result = read_mzid(
            "data/scout2/Cas9_HeLa_Cyt_r2/Cas9_HeLa_Cyt_r2-v1.2.mzid", verbose=0
        )
    csms = result["crosslink-spectrum-matches"]
    assert len(csms) > 0
    for csm in csms:
        assert csm.additional_information["pass_threshold"] is True
