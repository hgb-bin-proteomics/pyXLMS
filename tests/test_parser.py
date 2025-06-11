#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

FILES = {
    "data/maxquant/run1/crosslinkMsms.txt": {
        "engine": "MaxQuant",
        "crosslinker": "DSS",
    },
    "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.mzid": {
        "engine": "mzIdentML",
        "crosslinker": "DSS",
    },
    "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult": {
        "engine": "MS Annika",
        "crosslinker": "DSS",
    },
    "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt": {
        "engine": "MS Annika",
        "crosslinker": "DSS",
    },
    "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx": {
        "engine": "MS Annika",
        "crosslinker": "DSS",
    },
    "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt": {
        "engine": "MS Annika",
        "crosslinker": "DSS",
    },
    "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx": {
        "engine": "MS Annika",
        "crosslinker": "DSS",
    },
    "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv": {
        "engine": "pLink",
        "crosslinker": "DSS",
    },
    "data/plink3/Cas10_plus10_2025.04.07.filtered_cross-linked_spectra.csv": {
        "engine": "pLink",
        "crosslinker": "DSS",
    },
    "data/pyxlms/csm.txt": {"engine": "Custom", "crosslinker": "DSS"},
    "data/pyxlms/xl.txt": {"engine": "Custom", "crosslinker": "DSS"},
    "data/scout/Cas9_Filtered_CSMs.csv": {"engine": "Scout", "crosslinker": "DSSO"},
    "data/scout/Cas9_Residue_Pairs.csv": {"engine": "Scout", "crosslinker": "DSSO"},
    "data/scout/Cas9_Unfiltered_CSMs.csv": {"engine": "Scout", "crosslinker": "DSSO"},
    "data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv": {
        "engine": "xiSearch/xiFDR",
        "crosslinker": "BS3",
    },
    "data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv": {
        "engine": "xiSearch/xiFDR",
        "crosslinker": "BS3",
    },
    "data/xi/r1_Xi1.7.6.7.csv": {"engine": "xiSearch/xiFDR", "crosslinker": "BS3"},
    "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3.mzid": {
        "engine": "mzIdentML",
        "crosslinker": "DSSO",
    },
    "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3.pdResult": {
        "engine": "XlinkX",
        "crosslinker": "DSSO",
    },
    "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_Crosslinks.txt": {
        "engine": "XlinkX",
        "crosslinker": "DSSO",
    },
    "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_Crosslinks.xlsx": {
        "engine": "XlinkX",
        "crosslinker": "DSSO",
    },
    "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_CSMs.txt": {
        "engine": "XlinkX",
        "crosslinker": "DSSO",
    },
    "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3_CSMs.xlsx": {
        "engine": "XlinkX",
        "crosslinker": "DSSO",
    },
}


def test1():
    from pyXLMS import parser as p

    with pytest.raises(ValueError):
        _r = p.read(
            "data/maxquant/run1/crosslinkMsms.txt",
            engine="MS Amanda",
            crosslinker="DSS",
        )


@pytest.mark.slow
def test2():
    from pyXLMS import parser as p

    for k, v in FILES.items():
        pr = p.read(k, engine=v["engine"], crosslinker=v["crosslinker"], verbose=0)
        assert pr["search_engine"] == v["engine"]
