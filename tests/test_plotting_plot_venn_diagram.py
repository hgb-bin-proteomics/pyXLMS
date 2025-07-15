#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test1():
    from pyXLMS import parser
    from pyXLMS import plotting

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]
    fig, ax = plotting.plot_venn_diagram(
        a, b, labels=["MS Annika", "MaxQuant"], colors=["orange", "blue"]
    )
    assert fig is not None
    assert ax is not None


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test2():
    from pyXLMS import parser
    from pyXLMS import plotting

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]
    c = parser.read(
        "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv",
        engine="pLink",
        crosslinker="DSS",
    )
    c = c["crosslink-spectrum-matches"]
    fig, ax = plotting.plot_venn_diagram(
        a, b, c, labels=["MS Annika", "MaxQuant", "pLink"], contour=True
    )
    assert fig is not None
    assert ax is not None


def test3():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]

    err_str = r"Parameter 'by' has to be one of 'peptide' or 'protein'! Option 'peptide' will group by peptide sequence and "
    err_str += r"peptide crosslink position while option 'protein' will group by protein identifier and protein crosslink position."
    with pytest.raises(TypeError, match=err_str):
        _plot = plot_venn_diagram(a, b, by="sequence")


def test4():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]

    with pytest.raises(
        ValueError,
        match=r"Can't plot venn diagram if no crosslink-spectrum-matches or crosslinks are given in data_1!",
    ):
        _plot = plot_venn_diagram([], a)


def test5():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]

    with pytest.raises(
        ValueError,
        match=r"Can't plot venn diagram if no crosslink-spectrum-matches or crosslinks are given in data_2!",
    ):
        _plot = plot_venn_diagram(a, [])


def test6():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]

    with pytest.raises(
        ValueError,
        match=r"Can't plot 3-set venn diagram if no crosslink-spectrum-matches or crosslinks are given in data_3!",
    ):
        _plot = plot_venn_diagram(a, b, [])


def test7():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]
    a[0]["data_type"] = "peptide-spectrum-match"

    with pytest.raises(
        TypeError,
        match=r"Unsupported data type for input data! Parameter data_1 has to be a list of crosslink or crosslink-spectrum-match!",
    ):
        _plot = plot_venn_diagram(a, b)


def test8():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]
    a[0]["data_type"] = "peptide-spectrum-match"

    with pytest.raises(
        TypeError,
        match=r"Unsupported data type for input data! Parameter data_2 has to be a list of crosslink or crosslink-spectrum-match!",
    ):
        _plot = plot_venn_diagram(b, a)


def test9():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]
    c = parser.read(
        "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv",
        engine="pLink",
        crosslinker="DSS",
    )
    c = c["crosslink-spectrum-matches"]
    a[0]["data_type"] = "peptide-spectrum-match"

    with pytest.raises(
        TypeError,
        match=r"Unsupported data type for input data! Parameter data_3 has to be a list of crosslink or crosslink-spectrum-match, or None!",
    ):
        _plot = plot_venn_diagram(c, b, a)


def test10():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]
    c = parser.read(
        "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv",
        engine="pLink",
        crosslinker="DSS",
    )
    c = c["crosslink-spectrum-matches"]
    a[0]["completeness"] = "partial"
    a[0]["alpha_proteins"] = None

    err_str = r"Grouping by protein crosslink position is only available if all data have defined protein crosslink positions!\n"
    err_str += (
        r"This error might be fixable with 'transform\.reannotate_positions\(\)'!"
    )
    with pytest.raises(
        ValueError,
        match=err_str,
    ):
        _plot = plot_venn_diagram(c, b, a, by="protein")


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test11():
    from pyXLMS import parser
    from pyXLMS.plotting import plot_venn_diagram

    a = parser.read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
        engine="MS Annika",
        crosslinker="DSS",
    )
    a = a["crosslink-spectrum-matches"]
    b = parser.read(
        "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    )
    b = b["crosslink-spectrum-matches"]
    c = parser.read(
        "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv",
        engine="pLink",
        crosslinker="DSS",
    )
    c = c["crosslink-spectrum-matches"]

    plot = plot_venn_diagram(c, b, a, by="protein")
    assert plot is not None
