#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest
import pandas as pd


CROSSLINKS_ALL = "data/_test/exporter/pyxlinkviewer/unique_links_all.xlsx"
CROSSLINKS_NSP8 = "data/_test/exporter/pyxlinkviewer/unique_links_nsp8.xlsx"
YHU_PDB = "data/_test/exporter/pyxlinkviewer/6YHU/6yhu.pdb"
YHU_RESULT_ALL = (
    "data/_test/exporter/pyxlinkviewer/6YHU/all/unique_links_all_crosslinks.txt"
)
YHU_RESULT_NSP8 = (
    "data/_test/exporter/pyxlinkviewer/6YHU/nsp8/unique_links_nsp8_crosslinks.txt"
)
JLT_PDB = "data/_test/exporter/pyxlinkviewer/7JLT/7jlt.pdb"
JLT_RESULT_NSP8 = (
    "data/_test/exporter/pyxlinkviewer/7JLT/nsp8/unique_links_nsp8_7jlt_crosslinks.txt"
)


def helper_get_crosslink_position(pep_seq: str) -> int:
    for i, aa in enumerate(pep_seq.strip()):
        if aa == "[":
            return i


def helper_get_crosslinks(f: str) -> pd.DataFrame:
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.parser.util import format_sequence

    df = pd.read_excel(f)
    xls = list()
    for i, row in df.iterrows():
        xls.append(
            create_crosslink_min(
                format_sequence(row["Sequence A"]),
                helper_get_crosslink_position(row["Sequence A"]) + 1,
                format_sequence(row["Sequence B"]),
                helper_get_crosslink_position(row["Sequence B"]) + 1,
            )
        )
    return xls


def helper_get_dataframe_from_pyxlinkviewer(filename: str) -> pd.DataFrame:
    pyxlinkviewer = ""
    with open(filename, "r", encoding="utf-8") as f:
        pyxlinkviewer = f.read()
        f.close()
    pos_a = list()
    chain_a = list()
    pos_b = list()
    chain_b = list()
    for line in pyxlinkviewer.split("\n"):
        if len(line.strip()) > 0:
            xl = line.split("|")
            pos_a.append(int(xl[0]))
            chain_a.append(xl[1].strip())
            pos_b.append(int(xl[2]))
            chain_b.append(xl[3].strip())
    return pd.DataFrame(
        {"residue 1": pos_a, "chain 1": chain_a, "residue 2": pos_b, "chain 2": chain_b}
    )


def test1():
    import os
    from pyXLMS.exporter import to_pyxlinkviewer
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/pyxlinkviewer/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    pyxlinkviewer_result = to_pyxlinkviewer(
        crosslinks, pdb_file="6YHU", filename_prefix="6YHU"
    )
    assert pyxlinkviewer_result is not None
    pyxlinkviewer_output_file_str = pyxlinkviewer_result["PyXlinkViewer"]
    assert pyxlinkviewer_output_file_str is not None
    pyxlinkviewer_dataframe = pyxlinkviewer_result["PyXlinkViewer DataFrame"]
    assert pyxlinkviewer_dataframe is not None
    nr_mapped_crosslinks = pyxlinkviewer_result["Number of mapped crosslinks"]
    assert nr_mapped_crosslinks == 28
    crosslink_mapping = pyxlinkviewer_result["Mapping"]
    assert crosslink_mapping is not None
    parsed_pdb_sequenece = pyxlinkviewer_result["Parsed PDB sequence"]
    assert parsed_pdb_sequenece is not None
    parsed_pdb_chains = pyxlinkviewer_result["Parsed PDB chains"]
    assert parsed_pdb_chains is not None
    parsed_pdb_residue_numbers = pyxlinkviewer_result["Parsed PDB residue numbers"]
    assert parsed_pdb_residue_numbers is not None
    exported_files = pyxlinkviewer_result["Exported files"]
    assert len(exported_files) == 4
    for exported_file in exported_files:
        assert os.path.isfile(exported_file)


def test2():
    from pyXLMS.exporter import to_pyxlinkviewer
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/pyxlinkviewer/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    pyxlinkviewer_result = to_pyxlinkviewer(crosslinks, pdb_file="6YHU")
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_ALL)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test3():
    from pyXLMS.exporter import to_pyxlinkviewer
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/pyxlinkviewer/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    pyxlinkviewer_result = to_pyxlinkviewer(crosslinks, pdb_file=YHU_PDB)
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_ALL)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test4():
    from pyXLMS.exporter import to_pyxlinkviewer

    xls = helper_get_crosslinks(CROSSLINKS_ALL)
    pyxlinkviewer_result = to_pyxlinkviewer(xls, pdb_file=YHU_PDB)
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_ALL)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test5():
    from pyXLMS.exporter import to_pyxlinkviewer

    xls = helper_get_crosslinks(CROSSLINKS_ALL)
    pyxlinkviewer_result = to_pyxlinkviewer(xls, pdb_file="6YHU.pdb")
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_ALL)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test6():
    from pyXLMS.exporter import to_pyxlinkviewer

    xls = helper_get_crosslinks(CROSSLINKS_ALL)
    pyxlinkviewer_result = to_pyxlinkviewer(xls, pdb_file="6YHU")
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_ALL)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test7():
    from pyXLMS.exporter import to_pyxlinkviewer

    xls = helper_get_crosslinks(CROSSLINKS_ALL)
    pyxlinkviewer_result = to_pyxlinkviewer(xls, pdb_file="6yhu")
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_ALL)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test8():
    from pyXLMS.exporter import to_pyxlinkviewer

    xls = helper_get_crosslinks(CROSSLINKS_NSP8)
    pyxlinkviewer_result = to_pyxlinkviewer(xls, pdb_file="6yhu")
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_NSP8)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test9():
    from pyXLMS.exporter import to_pyxlinkviewer

    xls = helper_get_crosslinks(CROSSLINKS_NSP8)
    pyxlinkviewer_result = to_pyxlinkviewer(xls, pdb_file=JLT_PDB)
    should_result = helper_get_dataframe_from_pyxlinkviewer(JLT_RESULT_NSP8)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test10():
    from pyXLMS.exporter import to_pyxlinkviewer

    xls = helper_get_crosslinks(CROSSLINKS_NSP8)
    pyxlinkviewer_result = to_pyxlinkviewer(xls, pdb_file="7jlt")
    should_result = helper_get_dataframe_from_pyxlinkviewer(JLT_RESULT_NSP8)
    assert should_result.equals(pyxlinkviewer_result["PyXlinkViewer DataFrame"])


def test11():
    from pyXLMS.exporter import to_pyxlinkviewer
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/pyxlinkviewer/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    with pytest.raises(
        ValueError,
        match=r"Minimum sequence identity should be given as a fraction, e\.g\. 0\.8 for 80\% minimum sequence identity!",
    ):
        _pyxlinkviewer_result = to_pyxlinkviewer(
            crosslinks, pdb_file="6YHU", min_sequence_identity=80.0
        )


def test12():
    from pyXLMS.exporter import to_pyxlinkviewer

    with pytest.raises(
        ValueError, match=r"Provided data does not contain any crosslinks!"
    ):
        _pyxlinkviewer_result = to_pyxlinkviewer([], pdb_file="6YHU")


def test13():
    from pyXLMS.exporter import to_pyxlinkviewer
    from pyXLMS.data import create_csm_min

    csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    csms = [csm1]
    with pytest.raises(
        TypeError,
        match="Unsupported data type for input data! Parameter data has to be a list of crosslinks!",
    ):
        _pyxlinkviewer_result = to_pyxlinkviewer(csms, pdb_file="6YHU")


def test14():
    from pyXLMS.exporter import to_pyxlinkviewer
    from pyXLMS.data import create_csm_min
    from pyXLMS.data import create_crosslink_min

    csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    xl1 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    mixed = [xl1, csm1]
    with pytest.raises(
        TypeError, match="Not all elements in data have the same data type!"
    ):
        _pyxlinkviewer_result = to_pyxlinkviewer(mixed, pdb_file="6YHU")
