#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest
import pandas as pd


YHU_PDB = "data/_test/exporter/pyxlinkviewer/6YHU/6yhu.pdb"
YHU_RESULT_ALL = (
    "data/_test/exporter/pyxlinkviewer/6YHU/all/unique_links_all_crosslinks.txt"
)


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
    from pyXLMS.exporter import to_xlmstools
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/xlms-tools/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    xlmstools_result = to_xlmstools(crosslinks, pdb_file="6YHU", filename_prefix="6YHU")
    xlmstools_output_file_str = xlmstools_result["xlms-tools"]
    xlmstools_dataframe = xlmstools_result["xlms-tools DataFrame"]
    nr_mapped_crosslinks = xlmstools_result["Number of mapped crosslinks"]
    crosslink_mapping = xlmstools_result["Mapping"]
    parsed_pdb_sequenece = xlmstools_result["Parsed PDB sequence"]
    parsed_pdb_chains = xlmstools_result["Parsed PDB chains"]
    parsed_pdb_residue_numbers = xlmstools_result["Parsed PDB residue numbers"]
    exported_files = xlmstools_result["Exported files"]
    assert exported_files[0] == "6YHU_xlms-tools.txt"
    assert xlmstools_output_file_str is not None
    assert xlmstools_dataframe is not None
    assert nr_mapped_crosslinks is not None
    assert crosslink_mapping is not None
    assert parsed_pdb_sequenece is not None
    assert parsed_pdb_chains is not None
    assert parsed_pdb_residue_numbers is not None


def test2():
    import os
    from pyXLMS.exporter import to_xlmstools
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/xlms-tools/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    xlmstools_result = to_xlmstools(crosslinks, pdb_file="6YHU", filename_prefix="test")
    assert xlmstools_result is not None
    xlmstools_output_file_str = xlmstools_result["xlms-tools"]
    assert xlmstools_output_file_str is not None
    xlmstools_dataframe = xlmstools_result["xlms-tools DataFrame"]
    assert xlmstools_dataframe is not None
    nr_mapped_crosslinks = xlmstools_result["Number of mapped crosslinks"]
    assert nr_mapped_crosslinks == 28
    crosslink_mapping = xlmstools_result["Mapping"]
    assert crosslink_mapping is not None
    parsed_pdb_sequenece = xlmstools_result["Parsed PDB sequence"]
    assert parsed_pdb_sequenece is not None
    parsed_pdb_chains = xlmstools_result["Parsed PDB chains"]
    assert parsed_pdb_chains is not None
    parsed_pdb_residue_numbers = xlmstools_result["Parsed PDB residue numbers"]
    assert parsed_pdb_residue_numbers is not None
    exported_files = xlmstools_result["Exported files"]
    assert len(exported_files) == 4
    assert exported_files[0] == "test_xlms-tools.txt"
    for exported_file in exported_files:
        assert os.path.isfile(exported_file)


def test3():
    from pyXLMS.exporter import to_xlmstools
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/xlms-tools/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    xlmstools_result = to_xlmstools(crosslinks, pdb_file="6YHU")
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_ALL)
    assert should_result.equals(xlmstools_result["xlms-tools DataFrame"])


def test4():
    from pyXLMS.exporter import to_xlmstools
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/xlms-tools/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    xlmstools_result = to_xlmstools(crosslinks, pdb_file=YHU_PDB)
    should_result = helper_get_dataframe_from_pyxlinkviewer(YHU_RESULT_ALL)
    assert should_result.equals(xlmstools_result["xlms-tools DataFrame"])


def test5():
    from pyXLMS.exporter import to_xlmstools
    from pyXLMS.parser import read_custom

    pr = read_custom("data/_test/exporter/xlms-tools/unique_links_all_pyxlms.csv")
    crosslinks = pr["crosslinks"]
    with pytest.raises(
        ValueError,
        match=r"Minimum sequence identity should be given as a fraction, e\.g\. 0\.8 for 80\% minimum sequence identity!",
    ):
        _xlmstools_result = to_xlmstools(
            crosslinks, pdb_file="6YHU", min_sequence_identity=80.0
        )


def test6():
    from pyXLMS.exporter import to_xlmstools

    with pytest.raises(
        ValueError, match=r"Provided data does not contain any crosslinks!"
    ):
        _xlmstools_result = to_xlmstools([], pdb_file="6YHU")


def test7():
    from pyXLMS.exporter import to_xlmstools
    from pyXLMS.data import create_csm_min

    csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    csms = [csm1]
    with pytest.raises(
        TypeError,
        match="Unsupported data type for input data! Parameter data has to be a list of crosslinks!",
    ):
        _xlmstools_result = to_xlmstools(csms, pdb_file="6YHU")


def test14():
    from pyXLMS.exporter import to_xlmstools
    from pyXLMS.data import create_csm_min
    from pyXLMS.data import create_crosslink_min

    csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    xl1 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    mixed = [xl1, csm1]
    with pytest.raises(
        TypeError, match="Not all elements in data have the same data type!"
    ):
        _xlmstools_result = to_xlmstools(mixed, pdb_file="6YHU")
