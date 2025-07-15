#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test1():
    from pyXLMS.plotting import venn

    fig, ax = venn(
        {"A", "B", "C"},
        {"B", "C", "D", "E", "F"},
        labels=["A", "F"],
        colors=["orange", "blue"],
    )
    assert fig is not None
    assert ax is not None


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test2():
    from pyXLMS.plotting import venn

    fig, ax = venn({"A", "B", "C"}, {"B", "C", "D", "E", "F"}, {"F", "G"})
    assert fig is not None
    assert ax is not None


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test3():
    from pyXLMS.plotting import venn

    with pytest.warns(
        RuntimeWarning,
        match=r"More than two labels supplied for two sets\. Using first two\.\.\.",
    ):
        _fig, _ax = venn(
            {"A", "B", "C"}, {"B", "C", "D", "E", "F"}, colors=["orange", "blue"]
        )


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test4():
    from pyXLMS.plotting import venn

    with pytest.warns(
        RuntimeWarning,
        match=r"More than two colors supplied for two sets\. Using first two\.\.\.",
    ):
        _fig, _ax = venn({"A", "B", "C"}, {"B", "C", "D", "E", "F"}, labels=["A", "F"])


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test5():
    from pyXLMS.plotting import venn

    with pytest.raises(
        IndexError,
        match=r"At least two labels have to be given if two sets are supplied!",
    ):
        _fig, _ax = venn(
            {"A", "B", "C"},
            {"B", "C", "D", "E", "F"},
            labels=["A"],
            colors=["orange", "blue"],
        )


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test6():
    from pyXLMS.plotting import venn

    with pytest.raises(
        IndexError,
        match=r"At least two colors have to be given if two sets are supplied!",
    ):
        _fig, _ax = venn(
            {"A", "B", "C"},
            {"B", "C", "D", "E", "F"},
            labels=["A", "F"],
            colors=["orange"],
        )


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test7():
    from pyXLMS.plotting import venn

    with pytest.warns(
        RuntimeWarning,
        match=r"More than three colors supplied for three sets\. Using first three\.\.\.",
    ):
        _fig, _ax = venn(
            {"A", "B", "C"},
            {"B", "C", "D", "E", "F"},
            {"F", "G"},
            colors=["orange", "blue", "green", "red"],
        )


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test8():
    from pyXLMS.plotting import venn

    with pytest.warns(
        RuntimeWarning,
        match=r"More than three labels supplied for three sets\. Using first three\.\.\.",
    ):
        _fig, _ax = venn(
            {"A", "B", "C"},
            {"B", "C", "D", "E", "F"},
            {"F", "G"},
            labels=["A", "B", "C", "D"],
        )


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test9():
    from pyXLMS.plotting import venn

    with pytest.raises(
        IndexError,
        match=r"At least three labels have to be given if three sets are supplied!",
    ):
        _fig, _ax = venn(
            {"A", "B", "C"}, {"B", "C", "D", "E", "F"}, {"F", "G"}, labels=["A", "B"]
        )


@pytest.mark.filterwarnings("ignore:'mode' parameter is deprecated")
def test10():
    from pyXLMS.plotting import venn

    with pytest.raises(
        IndexError,
        match=r"At least three colors have to be given if three sets are supplied!",
    ):
        _fig, _ax = venn(
            {"A", "B", "C"},
            {"B", "C", "D", "E", "F"},
            {"F", "G"},
            colors=["orange", "blue"],
        )
