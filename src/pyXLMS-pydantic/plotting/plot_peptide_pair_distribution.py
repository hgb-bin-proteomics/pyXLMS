#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from ..data import check_input
from ..transform.filter import filter_peptide_pair_distribution

from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any


def plot_peptide_pair_distribution(
    data: List[Dict[str, Any]],
    top_n: int = 25,
    color: str = "#6d4bff",
    title: str = "Peptide Pair Distribution",
    figsize: Tuple[float, float] = (16.0, 9.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Plot the peptide pair distribution for a set of crosslink-spectrum-matches.

    Plot the peptide pair distribution as a barplot for a set of crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches.
    top_n : int, default = 25
        Number of peptide pairs to plot. Peptide pairs are sorted by number of
        crosslink-spectrum-matches.
    color : str, default = "#6d4bff"
        Color of the bars.
    title : str, default = "Peptide Pair Distribution"
        The title of the barplot.
    figsize : tuple of float, float, default = (16.0, 9.0)
        Width, height in inches.
    filename_prefix : str, or None
        If given, plot will be saved with and without title in .png and .svg format with the given
        prefix.

    Returns
    -------
    tuple of matplotlib.figure.Figure, any
        The created figure and axis ``from matplotlib.pyplot.subplots()``.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    ValueError
        If parameter data does not contain any crosslink-spectrum-matches.

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS import plotting
    >>> pr = parser.read_msannika(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
    ... )
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> fig, ax = plotting.plot_peptide_pair_distribution(csms)
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(top_n, "top_n", int)
    _ok = check_input(color, "color", str)
    _ok = check_input(title, "title", str)
    _ok = check_input(figsize, "figsize", tuple)
    _ok = (
        check_input(filename_prefix, "filename_prefix", str)
        if filename_prefix is not None
        else True
    )
    if len(data) == 0:
        raise ValueError(
            "Can't plot peptide pair distribution if no crosslink-spectrum-matches are given!"
        )
    if "data_type" not in data[0] or data[0]["data_type"] != "crosslink-spectrum-match":
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink-spectrum-match!"
        )
    pps = filter_peptide_pair_distribution(data)
    pp_names = list()
    pp_total = list()
    for pp in pps:
        pp_names.append(pp)
        pp_total.append(len(pps[pp]))

    sorted = pd.DataFrame(
        {
            "peptide_pair": pp_names,
            "total": pp_total,
        }
    ).sort_values(by="total", axis=0, ascending=False)
    pp_names = sorted["peptide_pair"].values.tolist()[:top_n]
    pp_total = sorted["total"].values.tolist()[:top_n]

    fig, ax = plt.subplots(figsize=figsize)

    bar = ax.bar(pp_names, pp_total, color=color)
    ax.bar_label(bar, padding=3.0)

    ax.set_xticks(range(len(pp_names)), pp_names, rotation=45, ha="right")
    ax.set_ylabel("Number of crosslink-spectrum-matches")
    ax.set_xlabel("Peptide Pair")

    if filename_prefix is not None:
        plt.savefig(
            filename_prefix + "_notitle.png",
            dpi=300,
            transparent=True,
            bbox_inches="tight",
        )
        plt.savefig(
            filename_prefix + "_notitle.svg",
            dpi=300,
            transparent=True,
            bbox_inches="tight",
        )
        ax.set_title(title)
        plt.savefig(
            filename_prefix + ".png", dpi=300, transparent=True, bbox_inches="tight"
        )
        plt.savefig(
            filename_prefix + ".svg", dpi=300, transparent=True, bbox_inches="tight"
        )
    else:
        ax.set_title(title)

    return (fig, ax)
