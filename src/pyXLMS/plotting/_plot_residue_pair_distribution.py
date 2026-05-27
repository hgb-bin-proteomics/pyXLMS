#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from ..data._csm import CrosslinkSpectrumMatch
from ..data._util import check_input
from ..transform._filter import filter_residue_pair_distribution

from typing import Optional
from typing import List
from typing import Tuple
from typing import Any


def plot_residue_pair_distribution(
    data: List[CrosslinkSpectrumMatch],
    top_n: int = 25,
    color: str = "#6d4bff",
    title: str = "Residue Pair Distribution",
    figsize: Tuple[float, float] = (16.0, 9.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Plot the residue pair distribution for a set of crosslink-spectrum-matches.

    Plot the residue pair distribution as a barplot for a set of crosslink-spectrum-matches.
    Requires that ``alpha_proteins``, ``beta_proteins``, ``alpha_proteins_crosslink_positions``, and
    ``beta_proteins_crosslink_positions`` fields are set for all crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch
        A list of crosslink-spectrum-matches.
    top_n : int, default = 25
        Number of residue pairs to plot. Residue pairs are sorted by number of
        crosslink-spectrum-matches.
    color : str, default = "#6d4bff"
        Color of the bars.
    title : str, default = "Residue Pair Distribution"
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
    >>> fig, ax = plotting.plot_residue_pair_distribution(csms)
    """
    _ok = check_input(data, "data", list, CrosslinkSpectrumMatch)
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
            "Can't plot residue pair distribution if no crosslink-spectrum-matches are given!"
        )
    rps = filter_residue_pair_distribution(data)
    rp_names = list()
    rp_total = list()
    for rp in rps:
        rp_names.append(rp)
        rp_total.append(len(rps[rp]))

    sorted = pd.DataFrame(
        {
            "residue_pair": rp_names,
            "total": rp_total,
        }
    ).sort_values(by="total", axis=0, ascending=False)
    rp_names = sorted["residue_pair"].values.tolist()[:top_n]
    rp_total = sorted["total"].values.tolist()[:top_n]

    fig, ax = plt.subplots(figsize=figsize)

    bar = ax.bar(rp_names, rp_total, color=color)
    ax.bar_label(bar, padding=3.0)

    ax.set_xticks(range(len(rp_names)), rp_names, rotation=45, ha="right")
    ax.set_ylabel("Number of crosslink-spectrum-matches")
    ax.set_xlabel("Residue Pair")

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
