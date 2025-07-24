#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from ..data import check_input
from ..transform.filter import filter_crosslink_type

from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any


def plot_crosslink_type_distribution(
    data: List[Dict[str, Any]],
    colors: List[str] = ["#6d4bff", "#ac99ff"],
    title: str = "Crosslink Type Distribution",
    figsize: Tuple[float, float] = (10.0, 10.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Plot the crosslink type distribution for a set of crosslink-spectrum-matches or crosslinks.

    Plot the crosslink type distribution (intra- and inter-links) as a pie chart for a set of
    crosslink-spectrum-matches or crosslinks.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks.
    colors : list of str, default = ["#6d4bff", "#ac99ff"]
        Colors of the pie slices (intra-link and inter-link).
    title : str, default = "Crosslink Type Distribution"
        The title of the pie chart.
    figsize : tuple of float, float, default = (10.0, 10.0)
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
        If parameter data does not contain any crosslink-spectrum-matches or crosslinks.
    IndexError
        If not enough colors where specified.

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS import plotting
    >>> pr = parser.read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> fig, ax = plotting.plot_crosslink_type_distribution(csms)
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(colors, "colors", list, str)
    _ok = check_input(title, "title", str)
    _ok = check_input(figsize, "figsize", tuple)
    _ok = (
        check_input(filename_prefix, "filename_prefix", str)
        if filename_prefix is not None
        else True
    )
    if len(colors) < 2:
        raise IndexError("At least two colors need to be given for the plot!")
    if len(data) == 0:
        raise ValueError(
            "Can't plot crosslink type distribution if no crosslink-spectrum-matches or crosslinks are given!"
        )
    if "data_type" not in data[0] or data[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
        )
    xlabel = (
        "crosslink-spectrum-matches"
        if data[0]["data_type"] == "crosslink-spectrum-match"
        else "crosslinks"
    )
    intra_inter = filter_crosslink_type(data)

    fig, ax = plt.subplots(figsize=figsize)

    ax.pie(
        [len(intra_inter["Intra"]), len(intra_inter["Inter"])],
        labels=["intra-links", "inter-links"],
        colors=colors,
        autopct="%1.1f%%",
    )

    ax.set_xlabel(
        f"Total number of {xlabel}: {sum([len(intra_inter['Intra']), len(intra_inter['Inter'])])}"
    )

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
