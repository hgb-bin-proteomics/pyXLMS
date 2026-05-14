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

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def plot_crosslink_type_distribution(
    data: List[Dict[str, Any]],
    plot_type: Literal["bar", "pie"] = "bar",
    colors: List[str] = ["#6d4bff", "#ac99ff"],
    title: str = "Crosslink Type Distribution",
    figsize: Tuple[float, float] = (16.0, 9.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Plot the crosslink type distribution for a set of crosslink-spectrum-matches or crosslinks.

    Plot the crosslink type distribution (intra- and inter-links) as a bar or pie chart for a set of
    crosslink-spectrum-matches or crosslinks.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks.
    plot_type : str, one of "bar" or "pie", default = "bar"
        Plot type, whether to plot as a bar or pie chart.
    colors : list of str, default = ["#6d4bff", "#ac99ff"]
        Colors of the bars/pie slices (intra-link and inter-link).
    title : str, default = "Crosslink Type Distribution"
        The title of the plot.
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
        If parameter data does not contain any crosslink-spectrum-matches or crosslinks.
    ValueError
        If parameter plot type was set incorrectly.
    IndexError
        If not enough colors where specified.

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS import plotting
    >>> pr = parser.read_msannika(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
    ... )
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> fig, ax = plotting.plot_crosslink_type_distribution(csms)
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(plot_type, "plot_type", str)
    _ok = check_input(colors, "colors", list, str)
    _ok = check_input(title, "title", str)
    _ok = check_input(figsize, "figsize", tuple)
    _ok = (
        check_input(filename_prefix, "filename_prefix", str)
        if filename_prefix is not None
        else True
    )
    if plot_type not in ["bar", "pie"]:
        raise ValueError("Plot type needs to be one of 'bar' or 'pie'!")
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
    axis_label = (
        "crosslink-spectrum-matches"
        if data[0]["data_type"] == "crosslink-spectrum-match"
        else "crosslinks"
    )
    intra_inter = filter_crosslink_type(data)
    values = [len(intra_inter["Intra"]), len(intra_inter["Inter"])]
    labels = ["intra-links", "inter-links"]

    fig, ax = plt.subplots(figsize=figsize)

    if plot_type == "pie":
        ax.pie(
            values,
            labels=labels,
            colors=colors,
            autopct="%1.1f%%",
        )

        ax.set_xlabel(
            f"Total number of {axis_label}: {sum([len(intra_inter['Intra']), len(intra_inter['Inter'])])}"
        )
    else:
        bar = ax.bar(labels, values, color=colors)
        ax.bar_label(bar, padding=3.0)

        ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
        ax.set_ylabel(f"Number of {axis_label}")
        ax.set_xlabel("Crosslink Type")

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
