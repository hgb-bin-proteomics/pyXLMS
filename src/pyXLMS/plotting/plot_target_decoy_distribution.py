#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from ..data import check_input
from ..transform.util import get_available_keys
from ..transform.filter import filter_target_decoy

from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any


def plot_target_decoy_distribution(
    data: List[Dict[str, Any]],
    colors: List[str] = ["#00a087", "#3c5488", "#e64b35"],
    title: str = "Target and Decoy Distribution",
    figsize: Tuple[float, float] = (16.0, 9.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Plot the target-decoy distribution for a set of crosslink-spectrum-matches or crosslinks.

    Plot the target-target, target-decoy, and decoy-decoy distribution as a barplot for a
    set of crosslink-spectrum-matches or crosslinks.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks.
    colors : list of str, default = ["#00a087", "#3c5488", "#e64b35"]
        Colors of the bars.
    title : str, default = "Target and Decoy Distribution"
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
        If parameter data does not contain any crosslink-spectrum-matches or crosslinks.
    ValueError
        If attribute 'alpha_decoy', or 'beta_decoy' is not available for any of the data.

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS import plotting
    >>> pr = parser.read_msannika(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
    ... )
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> fig, ax = plotting.plot_target_decoy_distribution(csms)
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
    if len(data) == 0:
        raise ValueError(
            "Can't plot target-decoy distribution if no crosslink-spectrum-matches or crosslinks are given!"
        )
    if "data_type" not in data[0] or data[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
        )
    available_keys = get_available_keys(data)
    if not available_keys["alpha_decoy"] or not available_keys["beta_decoy"]:
        raise ValueError(
            "Can't plot target-decoy distribution if target/decoy labels are missing!"
        )
    ylabel = (
        "crosslink-spectrum-matches"
        if data[0]["data_type"] == "crosslink-spectrum-match"
        else "crosslinks"
    )
    filtered = filter_target_decoy(data)
    tt = len(filtered["Target-Target"])
    td = len(filtered["Target-Decoy"])
    dd = len(filtered["Decoy-Decoy"])

    fig, ax = plt.subplots(figsize=figsize)

    bar = ax.bar(
        ["Target-Target", "Target-Decoy", "Decoy-Decoy"],
        [tt, td, dd],
        color=colors,
    )
    ax.bar_label(bar, padding=3.0)
    ax.set_ylabel(f"Number of {ylabel}")
    ax.set_xlabel("Type")

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
