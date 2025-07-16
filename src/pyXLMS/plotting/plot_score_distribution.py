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


def plot_score_distribution(
    data: List[Dict[str, Any]],
    bins: int = 25,
    density: bool = False,
    colors: List[str] = ["#00a087", "#3c5488", "#e64b35"],
    title: str = "Target and Decoy Score Distribution",
    figsize: Tuple[float, float] = (16.0, 9.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Plot the score distribution for a set of crosslink-spectrum-matches or crosslinks.

    Plot the target-target, target-decoy, and decoy-decoy score distribution as a histogram for a
    set of crosslink-spectrum-matches or crosslinks.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks.
    bins : int, default = 25
        The number of equal-width bins in the histogram.
    density : bool, default = False
        If True, draw and return a probability density: each bin will display the bin's raw count
        divided by the total number of counts and the bin width, so that the area under the histogram
        integrates to 1.
    colors : list of str, default = ["#00a087", "#3c5488", "#e64b35"]
        Colors of the histogram lines.
    title : str, default = "Target and Decoy Score Distribution"
        The title of the histogram.
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
        If attribute 'score', 'alpha_decoy', or 'beta_decoy' is not available for any of the data.

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS import plotting
    >>> pr = parser.read_msannika("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx")
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> fig, ax = plotting.plot_score_distribution(csms)
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(bins, "bins", int)
    _ok = check_input(density, "density", bool)
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
            "Can't plot score distribution if no crosslink-spectrum-matches or crosslinks are given!"
        )
    if "data_type" not in data[0] or data[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
        )
    available_keys = get_available_keys(data)
    if (
        not available_keys["score"]
        or not available_keys["alpha_decoy"]
        or not available_keys["beta_decoy"]
    ):
        raise ValueError(
            "Can't plot score distribution if 'score' or target/decoy labels are missing!"
        )
    ylabel = (
        "crosslink-spectrum-matches"
        if data[0]["data_type"] == "crosslink-spectrum-match"
        else "crosslinks"
    )
    filtered = filter_target_decoy(data)
    tt = [item["score"] for item in filtered["Target-Target"]]
    td = [item["score"] for item in filtered["Target-Decoy"]]
    dd = [item["score"] for item in filtered["Decoy-Decoy"]]

    fig, ax = plt.subplots(figsize=figsize)

    ax.hist(
        [tt, td, dd],
        bins=25,
        density=False,
        histtype="step",
        fill=False,
        color=colors,
        label=["Target-Target", "Target-Decoy", "Decoy-Decoy"],
    )
    ax.legend()
    ax.set_ylabel(f"Number of {ylabel}")
    ax.set_xlabel("Score")

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
