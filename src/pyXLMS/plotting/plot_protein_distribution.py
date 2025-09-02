#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from ..data import check_input
from ..transform.util import get_available_keys
from ..transform.filter import filter_protein_distribution
from ..transform.filter import filter_crosslink_type

from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any


def plot_protein_distribution(
    data: List[Dict[str, Any]],
    top_n: int = 25,
    colors: List[str] = ["#6d4bff", "#ac99ff"],
    title: str = "Protein Distribution",
    figsize: Tuple[float, float] = (16.0, 9.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Plot the protein distribution for a set of crosslink-spectrum-matches or crosslinks.

    Plot the protein distribution as a barplot for a set of crosslink-spectrum-matches or crosslinks.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks.
    top_n : int, default = 25
        Number of proteins to plot. Proteins are sorted by number of crosslinks
        or crosslink-spectrum-matches.
    colors : list of str, default = ["#6d4bff", "#ac99ff"]
        Colors of the bar-types (intra-link and inter-link).
    title : str, default = "Protein Distribution"
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
        If attribute 'alpha_proteins', or 'beta_proteins' is not available for any of the data.
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
    >>> fig, ax = plotting.plot_protein_distribution(csms)
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(top_n, "top_n", int)
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
            "Can't plot protein distribution if no crosslink-spectrum-matches or crosslinks are given!"
        )
    if "data_type" not in data[0] or data[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
        )
    available_keys = get_available_keys(data)
    if not available_keys["alpha_proteins"] or not available_keys["beta_proteins"]:
        raise ValueError(
            "Can't plot protein distribution if attributes 'alpha_proteins' or 'beta_proteins' are missing!"
        )
    ylabel = (
        "crosslink-spectrum-matches"
        if data[0]["data_type"] == "crosslink-spectrum-match"
        else "crosslinks"
    )
    proteins = filter_protein_distribution(data)
    protein_names = list()
    protein_intra = list()
    protein_inter = list()
    protein_total = list()
    for protein in proteins:
        protein_names.append(protein)
        intra_inter = filter_crosslink_type(proteins[protein])
        protein_intra.append(len(intra_inter["Intra"]))
        protein_inter.append(len(intra_inter["Inter"]))
        protein_total.append(len(proteins[protein]))

    sorted = pd.DataFrame(
        {
            "protein": protein_names,
            "intra": protein_intra,
            "inter": protein_inter,
            "total": protein_total,
        }
    ).sort_values(by="total", axis=0, ascending=False)
    protein_names = sorted["protein"].values.tolist()[:top_n]
    protein_intra = sorted["intra"].values.tolist()[:top_n]
    protein_inter = sorted["inter"].values.tolist()[:top_n]

    fig, ax = plt.subplots(figsize=figsize)
    bottom = [0.0 for i in protein_names]

    bar = ax.bar(
        protein_names, protein_intra, label="intra-link", bottom=bottom, color=colors[0]
    )
    bottom = protein_intra
    ax.bar_label(bar, label_type="center")

    bar = ax.bar(
        protein_names, protein_inter, label="inter-link", bottom=bottom, color=colors[1]
    )
    bottom = protein_inter
    ax.bar_label(bar, label_type="center")

    ax.set_xticks(range(len(protein_names)), protein_names, rotation=45, ha="right")
    ax.legend(loc="upper right")
    ax.set_ylabel(f"Number of {ylabel}")
    ax.set_xlabel("Protein")

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
