#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._util import check_input
from ..data._util import check_input_multi
from ..transform._util import assert_csms_or_xls
from ..transform._filter import filter_target_decoy
from ..transform._filter import filter_crosslink_type
from ..transform._annotate_string_scores import annotate_string_scores
from ..transform._annotate_string_scores import STRING_SCORES

from typing import Optional
from typing import List
from typing import Tuple
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def plot_string_score_distribution(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink],
    organism: Optional[str | int] = None,
    plot_type: Literal["bar", "hist"] = "bar",
    bins: int = 25,
    density: bool = False,
    zero_impute_nan: bool = True,
    colors: Optional[List[str]] = None,
    title: str = "STRING Score Distribution for Inter-Links",
    figsize: Tuple[float, float] = (16.0, 9.0),
    filename_prefix: Optional[str] = None,
    verbose: Literal[0, 1, 2] = 1,
) -> Tuple[Figure, Any]:
    r"""Plot the STRING score distribution for a set of inter-links.

    Plot the STRING score distribution as a barplot or histogram for inter-links of a
    set of crosslink-spectrum-matches or crosslinks.
    STRING is accessible via `string-db.org <https://string-db.org>`_.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslink-spectrum-matches or crosslinks.
    organism : str, or int, or None, default = None
        Organism name (e.g. Homo sapiens) or taxon identifier (e.g. 9606).
        Taxon identifiers are preferred. See also
        `string-db.org/cgi/organisms <https://string-db.org/cgi/organisms>`_.
        If ``None`` it is assumed that the input data is already annotated with STRING
        scores and will raise an error if that is not the case.
    plot_type : "bar", or "hist", default = "bar"
        If STRING scores should be plotted as a bar plot or as a histogram.
    bins : int, default = 25
        The number of equal-width bins in the histogram.
        Only applies to ``plot_type = "hist"``.
    density : bool, default = False
        If True, draw and return a probability density: each bin will display the bin's raw count
        divided by the total number of counts and the bin width, so that the area under the histogram
        integrates to 1.
        Only applies to ``plot_type = "hist"``.
    zero_impute_nan : bool, default = True
        If nan values should be imputed with zeros.
        Only applies to ``plot_type = "hist"``.
    colors : list of str, or None, default = None
        Colors of the bars. For ``plot_type = "bar"`` a total number of ``6`` colors have to be given.
        For ``plot_type = "hist"`` a total number of ``3`` colors have to be given. Uses the internal
        defaults if ``None`` (default) is given.
    title : str, default = "STRING Score Distribution for Inter-Links"
        The title of the plot.
    figsize : tuple of float, float, default = (16.0, 9.0)
        Width, height in inches.
    filename_prefix : str, or None
        If given, plot will be saved with and without title in .png and .svg format with the given
        prefix.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

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
        If parameter data does not contain any inter-links.
    ValueError
        If the number of given colors does not match the number of required colors for
        the plot type.
    ValueError
        If organism is None and data is not yet annotated with STRING scores.
    TypeError
        If parameter plot_type was not set correctly.
    TypeError
        If parameter verbose was not set correctly.

    Notes
    -----
    It is generally recommended to call ``transform.annotate_string_scores()`` before
    using this function to preemptively catch errors during annotation.

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS import plotting
    >>> pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    >>> xls = pr["crosslinks"]
    >>> fig, ax = plotting.plot_string_score_distribution(xls, organism="Homo sapiens")
    """
    _ok = check_input(data, "data", list)
    _ok = (
        check_input_multi(organism, "organism", [str, int])
        if organism is not None
        else True
    )
    _ok = check_input(plot_type, "plot_type", str)
    if plot_type not in ["bar", "hist"]:
        raise TypeError("Plot type has to be one of 'bar', or 'hist'!")
    _ok = check_input(colors, "colors", list, str) if colors is not None else True
    if colors is not None:
        if plot_type == "bar":
            if len(colors) != 6:
                raise ValueError("Six colors are required for plot type 'bar'!")
        else:
            if len(colors) != 3:
                raise ValueError("Three colors are required for plot type 'hist'!")
    else:
        if plot_type == "bar":
            colors = ["#e64b35", "#4dbbd5", "#00a087", "#3c5488", "#f39b7f", "#8491b4"]
        else:
            colors = ["#00a087", "#3c5488", "#e64b35"]
    _ok = check_input(title, "title", str)
    _ok = check_input(figsize, "figsize", tuple)
    _ok = (
        check_input(filename_prefix, "filename_prefix", str)
        if filename_prefix is not None
        else True
    )
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")
    if len(data) == 0:
        raise ValueError(
            "Can't plot STRING score distribution if no crosslink-spectrum-matches or crosslinks are given!"
        )
    data = assert_csms_or_xls(data)
    inter = filter_crosslink_type(data)["Inter"]
    if len(inter) == 0:
        raise ValueError(
            "Can't plot STRING score distribution because data does not contain inter-links!"
        )
    if (
        "additional_information" not in inter[0]
        or inter[0]["additional_information"] is None
        or "pyXLMS_annotated_STRING_score" not in inter[0]["additional_information"]
    ):
        if organism is None:
            raise ValueError(
                "Input data does not have annotated STRING scores! In this case a valid organism has to be given!"
            )
        _ = annotate_string_scores(inter, organism, verbose)
    ylabel = (
        "crosslink-spectrum-matches"
        if data[0]["data_type"] == "crosslink-spectrum-match"
        else "crosslinks"
    )

    fig, ax = plt.subplots(figsize=figsize)

    if plot_type == "bar":
        missing_scores = 0
        lowest_conf = 0
        low_conf = 0
        medium_conf = 0
        high_conf = 0
        highest_conf = 0
        for item in inter:
            score = item["additional_information"]["pyXLMS_annotated_STRING_score"]
            if pd.isna(score):
                missing_scores += 1
            elif score >= STRING_SCORES["highest confidence"]:
                highest_conf += 1
            elif score >= STRING_SCORES["high confidence"]:
                high_conf += 1
            elif score >= STRING_SCORES["medium confidence"]:
                medium_conf += 1
            elif score >= STRING_SCORES["low confidence"]:
                low_conf += 1
            else:
                lowest_conf += 1

        x = [
            "Missing STRING score",
            "Lowest Confidence",
            "Low Confidence",
            "Medium Confidence",
            "High Confidence",
            "Highest Confidence",
        ]
        y = [
            missing_scores,
            lowest_conf,
            low_conf,
            medium_conf,
            high_conf,
            highest_conf,
        ]
        bar = ax.bar(
            x,
            y,
            color=colors,
        )
        ax.bar_label(bar, label_type="center")
        ax.set_xticks(range(len(x)), x, rotation=45, ha="right")
        ax.set_ylabel(f"Number of {ylabel}")
        ax.set_xlabel("STRING Score Type")
    else:
        filtered = filter_target_decoy(inter)
        if zero_impute_nan:
            tt = [
                item["additional_information"]["pyXLMS_annotated_STRING_score"]
                if not pd.isna(
                    item["additional_information"]["pyXLMS_annotated_STRING_score"]
                )
                else 0.0
                for item in filtered["Target-Target"]
            ]
            td = [
                item["additional_information"]["pyXLMS_annotated_STRING_score"]
                if not pd.isna(
                    item["additional_information"]["pyXLMS_annotated_STRING_score"]
                )
                else 0.0
                for item in filtered["Target-Decoy"]
            ]
            dd = [
                item["additional_information"]["pyXLMS_annotated_STRING_score"]
                if not pd.isna(
                    item["additional_information"]["pyXLMS_annotated_STRING_score"]
                )
                else 0.0
                for item in filtered["Decoy-Decoy"]
            ]
        else:
            tt = [
                item["additional_information"]["pyXLMS_annotated_STRING_score"]
                for item in filtered["Target-Target"]
                if not pd.isna(
                    item["additional_information"]["pyXLMS_annotated_STRING_score"]
                )
            ]
            td = [
                item["additional_information"]["pyXLMS_annotated_STRING_score"]
                for item in filtered["Target-Decoy"]
                if not pd.isna(
                    item["additional_information"]["pyXLMS_annotated_STRING_score"]
                )
            ]
            dd = [
                item["additional_information"]["pyXLMS_annotated_STRING_score"]
                for item in filtered["Decoy-Decoy"]
                if not pd.isna(
                    item["additional_information"]["pyXLMS_annotated_STRING_score"]
                )
            ]

        ax.hist(
            [tt, td, dd],
            bins=bins,
            density=density,
            histtype="step",
            fill=False,
            color=colors,
            label=["Target-Target", "Target-Decoy", "Decoy-Decoy"],
        )
        ax.legend(loc="upper right")
        ax.set_ylabel(f"Number of {ylabel}")
        ax.set_xlabel("STRING Score")

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
