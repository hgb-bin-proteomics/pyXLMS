#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles

from ..data import check_input
from ..transform.util import get_available_keys

from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Set
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __get_key(data: Dict[str, Any], by: Literal["peptide", "protein"]) -> str:
    r"""Get the unique key for a crosslink-spectrum-match or crosslink.

    Parameters
    ----------
    data : dict of str, any
        A pyXLMS crosslink-spectrum-match or crosslink object.
    by : str, one of "peptide" or "protein"
        If peptide or protein crosslink position should be used for determining if a crosslink-spectrum-match
        or crosslink is unique.

    Returns
    -------
    str
        The unique key for the crosslink-spectrum-match or crosslink.

    Notes
    -----
    This function should not be called directly, it is called from ``plot_venn_diagram()``.
    """
    if by == "peptide":
        return f"{data['alpha_peptide']}_{data['alpha_peptide_crosslink_position']}-{data['beta_peptide']}_{data['beta_peptide_crosslink_position']}"
    prot_pos_a = "-".join(
        sorted(
            [
                f"{data['alpha_proteins'][i]}_{data['alpha_proteins_crosslink_positions'][i]}"
                for i in range(len(data["alpha_proteins"]))
            ]
        )
    )
    prot_pos_b = "-".join(
        sorted(
            [
                f"{data['beta_proteins'][i]}_{data['beta_proteins_crosslink_positions'][i]}"
                for i in range(len(data["beta_proteins"]))
            ]
        )
    )
    return ":".join(sorted([prot_pos_a, prot_pos_b]))


def venn(
    set_1: Set[Any],
    set_2: Set[Any],
    set_3: Optional[Set[Any]] = None,
    labels: List[str] = ["Set 1", "Set 2", "Set 3"],
    colors: List[str] = ["#4361EE", "#4CC9F0", "#F72585"],
    alpha: float = 0.6,
    contour: bool = False,
    linewidth: float = 0.5,
    title: str = "Venn Diagram",
    figsize: Tuple[float, float] = (10.0, 10.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Wrapper with pre-set defaults for creating venn diagrams with the matplotlib-venn package.

    Wrapper with pre-set defaults for creating venn diagrams with the matplotlib-venn package
    `github.com/konstantint/matplotlib-venn <https://github.com/konstantint/matplotlib-venn>`_.

    Parameters
    ----------
    set_1 : set
        First set of the venn diagram.
    set_2 : set
        Second set of the venn diagram.
    set_3 : set, or None, default = None
        If not None a three set venn diagram will be drawn, if None
        the two set venn diagram of ``set_1`` and ``set_2`` will be drawn.
    labels : List[str], default = ["Set 1", "Set 2", "Set 3"]
        List of labels for the sets.
    colors : List[str], default = ["#4361EE", "#4CC9F0", "#F72585"]
        List of valid colors to use for the venn circles.
    alpha : float, default = 0.6
        Color opacity.
    contour : bool, default = False
        If a contour should be drawn around venn circles.
    linewidth : float, default = 0.5
        Linewidth of the contour.
    title : str, default = "Venn Diagram"
        Title of the venn diagram.
    figsize : tuple of float, float, default = (10.0, 10.0)
        Width, height in inches.
    filename_prefix: str, or None, default = None
        If given, plot will be saved with and without title in .png and .svg format with the given
        prefix.

    Returns
    -------
    tuple of matplotlib.figure.Figure, any
        The created figure and axis ``from matplotlib.pyplot.subplots()``.

    Warns
    -----
    RuntimeWarning
        If more labels or colors than sets are supplied.

    Raises
    ------
    IndexError
        If less labels or colors than sets are supplied.

    Examples
    --------
    >>> from pyXLMS.plotting import venn
    >>> fig, ax = venn({"A", "B", "C"}, {"B", "C", "D", "E", "F"}, labels=["A", "F"], colors=["orange", "blue"])

    >>> from pyXLMS.plotting import venn
    >>> fig, ax = venn({"A", "B", "C"}, {"B", "C", "D", "E", "F"}, {"F", "G"})
    """

    fig, ax = plt.subplots(figsize=figsize)

    if set_3 is None:
        # checks
        if len(labels) > 2:
            warnings.warn(
                RuntimeWarning(
                    "More than two labels supplied for two sets. Using first two..."
                )
            )
            labels = labels[:2]

        if len(labels) < 2:
            raise IndexError(
                "At least two labels have to be given if two sets are supplied!"
            )

        if len(colors) > 2:
            warnings.warn(
                RuntimeWarning(
                    "More than two colors supplied for two sets. Using first two..."
                )
            )
            colors = colors[:2]

        if len(colors) < 2:
            raise IndexError(
                "At least two colors have to be given if two sets are supplied!"
            )

        # create 2 set venn diagram
        venn2(
            subsets=(
                len(set_1.difference(set_2)),
                len(set_2.difference(set_1)),
                len(set_1.intersection(set_2)),
            ),
            set_labels=tuple(labels),  # pyright: ignore[reportArgumentType]
            set_colors=tuple(colors),  # pyright: ignore[reportArgumentType]
            alpha=alpha,
        )

        # draw contour
        if contour:
            venn2_circles(
                subsets=(
                    len(set_1.difference(set_2)),
                    len(set_2.difference(set_1)),
                    len(set_1.intersection(set_2)),
                ),
                linewidth=linewidth,
            )

        # save file
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
            plt.title(title)
            plt.savefig(
                filename_prefix + ".png", dpi=300, transparent=True, bbox_inches="tight"
            )
            plt.savefig(
                filename_prefix + ".svg", dpi=300, transparent=True, bbox_inches="tight"
            )
        else:
            plt.title(title)

    else:
        # checks
        if len(labels) > 3:
            warnings.warn(
                RuntimeWarning(
                    "More than three labels supplied for three sets. Using first three..."
                )
            )
            labels = labels[:3]

        if len(labels) < 3:
            raise IndexError(
                "At least three labels have to be given if three sets are supplied!"
            )

        if len(colors) > 3:
            warnings.warn(
                RuntimeWarning(
                    "More than three colors supplied for three sets. Using first three..."
                )
            )
            colors = colors[:3]

        if len(colors) < 3:
            raise IndexError(
                "At least three colors have to be given if three sets are supplied!"
            )

        # create 3 set venn diagram
        venn3(
            subsets=(
                len(set_1.difference(set_2, set_3)),
                len(set_2.difference(set_1, set_3)),
                len(set_1.intersection(set_2).difference(set_3)),
                len(set_3.difference(set_1, set_2)),
                len(set_1.intersection(set_3).difference(set_2)),
                len(set_2.intersection(set_3).difference(set_1)),
                len(set_1.intersection(set_3).intersection(set_2)),
            ),
            set_labels=tuple(labels),  # pyright: ignore[reportArgumentType]
            set_colors=tuple(colors),  # pyright: ignore[reportArgumentType]
            alpha=alpha,
        )

        # draw contour
        if contour:
            venn3_circles(
                subsets=(
                    len(set_1.difference(set_2, set_3)),
                    len(set_2.difference(set_1, set_3)),
                    len(set_1.intersection(set_2).difference(set_3)),
                    len(set_3.difference(set_1, set_2)),
                    len(set_1.intersection(set_3).difference(set_2)),
                    len(set_2.intersection(set_3).difference(set_1)),
                    len(set_1.intersection(set_3).intersection(set_2)),
                ),
                linewidth=linewidth,  # pyright: ignore[reportArgumentType]
            )

        # save file
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
            plt.title(title)
            plt.savefig(
                filename_prefix + ".png", dpi=300, transparent=True, bbox_inches="tight"
            )
            plt.savefig(
                filename_prefix + ".svg", dpi=300, transparent=True, bbox_inches="tight"
            )
        else:
            plt.title(title)

    return (fig, ax)


def plot_venn_diagram(
    data_1: List[Dict[str, Any]],
    data_2: List[Dict[str, Any]],
    data_3: Optional[List[Dict[str, Any]]] = None,
    by: Literal["peptide", "protein"] = "peptide",
    labels: List[str] = ["Set 1", "Set 2", "Set 3"],
    colors: List[str] = ["#4361EE", "#4CC9F0", "#F72585"],
    alpha: float = 0.6,
    contour: bool = False,
    linewidth: float = 0.5,
    title: str = "Venn Diagram",
    figsize: Tuple[float, float] = (10.0, 10.0),
    filename_prefix: Optional[str] = None,
) -> Tuple[Figure, Any]:
    r"""Plot the venn diagram for two or three sets of crosslink-spectrum-matches or crosslinks.

    Plot the venn diagram for two or three sets of crosslink-spectrum-matches or crosslinks. Overlaps
    are calculated by either looking at peptide sequence and crosslink position in the peptide using
    parameter by = "peptide" or by looking at protein crosslink position by using parameter by = "protein".

    Parameters
    ----------
    data_1 : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks.
    data_2 : list of dict of str, any
        A list of crosslink-spectrum-matches or crosslinks.
    data_3 : list of dict of str, any, or None, default = None
        Optionally, a third list of crosslink-spectrum-matches or crosslinks.
    by : str, one of "peptide" or "protein"
        If peptide or protein crosslink position should be used for determining if a crosslink-spectrum-match
        or crosslink is unique.
    labels : List[str], default = ["Set 1", "Set 2", "Set 3"]
        List of labels for the sets.
    colors : List[str], default = ["#4361EE", "#4CC9F0", "#F72585"]
        List of valid colors to use for the venn circles.
    alpha : float, default = 0.6
        Color opacity.
    contour : bool, default = False
        If a contour should be drawn around venn circles.
    linewidth : float, default = 0.5
        Linewidth of the contour.
    title : str, default = "Venn Diagram"
        Title of the venn diagram.
    figsize : tuple of float, float, default = (10.0, 10.0)
        Width, height in inches.
    filename_prefix: str, or None, default = None
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
    TypeError
        If parameter by is not one of 'peptide' or 'protein'.
    ValueError
        If one of the data parameters does not contain any crosslink-spectrum-matches or crosslinks.
    ValueError
        If attribute 'alpha_proteins', 'alpha_proteins_crosslink_positions', 'beta_proteins', or
        'beta_proteins_crosslink_positions' is not available for any of the data and parameter 'by'
        was set to 'protein'.

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS import plotting
    >>> a = parser.read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt", engine="MS Annika", crosslinker="DSS")
    >>> a = a["crosslink-spectrum-matches"]
    >>> b = parser.read("data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS")
    >>> b = b["crosslink-spectrum-matches"]
    >>> fig, ax = plotting.plot_venn_diagram(a, b, labels=["MS Annika", "MaxQuant"], colors=["orange", "blue"])

    >>> from pyXLMS import parser
    >>> from pyXLMS import plotting
    >>> a = parser.read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt", engine="MS Annika", crosslinker="DSS")
    >>> a = a["crosslink-spectrum-matches"]
    >>> b = parser.read("data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS")
    >>> b = b["crosslink-spectrum-matches"]
    >>> c = parser.read("data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv", engine="pLink", crosslinker="DSS")
    >>> c = c["crosslink-spectrum-matches"]
    >>> fig, ax = plotting.plot_venn_diagram(a, b, c, labels=["MS Annika", "MaxQuant", "pLink"], contour=True)
    """
    _ok = check_input(data_1, "data_1", list, dict)
    _ok = check_input(data_2, "data_2", list, dict)
    _ok = check_input(data_3, "data_3", list, dict) if data_3 is not None else True
    _ok = check_input(by, "by", str)
    _ok = check_input(labels, "labels", list, str)
    _ok = check_input(colors, "colors", list, str)
    _ok = check_input(alpha, "alpha", float)
    _ok = check_input(contour, "contour", bool)
    _ok = check_input(linewidth, "linewidth", float)
    _ok = check_input(title, "title", str)
    _ok = check_input(figsize, "figsize", tuple)
    _ok = (
        check_input(filename_prefix, "filename_prefix", str)
        if filename_prefix is not None
        else True
    )
    if by not in ["peptide", "protein"]:
        raise TypeError(
            "Parameter 'by' has to be one of 'peptide' or 'protein'! Option 'peptide' will group by peptide sequence and "
            "peptide crosslink position while option 'protein' will group by protein identifier and protein crosslink position."
        )
    if len(data_1) == 0:
        raise ValueError(
            "Can't plot venn diagram if no crosslink-spectrum-matches or crosslinks are given in data_1!"
        )
    if len(data_2) == 0:
        raise ValueError(
            "Can't plot venn diagram if no crosslink-spectrum-matches or crosslinks are given in data_2!"
        )
    if data_3 is not None and len(data_3) == 0:
        raise ValueError(
            "Can't plot 3-set venn diagram if no crosslink-spectrum-matches or crosslinks are given in data_3!"
        )
    if "data_type" not in data_1[0] or data_1[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data_1 has to be a list of crosslink or crosslink-spectrum-match!"
        )
    if "data_type" not in data_2[0] or data_2[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data_2 has to be a list of crosslink or crosslink-spectrum-match!"
        )
    if data_3 is not None:
        if "data_type" not in data_3[0] or data_3[0]["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data_3 has to be a list of crosslink or crosslink-spectrum-match, or None!"
            )

    set_1 = set()
    set_2 = set()
    set_3 = set()
    if by == "protein":
        available_keys_1 = get_available_keys(data_1)
        available_keys_2 = get_available_keys(data_2)
        if data_3 is None:
            alpha_proteins = (
                available_keys_1["alpha_proteins"]
                and available_keys_2["alpha_proteins"]
            )
            alpha_proteins_crosslink_positions = (
                available_keys_1["alpha_proteins_crosslink_positions"]
                and available_keys_2["alpha_proteins_crosslink_positions"]
            )
            beta_proteins = (
                available_keys_1["beta_proteins"] and available_keys_2["beta_proteins"]
            )
            beta_proteins_crosslink_positions = (
                available_keys_1["beta_proteins_crosslink_positions"]
                and available_keys_2["beta_proteins_crosslink_positions"]
            )
            if (
                not alpha_proteins
                or not alpha_proteins_crosslink_positions
                or not beta_proteins
                or not beta_proteins_crosslink_positions
            ):
                raise ValueError(
                    "Grouping by protein crosslink position is only available if all data have defined protein crosslink positions!\n"
                    "This error might be fixable with 'transform.reannotate_positions()'!"
                )
            for item in data_1:
                set_1.add(__get_key(item, by))
            for item in data_2:
                set_2.add(__get_key(item, by))
        else:
            available_keys_3 = get_available_keys(data_3)
            alpha_proteins = (
                available_keys_1["alpha_proteins"]
                and available_keys_2["alpha_proteins"]
                and available_keys_3["alpha_proteins"]
            )
            alpha_proteins_crosslink_positions = (
                available_keys_1["alpha_proteins_crosslink_positions"]
                and available_keys_2["alpha_proteins_crosslink_positions"]
                and available_keys_3["alpha_proteins_crosslink_positions"]
            )
            beta_proteins = (
                available_keys_1["beta_proteins"]
                and available_keys_2["beta_proteins"]
                and available_keys_3["beta_proteins"]
            )
            beta_proteins_crosslink_positions = (
                available_keys_1["beta_proteins_crosslink_positions"]
                and available_keys_2["beta_proteins_crosslink_positions"]
                and available_keys_3["beta_proteins_crosslink_positions"]
            )
            if (
                not alpha_proteins
                or not alpha_proteins_crosslink_positions
                or not beta_proteins
                or not beta_proteins_crosslink_positions
            ):
                raise ValueError(
                    "Grouping by protein crosslink position is only available if all data have defined protein crosslink positions!\n"
                    "This error might be fixable with 'transform.reannotate_positions()'!"
                )
            for item in data_1:
                set_1.add(__get_key(item, by))
            for item in data_2:
                set_2.add(__get_key(item, by))
            for item in data_3:
                set_3.add(__get_key(item, by))
    else:
        if data_3 is None:
            for item in data_1:
                set_1.add(__get_key(item, by))
            for item in data_2:
                set_2.add(__get_key(item, by))
        else:
            for item in data_1:
                set_1.add(__get_key(item, by))
            for item in data_2:
                set_2.add(__get_key(item, by))
            for item in data_3:
                set_3.add(__get_key(item, by))

    return venn(
        set_1,
        set_2,
        set_3 if data_3 is not None else None,
        labels=labels,
        colors=colors,
        alpha=alpha,
        contour=contour,
        linewidth=linewidth,
        title=title,
        figsize=figsize,
        filename_prefix=filename_prefix,
    )
