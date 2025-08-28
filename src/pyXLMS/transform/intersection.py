#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import warnings

from ..data import check_input
from .aggregate import __score_better
from .aggregate import __get_xl_key
from .aggregate import unique
from .util import get_available_keys

from typing import Dict
from typing import List
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __get_csm_key(csm: Dict[str, Any]) -> str:
    r"""Get the unique intersection key for a crosslink-spectrum-match.

    Parameters
    ----------
    csm : dict of str, any
        A pyXLMS crosslink-spectrum-match object.

    Returns
    -------
    str
        The unique intersection key for the crosslink-spectrum-match.
    """
    k1 = csm["alpha_peptide"]
    k2 = csm["alpha_peptide_crosslink_position"]
    k3 = csm["beta_peptide"]
    k4 = csm["beta_peptide_crosslink_position"]
    k5 = csm["spectrum_file"]
    k6 = csm["scan_nr"]
    return f"{k1}{k2}-{k3}{k4}_{k5}_{k6}"


def intersection(
    data_a: List[Dict[str, Any]],
    data_b: List[Dict[str, Any]],
    use: Literal["better_score", "data_a", "data_b"] = "better_score",
    by: Literal["peptide", "protein"] = "peptide",
    score: Literal["higher_better", "lower_better"] = "higher_better",
    verbose: Literal[0, 1, 2] = 1,
) -> List[Dict[str, Any]]:
    r"""Get the intersection of two lists of crosslinks.

    Returns the intersection of two lists of crosslinks (or crosslink-spectrum-matches). Crosslink intersection is calculated
    by creating unique keys based on peptide sequence and peptide crosslink position or based on protein crosslink position.
    For crosslink-spectrum-match intersection the intersection is calculated by creating unique keys based on peptide sequence
    and peptide crosslink position as well as the corresponding spectrum file and scan number. For any unique key in the intersection,
    either the best (by score) crosslink or crosslink-spectrum-match is returned, or the one from the first or second list, based
    on user preference.

    Parameters
    ----------
    data_a : list of dict of str, any
        List of crosslinks (or crosslink-spectrum-matches).
    data_b : list of dict of str, any
        List of crosslinks (or crosslink-spectrum-matches) to intersect with. Note that the data types for ``data_a`` and
        ``data_b`` have to be the same.
    use : str, one of "better_score", "data_a", or "data_b", default = "better_score"
        Which element to use for the returned intersection. Option "better_score" will return the element (crosslink or
        crosslink-spectrum-match) with the higher score, option "data_a" will return the element from ``data_a``, and
        option "data_b" will return the element from ``data_b``. Please note that attribute ``score`` needs to be
        available for all elements if "better_score" is selected.
    by : str, one of "peptide" or "protein", default = "peptide"
        If peptide or protein crosslink position should be used for determining if a crosslink is unique.
        Only affects filtering for unique crosslinks and not crosslink-spectrum-matches. If protein crosslink
        position is not available for all crosslinks a ``ValueError`` will be raised. Make sure that all
        crosslinks have the ``_proteins`` and ``_proteins_crosslink_positions`` fields set. If this is not
        already done by the parser, this can be achieved with ``transform.reannotate_positions()``.
    score: str, one of "higher_better" or "lower_better", default = "higher_better"
        If a higher score is considered better, or a lower score is considered better.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    list of dict of str, any
        The list of crosslinks or crosslink-spectrum-matches in the intersection.

    Raises
    ------
    TypeError
        If a wrong data type is provided for 'data_a' or 'data_b'.
    TypeError
        If 'data_a' and 'data_b' are not of the same data type.
    TypeError
        If parameter use is not one of 'better_score', 'data_a' or 'data_b'.
    TypeError
        If parameter by is not one of 'peptide' or 'protein'.
    TypeError
        If parameter score is not one of 'higher_better' or 'lower_better'.
    TypeError
        If parameter verbose was not set correctly.
    ValueError
        If parameter use is set to 'better_score' but scores are not available.
    ValueError
        If parameter by is set to 'protein' but protein crosslink positions are not available.
    RuntimeError
        If crosslink-spectrum-matches are provided as data but verbose level is set to 2.

    Warns
    -----
    RuntimeWarning
        If crosslink-spectrum-matches are provided as data but verbose level is set to 1.

    Notes
    -----
    While technically the intersection of crosslink-spectrum-matches is supported, please be aware that this mostly only makes
    sense when searching the same mass spectrometry file several times, as the crosslink-spectrum-match intersection takes into
    account the spectrum files and scan numbers, e.g. two crosslinks with the same peptide sequences but different scan numbers will
    not be considered the same and therefore will not be in the intersection. If you want to get the intersection of crosslink-spectrum-matches
    using their crosslinked peptide sequences or protein positions, you can aggregate them first with ``transform.aggregate()`` or create
    crosslinks with ``data.create_crosslink_from_csm()``. Please also note that when intersecting all duplicates are removed and
    aggregated the same way as in ``transform.unique()``. To intersect more than two lists please repeatedly call ``intersection()``
    on itself.

    Examples
    --------
    >>> from pyXLMS.pipelines import pipeline
    >>> from pyXLMS.transform import aggregate
    >>> from pyXLMS.transform import intersection
    >>> msannika = pipeline(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> msannika = aggregate(msannika["crosslink-spectrum-matches"])
    >>> len(msannika)
    235
    >>> maxquant = pipeline(
    ...     "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    ... )
    >>> maxquant = aggregate(maxquant["crosslink-spectrum-matches"])
    >>> len(maxquant)
    226
    >>> crosslinks_intersection = intersection(msannika, maxquant)
    >>> len(crosslinks_intersection)
    206

    >>> from pyXLMS.pipelines import pipeline
    >>> from pyXLMS.transform import aggregate
    >>> from pyXLMS.transform import intersection
    >>> msannika = pipeline(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> msannika = aggregate(msannika["crosslink-spectrum-matches"])
    >>> len(msannika)
    235
    >>> maxquant = pipeline(
    ...     "data/maxquant/run1/crosslinkMsms.txt", engine="MaxQuant", crosslinker="DSS"
    ... )
    >>> maxquant = aggregate(maxquant["crosslink-spectrum-matches"])
    >>> len(maxquant)
    226
    >>> crosslinks_intersection = intersection(msannika, maxquant)
    >>> len(crosslinks_intersection)
    206
    >>> plink = pipeline(
    ...     "data/plink2/Cas9_plus10_2024.06.20.filtered_cross-linked_spectra.csv",
    ...     engine="pLink",
    ...     crosslinker="DSS",
    ... )
    >>> plink = aggregate(plink["crosslink-spectrum-matches"])
    >>> len(plink)
    252
    >>> crosslinks_intersection = intersection(crosslinks_intersection, plink)
    >>> len(crosslinks_intersection)
    203
    """
    _ok = check_input(data_a, "data_a", list, dict)
    _ok = check_input(data_b, "data_b", list, dict)
    _ok = check_input(use, "use", str)
    _ok = check_input(by, "by", str)
    _ok = check_input(score, "score", str)
    if use not in ["better_score", "data_a", "data_b"]:
        raise TypeError(
            "Parameter 'use' has to be one 'better_score', 'data_a', or 'data_b'! Option 'better_score' will return the intersection "
            "with the better score (if scores are available). Option 'data_a' will return the intersection using elements of 'data_a'. "
            "Option 'data_b' will return the intersection using elements of 'data_b'."
        )
    if by not in ["peptide", "protein"]:
        raise TypeError(
            "Parameter 'by' has to be one of 'peptide' or 'protein'! Option 'peptide' will group by peptide sequence and "
            "peptide crosslink position while option 'protein' will group by protein identifier and protein crosslink position."
        )
    if score not in ["higher_better", "lower_better"]:
        raise TypeError(
            "Parameter 'score' has to be one of 'higher_better' or 'lower_better'! If two identical crosslinks or crosslink-spectrum"
            "-matches are found, the one with the higher score is kept if 'higher_better' is selected, and vice versa."
        )
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")
    if len(data_a) == 0 or len(data_b) == 0:
        return []
    if (
        "data_type" not in data_a[0]
        or data_a[0]["data_type"]
        not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]
        or "data_type" not in data_b[0]
        or data_b[0]["data_type"]
        not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]
    ):
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
        )
    if data_a[0]["data_type"] != data_b[0]["data_type"]:
        raise TypeError(
            "Parameters 'data_a' and 'data_b' have to be of the same data type!"
        )
    available_keys_a = get_available_keys(data_a)
    available_keys_b = get_available_keys(data_b)
    if use == "better_score" and (
        not available_keys_a["score"] or not available_keys_b["score"]
    ):
        raise ValueError(
            "Can't intersect based on score because not all data have associated scores!"
        )
    unique_a = unique(data_a, by=by, score=score)
    unique_b = unique(data_b, by=by, score=score)
    if not isinstance(unique_a, list) or not isinstance(unique_b, list):
        raise RuntimeError(
            "Something went wrong while getting unique data.\n"
            f"Expected data type: list. Got: {type(unique_a)} and {type(unique_b)}."
        )
    if unique_a[0]["data_type"] == "crosslink":
        crosslinks_a = dict()
        for xl in unique_a:
            crosslinks_a[__get_xl_key(xl, by=by)] = xl
        crosslinks_b = dict()
        for xl in unique_b:
            crosslinks_b[__get_xl_key(xl, by=by)] = xl
        keys_intersection = set(crosslinks_a.keys()).intersection(
            set(crosslinks_b.keys())
        )
        crosslinks_intersection = list()
        for key in keys_intersection:
            if use == "data_a":
                crosslinks_intersection.append(crosslinks_a[key])
            elif use == "data_b":
                crosslinks_intersection.append(crosslinks_b[key])
            else:
                if __score_better(
                    score=crosslinks_a[key]["score"],
                    reference=crosslinks_b[key]["score"],
                    function=score,
                ):
                    crosslinks_intersection.append(crosslinks_a[key])
                else:
                    crosslinks_intersection.append(crosslinks_b[key])
        return crosslinks_intersection
    if verbose == 1:
        warnings.warn(
            RuntimeWarning(
                "Creating intersection of crosslink-spectrum-matches. Be sure that this makes sense for your data!"
            )
        )
    elif verbose == 2:
        raise RuntimeError(
            "Can't create intersection of crosslink-spectrum-matches for verbose level 2!"
        )
    csms_a = dict()
    for csm in unique_a:
        csms_a[__get_csm_key(csm)] = csm
    csms_b = dict()
    for csm in unique_b:
        csms_b[__get_csm_key(csm)] = csm
    keys_intersection = set(csms_a.keys()).intersection(set(csms_b.keys()))
    csms_intersection = list()
    for key in keys_intersection:
        if use == "data_a":
            csms_intersection.append(csms_a[key])
        elif use == "data_b":
            csms_intersection.append(csms_b[key])
        else:
            if __score_better(
                score=csms_a[key]["score"],
                reference=csms_b[key]["score"],
                function=score,
            ):
                csms_intersection.append(csms_a[key])
            else:
                csms_intersection.append(csms_b[key])
    return csms_intersection
