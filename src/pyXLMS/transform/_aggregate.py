#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._util import check_input_multi
from ..data._csm import create_crosslink_from_csm
from ..data._parser_result import create_parser_result
from ._util import get_available_keys
from ._util import check_available_keys
from ._util import assert_csms
from ._util import assert_xls
from ._util import assert_csms_or_xls

from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def __score_better(
    score: float, reference: float, function: Literal["higher_better", "lower_better"]
) -> bool:
    r"""Checks if the score is better than the provided reference score.

    Checks if the score is better than the provided reference score using the given scoring scheme.

    Parameters
    ----------
    score : float
        The score that should be compared.
    reference : float
        The reference score to compare to.
    function : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    bool
        If the given score is better than the reference score.
    """
    if function == "higher_better":
        return score > reference
    return score < reference


def __get_csm_key(csm: CrosslinkSpectrumMatch) -> str:
    r"""Get the unique key for a crosslink-spectrum-match.

    Parameters
    ----------
    csm : CrosslinkSpectrumMatch
        A pyXLMS crosslink-spectrum-match object.

    Returns
    -------
    str
        The unique key for the crosslink-spectrum-match.
    """
    return f"{csm['spectrum_file']}_{csm['scan_nr']}"


def __get_xl_key(
    xl: Crosslink | CrosslinkSpectrumMatch, by: Literal["peptide", "protein"]
) -> str:
    r"""Get the unique key for a crosslink.

    Parameters
    ----------
    xl : Crosslink, or CrosslinkSpectrumMatch
        A pyXLMS crosslink object.
        Technically, a crosslink-spectrum-match is also allowed for type support in some functions.
    by : str, one of "peptide" or "protein"
        If peptide or protein crosslink position should be used for determining if a crosslink is unique.

    Returns
    -------
    str
        The unique key for the crosslink.

    Notes
    -----
    This function should not be called directly, it is called from ``__unique_xls()``.
    """
    if by == "peptide":
        return (
            f"{xl['alpha_peptide']}_{xl['alpha_peptide_crosslink_position']}_{xl['alpha_decoy']}-"
            f"{xl['beta_peptide']}_{xl['beta_peptide_crosslink_position']}_{xl['beta_decoy']}"
        )
    _ok = check_available_keys(
        [
            "alpha_proteins",
            "alpha_proteins_crosslink_positions",
            "beta_proteins",
            "beta_proteins_crosslink_positions",
        ],
        assert_csms_or_xls([xl]),
    )
    prot_pos_a = (
        "-".join(
            sorted(
                [
                    f"{xl['alpha_proteins'][i]}_{xl['alpha_proteins_crosslink_positions'][i]}"
                    for i in range(len(xl["alpha_proteins"]))
                ]
            )
        )
        + f"_{xl['alpha_decoy']}"
    )
    prot_pos_b = (
        "-".join(
            sorted(
                [
                    f"{xl['beta_proteins'][i]}_{xl['beta_proteins_crosslink_positions'][i]}"
                    for i in range(len(xl["beta_proteins"]))
                ]
            )
        )
        + f"_{xl['beta_decoy']}"
    )
    return ":".join(sorted([prot_pos_a, prot_pos_b]))


def __unique_csms(
    csms: List[CrosslinkSpectrumMatch],
    has_scores: bool,
    score: Literal["higher_better", "lower_better"],
) -> List[CrosslinkSpectrumMatch]:
    r"""Filter for unique crosslink-spectrum-matches from a list on non-unique crosslink-spectrum-matches.

    Filters for unique crosslink-spectrum-matches from a list on non-unique crosslink-spectrum-matches. A crosslink-
    spectrum-match is considered unique if there is no other crosslink-spectrum-match from the same spectrum file and
    with the same scan number. If more than one crosslink-spectrum-match exists per spectrum file and scan number, the
    one with the better/best score is kept and the rest is filtered out. If crosslink-spectrum-matches without scores
    are provided, the first crosslink-spectrum-match in the list is kept instead.

    Parameters
    ----------
    csms : list of CrosslinkSpectrumMatch
        A list of pyXLMS crosslink-spectrum-match objects.
    has_scores : bool
        If the crosslink-spectrum-match objects contain scores.
    score : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    list of CrosslinkSpectrumMatch
        List of unique crosslink-spectrum-matches.

    Notes
    -----
    This function should not be called directly, it is called from ``unique()``.
    """
    unique_csms = dict()
    for csm in csms:
        key = __get_csm_key(csm)
        if key not in unique_csms:
            unique_csms[key] = csm
        elif has_scores and __score_better(
            csm["score"], unique_csms[key]["score"], score
        ):
            unique_csms[key] = csm
        else:
            # do nothing
            pass
    return list(unique_csms.values())


def __unique_xls(
    xls: List[Crosslink],
    by: Literal["peptide", "protein"],
    has_scores: bool,
    score: Literal["higher_better", "lower_better"],
) -> List[Crosslink]:
    r"""Filter for unique crosslinks from a list on non-unique crosslinks.

    Filters for unique crosslinks from a list on non-unique crosslinks. A crosslink is considered unique if there is no
    other crosslink with the same peptide sequence and crosslink position if ``by = "peptide"``, otherwise it is considered
    unique if there are no other crosslinks with the same protein crosslink position (residue pair). If more than one
    crosslink exists per peptide sequence/residue pair, the one with the better/best score is kept and the rest is filtered
    out. If crosslinks without scores are provided, the first crosslink in the list is kept instead.

    Parameters
    ----------
    xls : list of Crosslink
        A list of pyXLMS crosslink objects.
    by : str, one of "peptide" or "protein"
        If peptide or protein crosslink position should be used for determining if a crosslink is unique.
    has_scores : bool
        If the crosslink objects contain scores.
    score : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    list of Crosslink
        List of unique crosslinks.

    Notes
    -----
    This function should not be called directly, it is called from ``unique()``.
    """
    unique_xls = dict()
    for xl in xls:
        key = __get_xl_key(xl, by)
        if key not in unique_xls:
            unique_xls[key] = xl
        elif has_scores and __score_better(
            xl["score"], unique_xls[key]["score"], score
        ):
            unique_xls[key] = xl
        else:
            # do nothing
            pass
    return list(unique_xls.values())


def unique(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult,
    by: Literal["peptide", "protein"] = "peptide",
    score: Literal["higher_better", "lower_better"] = "higher_better",
) -> List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult:
    r"""Filter for unique crosslinks or crosslink-spectrum-matches.

    Filters for unique crosslinks from a list on non-unique crosslinks. A crosslink is considered unique if there is no
    other crosslink with the same peptide sequence and crosslink position if ``by = "peptide"``, otherwise it is considered
    unique if there are no other crosslinks with the same protein crosslink position (residue pair). If more than one
    crosslink exists per peptide sequence/residue pair, the one with the better/best score is kept and the rest is filtered
    out. If crosslinks without scores are provided, the first crosslink in the list is kept instead.

    *or*

    Filters for unique crosslink-spectrum-matches from a list on non-unique crosslink-spectrum-matches. A crosslink-
    spectrum-match is considered unique if there is no other crosslink-spectrum-match from the same spectrum file and
    with the same scan number. If more than one crosslink-spectrum-match exists per spectrum file and scan number, the
    one with the better/best score is kept and the rest is filtered out. If crosslink-spectrum-matches without scores
    are provided, the first crosslink-spectrum-match in the list is kept instead.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        A list of crosslink-spectrum-matches or crosslinks to filter, or a parser_result.
    by : str, one of "peptide" or "protein", default = "peptide"
        If peptide or protein crosslink position should be used for determining if a crosslink is unique.
        Only affects filtering for unique crosslinks and not crosslink-spectrum-matches. If protein crosslink
        position is not available for all crosslinks a ``ValueError`` will be raised. Make sure that all
        crosslinks have the ``_proteins`` and ``_proteins_crosslink_positions`` fields set. If this is not
        already done by the parser, this can be achieved with ``transform.reannotate_positions()``.
    score : str, one of "higher_better" or "lower_better", default = "higher_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        If a list of crosslink-spectrum-matches or crosslinks was provided, a list of unique
        crosslink-spectrum-matches or crosslinks is returned. If a parser_result was provided,
        a parser_result with unique crosslink-spectrum-matches and/or unique crosslinks will
        be returned.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If parameter by is not one of 'peptide' or 'protein'.
    TypeError
        If parameter score is not one of 'higher_better' or 'lower_better'.
    ValueError
        If parameter by is set to 'protein' but protein crosslink positions are not available.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import unique
    >>> pr = read(
    ...     ["data/_test/aggregate/csms.txt", "data/_test/aggregate/xls.txt"],
    ...     engine="custom",
    ...     crosslinker="DSS",
    ... )
    >>> len(pr["crosslink-spectrum-matches"])
    10
    >>> len(pr["crosslinks"])
    10
    >>> unique_peptide = unique(pr, by="peptide")
    >>> len(unique_peptide["crosslink-spectrum-matches"])
    5
    >>> len(unique_peptide["crosslinks"])
    3

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import unique
    >>> pr = read(
    ...     ["data/_test/aggregate/csms.txt", "data/_test/aggregate/xls.txt"],
    ...     engine="custom",
    ...     crosslinker="DSS",
    ... )
    >>> len(pr["crosslink-spectrum-matches"])
    10
    >>> len(pr["crosslinks"])
    10
    >>> unique_protein = unique(pr, by="protein")
    >>> len(unique_protein["crosslink-spectrum-matches"])
    5
    >>> len(unique_protein["crosslinks"])
    2
    """
    _ok = check_input_multi(data, "data", [ParserResult, list])
    _ok = check_input(by, "by", str)
    _ok = check_input(score, "score", str)
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
    if isinstance(data, list):
        if len(data) == 0:
            return data
        data = assert_csms_or_xls(data)
        available_keys = get_available_keys(data)
        # crosslink and by protein
        if isinstance(data[0], Crosslink) and by == "protein":
            _ok = check_available_keys(
                [
                    "alpha_proteins",
                    "alpha_proteins_crosslink_positions",
                    "beta_proteins",
                    "beta_proteins_crosslink_positions",
                ],
                data,
            )
            return __unique_xls(assert_xls(data), by, available_keys["score"], score)
        # crosslink but not by protein
        if isinstance(data[0], Crosslink):
            return __unique_xls(assert_xls(data), by, available_keys["score"], score)
        # csm
        return __unique_csms(assert_csms(data), available_keys["score"], score)
    new_csms = (
        assert_csms(unique(data["crosslink-spectrum-matches"], by, score))
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        assert_xls(unique(data["crosslinks"], by, score))
        if data["crosslinks"] is not None
        else None
    )
    return create_parser_result(
        search_engine=data["search_engine"],
        csms=new_csms,
        crosslinks=new_xls,
    )


def aggregate(
    csms: List[CrosslinkSpectrumMatch],
    by: Literal["peptide", "protein"] = "peptide",
    score: Literal["higher_better", "lower_better"] = "higher_better",
) -> List[Crosslink]:
    r"""Aggregate crosslink-spectrum-matches to crosslinks.

    Aggregates a list of crosslink-spectrum-matches to unique crosslinks. A crosslink is considered unique if there is no
    other crosslink with the same peptide sequence and crosslink position if ``by = "peptide"``, otherwise it is considered
    unique if there are no other crosslinks with the same protein crosslink position (residue pair). If more than one
    crosslink exists per peptide sequence/residue pair, the one with the better/best score is kept and the rest is filtered
    out. If crosslink-spectrum-matches without scores are provided, the crosslink of the first corresponding crosslink-spectrum
    -match in the list is kept instead.

    Parameters
    ----------
    csms : list of CrosslinkSpectrumMatch
        A list of crosslink-spectrum-matches.
    by : str, one of "peptide" or "protein", default = "peptide"
        If peptide or protein crosslink position should be used for determining if a crosslink is unique.
        If protein crosslink position is not available for all crosslink-spectrum-matches a ``ValueError``
        will be raised. Make sure that all crosslink-spectrum-matches have the ``_proteins`` and
        ``_proteins_crosslink_positions`` fields set. If this is not already done by the parser, this can
        be achieved with ``transform.reannotate_positions()``.
    score : str, one of "higher_better" or "lower_better", default = "higher_better"
        If a higher score is considered better, or a lower score is considered better.

    Returns
    -------
    list of Crosslink
        A list of aggregated, unique crosslinks.

    Warnings
    --------
    Aggregation will not conserve false discovery rate (FDR)! Aggregating crosslink-spectrum-matches that are
    validated for 1% FDR will not result in crosslinks validated for 1% FDR! Aggregated crosslinks should be
    validated with either external tools or with the built-in ``transform.validate()``!

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If parameter by is not one of 'peptide' or 'protein'.
    TypeError
        If parameter score is not one of 'higher_better' or 'lower_better'.
    ValueError
        If parameter by is set to 'protein' but protein crosslink positions are not available.

    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import aggregate
    >>> pr = read("data/_test/aggregate/csms.txt", engine="custom", crosslinker="DSS")
    >>> len(pr["crosslink-spectrum-matches"])
    10
    >>> aggregate_peptide = aggregate(pr["crosslink-spectrum-matches"], by="peptide")
    >>> len(aggregate_peptide)
    3

    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import aggregate
    >>> pr = read("data/_test/aggregate/csms.txt", engine="custom", crosslinker="DSS")
    >>> len(pr["crosslink-spectrum-matches"])
    10
    >>> aggregate_protein = aggregate(pr["crosslink-spectrum-matches"], by="protein")
    >>> len(aggregate_protein)
    2
    """
    _ok = check_input(csms, "csms", list, CrosslinkSpectrumMatch)
    _ok = check_input(by, "by", str)
    _ok = check_input(score, "score", str)
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
    if len(csms) == 0:
        return []
    if by == "protein":
        _ok = check_available_keys(
            [
                "alpha_proteins",
                "alpha_proteins_crosslink_positions",
                "beta_proteins",
                "beta_proteins_crosslink_positions",
            ],
            csms,
        )
    xls = [create_crosslink_from_csm(csm) for csm in csms]
    return assert_xls(unique(xls, by, score))
