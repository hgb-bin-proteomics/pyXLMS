#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import time
import requests
import warnings
from tqdm import tqdm

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._util import check_input_multi
from ..data._parser_result import create_parser_result
from ._filter import filter_crosslink_type
from ._filter import filter_protein_distribution
from ._util import get_available_keys
from ._util import assert_csms
from ._util import assert_xls
from ._util import assert_csms_or_xls

from typing import List
from typing import Dict
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

STRING_STABLE_URL = "https://version-12-0.string-db.org/api"
CALLER_IDENTITY = "https://github.com/hgb-bin-proteomics/pyXLMS"
STRING_ORGANISMS = {
    "Homo sapiens": 9606,
    "Mus musculus": 10090,
    "Arabidopsis thaliana": 3702,
    "Saccharomyces cerevisiae": 4932,
    "Drosophila melanogaster": 7227,
    "Danio rerio": 7955,
    "Caenorhabditis elegans": 6239,
    "Escherichia coli str. K-12 substr. MG1655": 511145,
    "Pseudomonas aeruginosa PAO1": 208964,
}
# from https://string-db.org/help/getting_started/
STRING_SCORES = {
    "low confidence": 0.15,
    "medium confidence": 0.4,
    "high confidence": 0.7,
    "highest confidence": 0.9,
}


def __float_or_none(value: Any) -> float | None:
    r"""Converts a value to float if possible and returns None if otherwise.

    Parameters
    ----------
    value : any
        The value to convert to float.

    Returns
    -------
    float, or None
        Returns a float if the value can be cast to float, otherwise None.
    """
    try:
        return float(value)
    except Exception as _e:
        pass
    return None


def get_string_ids(
    proteins: List[str], organism: str | int, verbose: Literal[0, 1, 2] = 1
) -> Dict[str, str | None]:
    r"""Map proteins to STRING IDs.

    Calls the STRING API to resolve common protein or gene names, synonyms, or UniProt identifiers
    and accession numbers to map them to the identifiers used internally by STRING.
    STRING is accessible via `string-db.org <https://string-db.org>`_.

    Parameters
    ----------
    proteins : list of str
        A list of protein/gene accessions, names, or symbols.
    organism : str, or int
        Organism name (e.g. Homo sapiens) or taxon identifier (e.g. 9606).
        Taxon identifiers are preferred. See also
        `string-db.org/cgi/organisms <https://string-db.org/cgi/organisms>`_.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    dict of str, str or None
        Returns a dictionary that maps the input proteins to their STRING IDs. If
        a protein could not be resolved its STRING ID is None.

    Raises
    ------
    TypeError
        If parameter verbose was not set correctly.
    KeyError
        If the organism could not be resolved to a taxon identifier.
    RuntimeError
        If ``verbose = 2`` and the API request failed.

    Warns
    -----
    RuntimeWarning
        If ``verbose = 1`` and the API request failed.

    Examples
    --------
    >>> from pyXLMS.transform import get_string_ids
    >>> get_string_ids(["p53", "BRCA1", "cdk2", "Q99835"], organism=9606)
    {'p53': '9606.ENSP00000269305', 'BRCA1': '9606.ENSP00000418960', 'cdk2': '9606.ENSP00000266970', 'Q99835': '9606.ENSP00000249373'}
    """
    _ok = check_input(proteins, "proteins", list, str)
    _ok = check_input_multi(organism, "organism", [str, int])
    if isinstance(organism, str):
        if organism not in STRING_ORGANISMS:
            raise KeyError(
                f"Could not resolve organism {organism}, please specify taxon identifier manually!"
            )
        organism = STRING_ORGANISMS[organism]
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")
    output_proteins: Dict[str, str | None] = dict()
    params = {
        "identifiers": "\r".join(proteins),
        "species": organism,
        "echo_query": 1,
        "caller_identity": CALLER_IDENTITY,
    }
    request_url = f"{STRING_STABLE_URL}/json/get_string_ids"
    response: requests.models.Response | None = None
    try:
        response = requests.post(request_url, data=params)
    except Exception as e:
        response = None
        if verbose == 1:
            warnings.warn(
                RuntimeWarning(f"Request to STRING API failed with error {e}!")
            )
        if verbose == 2:
            raise
    if response is None:
        return output_proteins
    # wait one second after request to delay subsequent requests - be polite
    time.sleep(1)
    if not response.ok:
        if verbose == 1:
            warnings.warn(RuntimeWarning(f"{response.text}"))
        if verbose == 2:
            raise RuntimeError(f"{response.text}")
        return output_proteins
    response_json = response.json()
    response_proteins: Dict[str, str] = dict()
    for item in response_json:
        if "queryItem" in item and "stringId" in item:
            response_proteins[str(item["queryItem"]).strip()] = str(
                item["stringId"]
            ).strip()
    for protein in proteins:
        if protein in response_proteins:
            output_proteins[protein] = response_proteins[protein]
        else:
            output_proteins[protein] = None
    return output_proteins


def get_string_network(
    string_ids: List[str], organism: str | int, verbose: Literal[0, 1, 2] = 1
) -> Dict[str, Dict[str, str | float | None]]:
    r"""Retrieve a STRING interaction network.

    Retrieves the STRING interaction network via the STRING API for the given STRING IDs
    with interactions including the combined score and all the channel specific scores.
    STRING is accessible via `string-db.org <https://string-db.org>`_.

    Parameters
    ----------
    string_ids : list of str
        A list of STRING IDs for which the network should be created.
    organism : str, or int
        Organism name (e.g. Homo sapiens) or taxon identifier (e.g. 9606).
        Taxon identifiers are preferred. See also
        `string-db.org/cgi/organisms <https://string-db.org/cgi/organisms>`_.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    dict of str, dict of str, str or float or None
        Returns a dictionary where every key is an interaction composed of the two
        STRING IDs sorted alphabetically and joined via an underscore. The values are
        dictionaries with keys ``A`` for STRING ID A, ``B`` for STRING ID B, and
        ``score``, ``nscore``, ``fscore``, ``pscore``, ``ascore``, ``escore``,
        ``dscore``, ``tscore`` for the various STRING scores. Scores may be None.

    Raises
    ------
    TypeError
        If parameter verbose was not set correctly.
    KeyError
        If the organism could not be resolved to a taxon identifier.
    KeyError
        If ``verbose = 2`` and more than one interaction per protein pair was found.
    RuntimeError
        If ``verbose = 2`` and the API request failed.

    Warns
    -----
    RuntimeWarning
        If ``verbose = 1`` and the API request failed.
    RuntimeWarning
        If ``verbose = 1`` and more than one interaction per protein pair was found.

    Notes
    -----
    STRING limits the number of nodes via the API two 2000 - exceeding that limit will return an empty
    network/no annotation or raise an error (if ``verbose = 2``).

    Examples
    --------
    >>> from pyXLMS.transform import get_string_network
    >>> get_string_network(
    ...     ["9606.ENSP00000269305", "9606.ENSP00000418960"], organism=9606
    ... )
    {'9606.ENSP00000269305_9606.ENSP00000418960': {'A': '9606.ENSP00000269305', 'B': '9606.ENSP00000418960', 'score': 0.999, 'nscore': 0.0, 'fscore': 0.0, 'pscore': 0.0, 'ascore': 0.067, 'escore': 0.895, 'dscore': 0.5, 'tscore': 0.999}}
    """
    _ok = check_input(string_ids, "string_ids", list, str)
    _ok = check_input_multi(organism, "organism", [str, int])
    if isinstance(organism, str):
        if organism not in STRING_ORGANISMS:
            raise KeyError(
                f"Could not resolve organism {organism}, please specify taxon identifier manually!"
            )
        organism = STRING_ORGANISMS[organism]
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")
    network: Dict[str, Dict[str, str | float | None]] = dict()
    params = {
        "identifiers": "\r".join(string_ids),
        "species": organism,
        "required_score": 0,
        "add_nodes": 0,
        "show_query_node_labels": 0,
        "caller_identity": CALLER_IDENTITY,
    }
    request_url = f"{STRING_STABLE_URL}/json/network"
    response: requests.models.Response | None = None
    try:
        response = requests.post(request_url, data=params)
    except Exception as e:
        response = None
        if verbose == 1:
            warnings.warn(
                RuntimeWarning(f"Request to STRING API failed with error {e}!")
            )
        if verbose == 2:
            raise
    if response is None:
        return network
    # wait one second after request to delay subsequent requests - be polite
    time.sleep(1)
    if not response.ok:
        if verbose == 1:
            warnings.warn(RuntimeWarning(f"{response.text}"))
        if verbose == 2:
            raise RuntimeError(f"{response.text}")
        return network
    response_json = response.json()
    for item in response_json:
        a = str(item["stringId_A"]).strip()
        b = str(item["stringId_B"]).strip()
        key = "_".join(sorted([a, b]))
        parsed_item: Dict[str, str | float | None] = dict()
        parsed_item["A"] = a
        parsed_item["B"] = b
        # combined score
        parsed_item["score"] = __float_or_none(item["score"])
        # gene neighborhood score
        parsed_item["nscore"] = __float_or_none(item["nscore"])
        # gene fusion score
        parsed_item["fscore"] = __float_or_none(item["fscore"])
        # phylogenetic profile score
        parsed_item["pscore"] = __float_or_none(item["pscore"])
        # coexpression score
        parsed_item["ascore"] = __float_or_none(item["ascore"])
        # experimental score
        parsed_item["escore"] = __float_or_none(item["escore"])
        # database score
        parsed_item["dscore"] = __float_or_none(item["dscore"])
        # textmining score
        parsed_item["tscore"] = __float_or_none(item["tscore"])
        if key not in network:
            network[key] = parsed_item
        else:
            if network[key]["score"] is None:
                if parsed_item["score"] is not None:
                    network[key] = parsed_item
                else:
                    # do nothing
                    pass
            else:
                if parsed_item["score"] is None:
                    # do nothing
                    pass
                else:
                    if parsed_item["score"] > network[key]["score"]:  # pyright: ignore[reportOperatorIssue] # ty: ignore[unsupported-operator]
                        network[key] = parsed_item
                    else:
                        # do nothing
                        pass
            if verbose == 1:
                warnings.warn(
                    RuntimeWarning(
                        f"Found more than one interaction for {key}. Using highest scoring one!"
                    )
                )
            if verbose == 2:
                raise KeyError(f"Found more than one interaction for {key}!")
    return network


def annotate_string_scores(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult,
    organism: str | int,
    verbose: Literal[0, 1, 2] = 1,
) -> List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult:
    r"""Annotates STRING interactions and scores for inter-links.

    Annotates STRING interactions and STRING scores for inter-links based on their associated proteins.
    Takes a list of crosslink-spectrum-matches or crosslinks, or a parser_result as input.
    STRING is accessible via `string-db.org <https://string-db.org>`_.

    Parameters
    ----------
    data : list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        A list of crosslink-spectrum-matches or crosslinks to annotate, or a parser_result.
    organism : str, or int
        Organism name (e.g. Homo sapiens) or taxon identifier (e.g. 9606).
        Taxon identifiers are preferred. See also
        `string-db.org/cgi/organisms <https://string-db.org/cgi/organisms>`_.
    verbose : 0, 1, or 2, default = 1
        - 0: All warnings are ignored.
        - 1: Warnings are printed to stdout.
        - 2: Warnings are treated as errors.

    Returns
    -------
    list of CrosslinkSpectrumMatch, list of Crosslink, or ParserResult
        If a list of crosslink-spectrum-matches or crosslinks was provided, a list of annotated
        crosslink-spectrum-matches or crosslinks is returned. If a parser_result was provided,
        an annotated parser_result will be returned. Please note that only inter-links are
        annotated. Annotated interactions and scores are available via ``additional_information``
        using keys ``pyXLMS_annotated_STRING_interactions`` and ``pyXLMS_annotated_STRING_score``.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If parameter verbose was not set correctly.
    ValueError
        If the input data does not contain inter-links.
    KeyError
        If the organism could not be resolved to a taxon identifier.
    RuntimeError
        If ``verbose = 2`` and not all of the provided data does have associated proteins.
    RuntimeError
        If ``verbose = 2`` and data with more than 2000 proteins/STRING IDs is provided.
    RuntimeError
        If ``verbose = 2`` and the API request failed.

    Warns
    -----
    RuntimeWarning
        If ``verbose = 1`` and not all of the provided data does have associated proteins.
    RuntimeWarning
        If ``verbose = 1`` and data with more than 2000 proteins/STRING IDs is provided.
    RuntimeWarning
        If ``verbose = 1`` and the API request failed.

    Notes
    -----
    STRING limits the number of nodes via the API two 2000 - exceeding that limit will return an empty
    network/no annotation or raise an error (if ``verbose = 2``).

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS.transform import filter_crosslink_type
    >>> from pyXLMS.transform import annotate_string_scores
    >>> pr = parser.read_custom("data/ms_annika/Nucleus_Rep1_Crosslinks.parquet")
    >>> xls = pr["crosslinks"]
    >>> xls = annotate_string_scores(xls, organism="Homo sapiens")
    >>> inter = filter_crosslink_type(xls)["Inter"]
    >>> example = inter[4]  # example link with STRING score
    >>> example["additional_information"]["pyXLMS_annotated_STRING_interactions"]
    [{'A': '9606.ENSP00000441875', 'B': '9606.ENSP00000479488', 'score': 0.999, 'nscore': 0.0, 'fscore': 0.0, 'pscore': 0.068, 'ascore': 0.923, 'escore': 0.973, 'dscore': 0.9, 'tscore': 0.988}]
    >>> example["additional_information"]["pyXLMS_annotated_STRING_score"]
    0.999
    """
    _ok = check_input_multi(data, "data", [list, ParserResult])
    _ok = check_input_multi(organism, "organism", [str, int])
    if isinstance(organism, str):
        if organism not in STRING_ORGANISMS:
            raise KeyError(
                f"Could not resolve organism {organism}, please specify taxon identifier manually!"
            )
        organism = STRING_ORGANISMS[organism]
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")
    if isinstance(data, list):
        if len(data) == 0:
            return data
        data = assert_csms_or_xls(data)
        available_keys = get_available_keys(data)
        if not available_keys["alpha_proteins"] or not available_keys["beta_proteins"]:
            if verbose == 1:
                warnings.warn(
                    RuntimeWarning(
                        "Some of your crosslink-spectrum-matches/crosslinks do not have associated proteins. Their STRING scores will be nan!"
                    )
                )
            if verbose == 2:
                raise RuntimeError(
                    "Some of your crosslink-spectrum-matches/crosslinks do not have associated proteins!"
                )
        # annotate STRING scores
        # this if clause is technically not needed anymore but kept for legacy
        if (
            data[0]["data_type"] == "crosslink"
            or data[0]["data_type"] == "crosslink-spectrum-match"
        ):
            inter = filter_crosslink_type(data)["Inter"]
            if len(inter) == 0:
                raise ValueError(
                    "Can't annotate STRING scores for input data because it does not contain inter-links!"
                )
            proteins = list(filter_protein_distribution(inter).keys())
            proteins_to_string_ids = get_string_ids(proteins, organism, verbose)
            string_ids: List[str] = list()
            for k, v in proteins_to_string_ids.items():
                if v is not None:
                    string_ids.append(v)
            print(
                f"Mapped {len(string_ids)} of {len(proteins)} proteins ({len(string_ids) / len(proteins) * 100}%) to STRING IDs."
            )
            # this is a hard limit as of 2026-04
            if len(string_ids) >= 2000:
                if verbose == 1:
                    warnings.warn(
                        RuntimeWarning(
                            f"More than 2000 proteins/STRING IDs specified: {len(string_ids)}. Please reduce the number of proteins for a successful request!"
                        )
                    )
                if verbose == 2:
                    raise RuntimeError(
                        f"More than 2000 proteins/STRING IDs specified: {len(string_ids)}. Please reduce the number of proteins for a successful request!"
                    )
                return data
            network = get_string_network(string_ids, organism, verbose)
            for item in tqdm(
                inter,
                total=len(inter),
                desc="Annotating STRING scores for inter-links...",
            ):
                string_items: List[Dict[str, str | float | None]] = list()
                string_scores: List[float] = list()
                if (
                    item["alpha_proteins"] is not None
                    and item["beta_proteins"] is not None
                ):
                    for alpha_protein in item["alpha_proteins"]:
                        for beta_protein in item["beta_proteins"]:
                            alpha_string_id: str | None = (
                                proteins_to_string_ids[alpha_protein]
                                if alpha_protein in proteins_to_string_ids
                                else None
                            )
                            beta_string_id: str | None = (
                                proteins_to_string_ids[beta_protein]
                                if beta_protein in proteins_to_string_ids
                                else None
                            )
                            if (
                                alpha_string_id is not None
                                and beta_string_id is not None
                            ):
                                key = "_".join([alpha_string_id, beta_string_id])
                                if key in network:
                                    string_items.append(network[key])
                                    if network[key]["score"] is not None:
                                        string_scores.append(network[key]["score"])  # pyright: ignore[reportArgumentType] # ty: ignore[invalid-argument-type]
                if item["additional_information"] is None:
                    item.additional_information = dict()
                item["additional_information"][
                    "pyXLMS_annotated_STRING_interactions"
                ] = string_items
                item["additional_information"]["pyXLMS_annotated_STRING_score"] = (
                    max(string_scores) if len(string_scores) > 0 else float("nan")
                )
            return data
        else:
            raise TypeError(
                f"Can't annotate STRING scores for data type {type(data[0])}. Valid data types are:\n"
                "CrosslinkSpectrumMatch, Crosslink, and ParserResult."
            )
        return data
    new_csms = (
        assert_csms(
            annotate_string_scores(
                data["crosslink-spectrum-matches"], organism=organism, verbose=verbose
            )
        )
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        assert_xls(
            annotate_string_scores(
                data["crosslinks"],
                organism=organism,
                verbose=verbose,
            )
        )
        if data["crosslinks"] is not None
        else None
    )
    return create_parser_result(
        search_engine=data["search_engine"],
        csms=new_csms,
        crosslinks=new_xls,
    )
