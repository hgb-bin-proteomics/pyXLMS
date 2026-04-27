#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import time
import requests
import warnings
from tqdm import tqdm

from ..data import check_input
from ..data import check_input_multi
from ..data import create_parser_result
from .filter import filter_crosslink_type
from .filter import filter_protein_distribution
from .util import assert_data_type_same

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
ORGANISMS = {
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
SCORES = {
    "low confidence": 0.15,
    "medium confidence": 0.4,
    "high confidence": 0.7,
    "highest confidence": 0.9,
}


def __float_or_none(value: Any) -> float | None:
    try:
        return float(value)
    except Exception as _e:
        pass
    return None


def get_string_ids(
    proteins: List[str], organism: str | int, verbose: Literal[0, 1, 2] = 1
) -> Dict[str, str | None]:
    _ok = check_input(proteins, "proteins", list, str)
    _ok = check_input_multi(organism, "organism", [str, int])
    if isinstance(organism, str):
        if organism not in ORGANISMS:
            raise KeyError(
                f"Could not resolve organism {organism}, please specify taxon identifier manually!"
            )
        organism = ORGANISMS[organism]
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
    _ok = check_input(string_ids, "string_ids", list, str)
    _ok = check_input_multi(organism, "organism", [str, int])
    if isinstance(organism, str):
        if organism not in ORGANISMS:
            raise KeyError(
                f"Could not resolve organism {organism}, please specify taxon identifier manually!"
            )
        organism = ORGANISMS[organism]
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
    data: List[Dict[str, Any]] | Dict[str, Any],
    organism: str | int,
    verbose: Literal[0, 1, 2] = 1,
) -> List[Dict[str, Any]] | Dict[str, Any]:
    r"""Annotates STRING interactions and scores for inter-links.

    Annotates STRING interactions and STRING scores for inter-links based on their associated proteins.
    Takes a list of crosslink-spectrum-matches or crosslinks, or a parser_result as input.

    Parameters
    ----------
    data : list of dict of str, any, or dict of str, any
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
    list of dict of str, any, or dict of str, any
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
    KeyError
        If the organism could not be resolved to a taxon identifier.

    Examples
    --------
    >>> from pyXLMS.transform import annotate_string_scores
    """
    _ok = check_input_multi(organism, "organism", [str, int])
    if isinstance(organism, str):
        if organism not in ORGANISMS:
            raise KeyError(
                f"Could not resolve organism {organism}, please specify taxon identifier manually!"
            )
        organism = ORGANISMS[organism]
    _ok = check_input(verbose, "verbose", int)
    if verbose not in [0, 1, 2]:
        raise TypeError("Verbose level has to be one of 0, 1, or 2!")
    if isinstance(data, list):
        _ok = check_input(data, "data", list, dict)
        if len(data) == 0:
            return data
        if "data_type" not in data[0]:
            raise TypeError(
                "Can't annotate STRING scores for input data. Input data has to be a list of crosslink-spectrum-matches or crosslinks "
                "or a 'parser_result'!"
            )
        _ok = assert_data_type_same(data)
        # annotate STRING scores
        if (
            data[0]["data_type"] == "crosslink"
            or data[0]["data_type"] == "crosslink-spectrum-match"
        ):
            inter = filter_crosslink_type(data)["Inter"]
            proteins = list(filter_protein_distribution(inter).keys())
            proteins_to_string_ids = get_string_ids(proteins, organism, verbose)
            string_ids: List[str] = list()
            for k, v in proteins_to_string_ids.items():
                if v is not None:
                    string_ids.append(v)
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
                item["additional_information"][
                    "pyXLMS_annotated_STRING_interactions"
                ] = string_items
                item["additional_information"]["pyXLMS_annotated_STRING_score"] = (
                    max(string_scores) if len(string_scores) > 0 else float("nan")
                )
            return data
        else:
            raise TypeError(
                f"Can't reannotate decoy labels for data type {data[0]['data_type']}. Valid data types are:\n"
                "'crosslink-spectrum-match', 'crosslink', and 'parser_result'."
            )
        return data
    _ok = check_input(data, "data", dict)
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError(
            "Can't annotate STRING scores for dict. Dict has to be a valid 'parser_result'!"
        )
    new_csms = (
        annotate_string_scores(
            data["crosslink-spectrum-matches"], organism=organism, verbose=verbose
        )
        if data["crosslink-spectrum-matches"] is not None
        else None
    )
    new_xls = (
        annotate_string_scores(
            data["crosslinks"],
            organism=organism,
            verbose=verbose,
        )
        if data["crosslinks"] is not None
        else None
    )
    if new_csms is not None:
        if not isinstance(new_csms, list):
            raise RuntimeError(
                "Something went wrong while annotating STRING scores.\n"
                f"Expected data type: list. Got: {type(new_csms)}."
            )
    if new_xls is not None:
        if not isinstance(new_xls, list):
            raise RuntimeError(
                "Something went wrong while annotating STRING scores.\n"
                f"Expected data type: list. Got: {type(new_xls)}."
            )
    return create_parser_result(
        search_engine=data["search_engine"], csms=new_csms, crosslinks=new_xls
    )
