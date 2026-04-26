#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import time
import requests

from typing import List
from typing import Dict
from typing import Any

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

# delete - this gives network len 28
EXAMPLES = ["CDC42", "CDK1", "KIF23", "PLK1", "RAC2", "RACGAP1", "RHOA", "RHOB"]


def __float_or_none(value: Any) -> float | None:
    try:
        return float(value)
    except Exception as _e:
        pass
    return None


def get_string_ids(proteins: List[str], organism: str | int) -> Dict[str, str | None]:
    if isinstance(organism, str):
        if organism not in ORGANISMS:
            raise KeyError()
        organism = ORGANISMS[organism]
    params = {
        "identifiers": "\r".join(proteins),
        "species": organism,
        "echo_query": 1,
        "caller_identity": CALLER_IDENTITY,
    }
    request_url = f"{STRING_STABLE_URL}/json/get_string_ids"
    response = requests.post(request_url, data=params)
    # wait one second after request to delay subsequent requests - be polite
    time.sleep(1)
    if not response.ok:
        raise RuntimeError()
    response_json = response.json()
    response_proteins: Dict[str, str] = dict()
    for item in response_json:
        if "queryItem" in item and "stringId" in item:
            response_proteins[str(item["queryItem"]).strip()] = str(
                item["stringId"]
            ).strip()
    output_proteins: Dict[str, str | None] = dict()
    for protein in proteins:
        if protein in response_proteins:
            output_proteins[protein] = response_proteins[protein]
        else:
            output_proteins[protein] = None
    return output_proteins


def get_string_network(string_ids: List[str], organism: str | int):
    if isinstance(organism, str):
        if organism not in ORGANISMS:
            raise KeyError()
        organism = ORGANISMS[organism]
    params = {
        "identifiers": "\r".join(string_ids),
        "species": organism,
        "required_score": 0,
        "add_nodes": 0,
        "show_query_node_labels": 0,
        "caller_identity": CALLER_IDENTITY,
    }
    request_url = f"{STRING_STABLE_URL}/json/network"
    response = requests.post(request_url, data=params)
    # wait one second after request to delay subsequent requests - be polite
    time.sleep(1)
    if not response.ok:
        raise RuntimeError()
    response_json = response.json()
    network: List[Dict[str, str | float | None]] = list()
    for item in response_json:
        parsed_item: Dict[str, str | float | None] = dict()
        parsed_item["A"] = str(item["stringId_A"]).strip()
        parsed_item["B"] = str(item["stringId_B"]).strip()
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
        network.append(parsed_item)
    return network
