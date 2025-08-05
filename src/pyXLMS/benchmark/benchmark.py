#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from tqdm import tqdm

from . import beveridge2020

from typing import Optional
from typing import Dict
from typing import Any
from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def info(
    library: Literal["beveridge2020", "matzinger2022"], dataset: Optional[int] = None
) -> Dict[str, Any]:
    # input checks
    # if dataset none, give overview of datasets and trigger manual selection
    info_dict = dict()
    if library == "beveridge2020":
        info_dict["pride_url"] = beveridge2020.PRIDE_URL
        info_dict["proteomexchange_url"] = beveridge2020.PROTEOMEXCHANGE_URL
        info_dict["dataset_description"] = beveridge2020.DATASETS[dataset][
            "description"
        ]
    # print info
    # return info
    return info_dict


def get_labels(
    data: List[Dict[str, Any]],
    library: Literal["beveridge2020", "matzinger2022"],
    dataset: Optional[int] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    # input checks
    tp = list()
    fp = list()
    # get correct library via dataset
    # can also just call info() ideally and retrieve library that way
    for item in tqdm(
        data, total=len(data), desc=f"Labelling {data[0]['data_type']}s..."
    ):
        # check if crosslink / CSM is true positive or false positive
        # assign to list
        pass
    return {"True-Positive": tp, "False-Positive": fp}


def get_metrics(labelled_data: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Any]:
    metrics = dict()
    # e.g.
    #   number of tp
    #   number of fp
    #   experimentally validated FDR
    #   if scores -> calculate number of crosslinks/CSMs at experimentally validated FDR of 5% and 1%
    return metrics


def benchmark(
    data: List[Dict[str, Any]],
    library: Literal["beveridge2020", "matzinger2022"],
    dataset: Optional[int] = None,
) -> Dict[str, Any]:
    # input checks
    labels = get_labels(data, library, dataset)
    metrics = get_metrics(labels)
    # print metrics
    return {"Labelled-Data": labels, "Metrics": metrics}
