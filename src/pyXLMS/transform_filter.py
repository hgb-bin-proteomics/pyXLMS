#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .data import check_input_multi

from typing import Optional
from typing import Dict
from typing import List
from typing import Set
from typing import Any


def filter_target_decoy(data: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    r"""
    """
    _ok = check_input(data, "data", list, dict)
    tt = list()
    td = list()
    dd = list()
    for item in data:
        if "data_type" not in item or item["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
            )
        if item["alpha_decoy"] is not None and item["beta_decoy"] is not None:
            if item["alpha_decoy"] and item["beta_decoy"]:
                dd.append(item)
            elif not item["alpha_decoy"] and not item["beta_decoy"]:
                tt.append(item)
            else:
                td.append(item)
    return {"Target-Target": tt, "Target-Decoy": td, "Decoy-Decoy": dd}


def filter_protein(data: List[Dict[str, Any]], proteins: Set[str] | List[str]) -> Dict[str, List[Any]]:
    r"""
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input_multi(proteins, "proteins", [set, list], str)
    proteins = set(proteins)
    intra = list()
    inter = list()
    for item in data:
        if "data_type" not in item or item["data_type"] not in [
            "crosslink",
            "crosslink-spectrum-match",
        ]:
            raise TypeError(
                "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
            )
        if item["alpha_proteins"] is not None and item["beta_proteins"] is not None:
            a = set(item["alpha_proteins"])
            b = set(item["beta_proteins"])
            if len(proteins.intersection(a)) > 0 and len(proteins.intersection(b)) > 0:
                intra.append(item)
            elif len(proteins.intersection(a)) == 0 and len(proteins.intersection(b)) == 0:
                continue
            else:
                inter.append(item)
    return {"Proteins": list(proteins), "Intra": intra, "Inter": inter}
