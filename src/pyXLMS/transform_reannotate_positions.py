#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from Bio import SeqIO

from .data import check_input
from .data import create_parser_result
from .transform_util import assert_data_type_same

from typing import Optional
from typing import Dict
from typing import Tuple
from typing import List
from typing import Any


def reannotate_positions(
    data: List[Dict[str, Any]] | Dict[str, Any],
    fasta: str
) -> List[Dict[str, Any]] | Dict[str, Any]:
    r"""
    """
    if type(data) == list:
        _ok = assert_data_type_same(data)
        reannoted = list()
        if data[0]["data_type"] == "crosslink":
            # todo
            pass
        elif data[0]["data_type"] == "crosslink-spectrum-match":
            # todo
            pass
        else:
            raise TypeError()
        return reannoted
    _ok = check_input(data, "data", dict)
    if "data_type" not in data or data["data_type"] != "parser_result":
        raise TypeError()
    return create_parser_result(
        search_engine = data["search_engine"],
        csms = reannotate_positions(data["crosslink-spectrum-matches"], fasta) if data["crosslink-spectrum-matches"] is not None else None,
        crosslinks= reannotate_positions(data["crosslinks"], fasta) if data["crosslinks"] is not None else None
    )
