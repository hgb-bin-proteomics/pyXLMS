#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import os
import json
from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._util import check_input_multi

from typing import Optional
from typing import BinaryIO
from typing import List


def to_json(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult,
    output_file: Optional[str | BinaryIO] = None,
    ensure_ascii: bool = False,
    indent: int = 4,
) -> str:
    _ok = check_input_multi(data, "data", [list, ParserResult])
    _ok = check_input(ensure_ascii, "ensure_ascii", bool)
    if output_file is not None:
        if isinstance(output_file, str):
            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=ensure_ascii, indent=indent)
        else:
            json.dump(data, output_file, ensure_ascii=ensure_ascii, indent=indent)  # ty: ignore[invalid-argument-type]
    return json.dumps(data, ensure_ascii=ensure_ascii, indent=indent)


def from_json(json_input: str | BinaryIO) -> None:
    if isinstance(json_input, str):
        if os.path.isfile(json_input):
            return
        else:
            return
    return
