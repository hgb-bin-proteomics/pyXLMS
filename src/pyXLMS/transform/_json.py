#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import os
import json

from pydantic import ValidationError

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._util import check_input_multi
from ._util import assert_csms_or_xls

from typing import Optional
from typing import BinaryIO
from typing import List
from typing import Dict
from typing import Any


def to_json(
    data: List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult,
    output_file: Optional[str | BinaryIO] = None,
    ensure_ascii: bool = False,
    indent: int = 4,
) -> str:
    r"""
    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import to_json, from_json
    >>> pr = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> pr.display()
    Data Type:                            parser_result
    Completeness:                         full
    Identifying Search Engine:            MS Annika
    Number of Crosslink-Spectrum-Matches: 826
    Number of Crosslinks:                 300
    >>> json_data_pr = to_json(pr)
    >>> json_data_csms = to_json(pr.csms())
    >>> json_data_xls = to_json(pr.xls(), output_file="xls.json")
    >>> pr = from_json(json_data_pr)
    >>> type(pr)
    <class 'pyXLMS.data._parser_result.ParserResult'>
    >>> csms = from_json(json_data_csms)
    >>> len(csms)
    826
    >>> type(csms[0])
    <class 'pyXLMS.data._csm.CrosslinkSpectrumMatch'>
    >>> xls = from_json("xls.json")
    >>> len(xls)
    300
    >>> type(xls[0])
    <class 'pyXLMS.data._crosslink.Crosslink'>
    """
    _ok = check_input_multi(data, "data", [list, ParserResult])
    _ok = check_input(ensure_ascii, "ensure_ascii", bool)
    json_data: List[Dict[str, Any]] | Dict[str, Any] | None = list()
    if isinstance(data, list):
        csms_or_xls = assert_csms_or_xls(data)
        for item in csms_or_xls:
            json_data.append(item.model_dump(mode="python"))
    else:
        json_data = data.model_dump(mode="python")
    if output_file is not None:
        if isinstance(output_file, str):
            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(json_data, f, ensure_ascii=ensure_ascii, indent=indent)
        else:
            json.dump(json_data, output_file, ensure_ascii=ensure_ascii, indent=indent)  # ty: ignore[invalid-argument-type]
    return json.dumps(json_data, ensure_ascii=ensure_ascii, indent=indent)


def from_json(
    json_input: str | BinaryIO,
) -> List[CrosslinkSpectrumMatch] | List[Crosslink] | ParserResult:
    r"""
    Examples
    --------
    >>> from pyXLMS.parser import read
    >>> from pyXLMS.transform import to_json, from_json
    >>> pr = read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> pr.display()
    Data Type:                            parser_result
    Completeness:                         full
    Identifying Search Engine:            MS Annika
    Number of Crosslink-Spectrum-Matches: 826
    Number of Crosslinks:                 300
    >>> json_data_pr = to_json(pr)
    >>> json_data_csms = to_json(pr.csms())
    >>> json_data_xls = to_json(pr.xls(), output_file="xls.json")
    >>> pr = from_json(json_data_pr)
    >>> type(pr)
    <class 'pyXLMS.data._parser_result.ParserResult'>
    >>> csms = from_json(json_data_csms)
    >>> len(csms)
    826
    >>> type(csms[0])
    <class 'pyXLMS.data._csm.CrosslinkSpectrumMatch'>
    >>> xls = from_json("xls.json")
    >>> len(xls)
    300
    >>> type(xls[0])
    <class 'pyXLMS.data._crosslink.Crosslink'>
    """
    json_data: List[Dict[str, Any]] | Dict[str, Any] | None = None
    if isinstance(json_input, str):
        if os.path.isfile(json_input):
            with open(json_input, "r", encoding="utf-8") as f:
                json_data = json.load(f)
        else:
            json_data = json.loads(json_input)
    else:
        json_input.seek(0)
        json_data = json.load(json_input)
    if json_data is None:
        raise ValueError()
    if isinstance(json_data, list):
        maybe_csms_or_xls = list()
        for item in json_data:
            maybe_csm_or_xl: CrosslinkSpectrumMatch | Crosslink | None = None
            try:
                maybe_csm_or_xl = CrosslinkSpectrumMatch.model_validate(
                    item, strict=False
                )
            except ValidationError as _e:
                maybe_csm_or_xl = Crosslink.model_validate(item, strict=False)
            if maybe_csm_or_xl is None:
                raise ValueError()
            maybe_csms_or_xls.append(maybe_csm_or_xl)
        return assert_csms_or_xls(maybe_csms_or_xls)
    return ParserResult.model_validate(json_data, strict=False)
