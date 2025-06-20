#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from .parser import read
from .transform_aggregate import unique as transform_unique

from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def pipeline(
    files: str | List[str] | BinaryIO,
    engine: Literal[
        "Custom",
        "MaxQuant",
        "MaxLynx",
        "MS Annika",
        "mzIdentML",
        "pLink",
        "Scout",
        "xiSearch/xiFDR",
        "XlinkX",
    ],
    crosslinker: str,
    unique: Optional[bool|Dict[str, Any]] = None,
    **kwargs,
) -> Dict[str, Any]:
    pr = read(files, engine, crosslinker, **kwargs)
    unique_params = {"by": "peptide", "score": "higher_better"}
    if unique is not None:
        if isinstance(unique, dict):
            unique_params.update(unique)
            pr = transform_unique(pr, by=unique_params["by"], score=unique_params["score"])
        elif isinstance(unique, bool):
            if unique:
                pr = transform_unique(pr, by=unique_params["by"], score=unique_params["score"])
        else:
            raise TypeError()
    return pr
