#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from .parser import read
from .transform_summary import summary as transform_summary
from .transform_aggregate import unique as transform_unique
from .transform_validate import validate as transform_validate
from .transform_targets_only import targets_only as transform_targets_only

from typing import Optional
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
    unique: Optional[bool | Dict[str, Any]] = None,
    validate: Optional[bool | Dict[str, Any]] = None,
    targets_only: Optional[bool] = None,
    **kwargs,
) -> Dict[str, Any]:
    # steps: reading
    pr = read(files, engine, crosslinker, **kwargs)
    # steps: summary (before)
    print("---- Summary statistics before pipeline ----")
    _ = transform_summary(pr)
    # steps: unique
    unique_params = {"by": "peptide", "score": "higher_better"}
    if unique is not None:
        if isinstance(unique, dict):
            unique_params.update(unique)
            pr = transform_unique(
                pr,
                by=str(unique_params["by"]),  # pyright: ignore[reportArgumentType]
                score=str(unique_params["score"]),  # pyright: ignore[reportArgumentType]
            )
        elif isinstance(unique, bool):
            if unique:
                pr = transform_unique(
                    pr,
                    by=str(unique_params["by"]),  # pyright: ignore[reportArgumentType]
                    score=str(unique_params["score"]),  # pyright: ignore[reportArgumentType]
                )
        else:
            raise TypeError(
                "Parameter unique has to be a dictionary of parameters for transform.unique(), a boolean or None!"
            )
    # steps: validate
    validate_params = {
        "fdr": 0.01,
        "formula": "D/T",
        "score": "higher_better",
        "separate_intra_inter": False,
        "ignore_missing_labels": False,
    }
    if validate is not None:
        if isinstance(validate, dict):
            validate_params.update(validate)
            pr = transform_validate(
                pr,
                fdr=float(validate_params["fdr"]),
                formula=str(validate_params["formula"]),  # pyright: ignore[reportArgumentType]
                score=str(validate_params["score"]),  # pyright: ignore[reportArgumentType]
                separate_intra_inter=bool(validate_params["separate_intra_inter"]),
                ignore_missing_labels=bool(validate_params["ignore_missing_labels"]),
            )
        elif isinstance(validate, bool):
            if validate:
                pr = transform_validate(
                    pr,
                    fdr=float(validate_params["fdr"]),
                    formula=str(validate_params["formula"]),  # pyright: ignore[reportArgumentType]
                    score=str(validate_params["score"]),  # pyright: ignore[reportArgumentType]
                    separate_intra_inter=bool(validate_params["separate_intra_inter"]),
                    ignore_missing_labels=bool(
                        validate_params["ignore_missing_labels"]
                    ),
                )
        else:
            raise TypeError(
                "Parameter validate has to be a dictionary of parameters for transform.validate(), a boolean or None!"
            )
    # steps: targets only
    if targets_only is not None:
        if isinstance(targets_only, bool):
            if targets_only:
                pr = transform_targets_only(pr)
        else:
            raise TypeError("Parameter targets_only has to be a boolean or None!")
    # steps: summary (after)
    print("---- Summary statistics after pipeline ----")
    _ = transform_summary(pr)
    if not isinstance(pr, dict):
        raise RuntimeError(
            "Something went wrong while running the pipeline.\n"
            f"Expected data type: dict. Got: {type(pr)}."
        )
    return pr
