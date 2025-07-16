#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from .parser import read
from .transform.summary import summary as transform_summary
from .transform.aggregate import unique as transform_unique
from .transform.validate import validate as transform_validate
from .transform.targets_only import targets_only as transform_targets_only

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
    unique: Optional[bool | Dict[str, Any]] = True,
    validate: Optional[bool | Dict[str, Any]] = True,
    targets_only: Optional[bool] = True,
    **kwargs,
) -> Dict[str, Any]:
    r"""Runs a standard down-stream analysis pipeline for crosslinks and crosslink-spectrum-matches.

    Runs a standard down-stream analysis pipeline for crosslinks and crosslink-spectrum-matches. The pipeline first reads
    a result file and subsequently optionally filters the the read data for unique crosslinks and crosslink-spectrum-matches,
    optionally the data is validated by false discovery rate estimation and - also optionally - only target-target matches
    are returned. Internally the pipeline calls ``parser.read()``, ``transform.unique()``, ``transform.validate()``, and
    ``transform.targets_only()``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the result file(s) or a file-like object/stream.
    engine : "Custom", "MaxQuant", "MaxLynx", "MS Annika", "mzIdentML", "pLink", "Scout", "xiSearch/xiFDR", or "XlinkX"
        Crosslink search engine or format of the result file.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    unique : dict of str, any, or bool, or None, default = True
        If ``transform.unique()`` should be run in the pipeline. If None or False this step is omitted.
        If True this step is run with default parameters. If a dictionary is given it should contain parameters for
        running ``transform.unique()``. Omitting a parameter in the dictionary will fall back to its default value.
    validate : dict of str, any, or bool, or None, default = True
        If ``transform.validate()`` should be run in the pipeline. If None or False this step is omitted.
        If True this step is run with default parameters. If a dictionary is given it should contain parameters for
        running ``transform.validate()``. Omitting a parameter in the dictionary will fall back to its default value.
    targets_only : bool, or None, default = True
        If ``transform.targets_only()`` should be run in the pipeline. If None or False this step is omitted.
    **kwargs
        Any additional parameters will be passed to the specific result file parsers.

    Returns
    -------
    dict of str, any
        The transformed parser_result after all pipeline steps are completed.

    Raises
    ------
    TypeError
        If any of the parameters do not have the correct type.

    Notes
    -----
    Various helpful pipeline information is also printed to ``stdout``.

    Examples
    --------
    >>> from pyXLMS.pipelines import pipeline
    >>> pr = pipeline("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
    ...               engine="MS Annika",
    ...               crosslinker="DSS",
    ...               unique=True,
    ...               validate={"fdr": 0.05, "formula":"(TD-DD)/TT"},
    ...               targets_only=True)
    Reading MS Annika CSMs...: 100%|██████████████████████████████████████████████████| 826/826 [00:00<00:00, 10337.98it/s]
    ---- Summary statistics before pipeline ----
    Number of CSMs: 826.0
    Number of unique CSMs: 826.0
    Number of intra CSMs: 803.0
    Number of inter CSMs: 23.0
    Number of target-target CSMs: 786.0
    Number of target-decoy CSMs: 39.0
    Number of decoy-decoy CSMs: 1.0
    Minimum CSM score: 1.11
    Maximum CSM score: 452.99
    Iterating over scores for FDR calculation...:   0%|                                            | 0/826 [00:00<?, ?it/s]
    ---- Summary statistics after pipeline ----
    Number of CSMs: 786.0
    Number of unique CSMs: 786.0
    Number of intra CSMs: 774.0
    Number of inter CSMs: 12.0
    Number of target-target CSMs: 786.0
    Number of target-decoy CSMs: 0.0
    Number of decoy-decoy CSMs: 0.0
    Minimum CSM score: 1.28
    Maximum CSM score: 452.99
    ---- Performed pipeline steps ----
    :: parser.read() ::
    :: parser.read() :: params :: <params omitted>
    :: transform.unique() ::
    :: transform.unique() :: params :: by=peptide
    :: transform.unique() :: params :: score=higher_better
    :: transform.validate() ::
    :: transform.validate() :: params :: fdr=0.05
    :: transform.validate() :: params :: formula=(TD-DD)/TT
    :: transform.validate() :: params :: score=higher_better
    :: transform.validate() :: params :: separate_intra_inter=False
    :: transform.validate() :: params :: ignore_missing_labels=False
    :: transform.targets_only() ::
    :: transform.targets_only() :: params :: no params
    """
    # steps: reading
    pr = read(files, engine=engine, crosslinker=crosslinker, **kwargs)
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
    # steps: pipeline summary
    print("---- Performed pipeline steps ----")
    print(":: parser.read() ::")
    print(":: parser.read() :: params :: <params omitted>")
    if unique is not None:
        if isinstance(unique, dict) or (isinstance(unique, bool) and unique):
            print(":: transform.unique() ::")
            for k, v in unique_params.items():
                print(f":: transform.unique() :: params :: {k}={v}")
    if validate is not None:
        if isinstance(validate, dict) or (isinstance(validate, bool) and validate):
            print(":: transform.validate() ::")
            for k, v in validate_params.items():
                print(f":: transform.validate() :: params :: {k}={v}")
    if targets_only is not None and targets_only:
        print(":: transform.targets_only() ::")
        print(":: transform.targets_only() :: params :: no params")
    # steps: finalize
    if not isinstance(pr, dict):
        raise RuntimeError(
            "Something went wrong while running the pipeline.\n"
            f"Expected data type: dict. Got: {type(pr)}."
        )
    return pr
