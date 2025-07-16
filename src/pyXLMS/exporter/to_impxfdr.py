#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from ..data import check_input
from ..data import create_crosslink_from_csm
from ..transform.util import get_available_keys
from ..transform.filter import filter_target_decoy
from .to_msannika import to_msannika

from typing import Optional
from typing import Dict
from typing import Any
from typing import List


def to_impxfdr(
    data: List[Dict[str, Any]],
    filename: Optional[str],
    targets_only: bool = True,
) -> pd.DataFrame:
    r"""Exports a list of crosslinks or crosslink-spectrum-matches to IMP-X-FDR format.

    Exports a list of crosslinks or crosslink-spectrum-matches to IMP-X-FDR format for benchmarking purposes.
    The tool IMP-X-FDR is available from
    `github.com/vbc-proteomics-org/imp-x-fdr <https://github.com/vbc-proteomics-org/imp-x-fdr>`_.
    We recommend using version 1.1.0 and selecting "MS Annika" as input file format for the here exported file.
    A slightly modified version is available from
    `github.com/hgb-bin-proteomics/MSAnnika_NC_Results <https://github.com/hgb-bin-proteomics/MSAnnika_NC_Results/blob/master/Peplib_Beveridge/MS_Annika/Tools/IMP-X-FDR.v1.1.0.zip>`_.
    This version contains a few bug fixes and was used for the MS Annika 2.0 and MS Annika 3.0 publications.
    Requires that ``alpha_proteins``, ``beta_proteins``, ``alpha_proteins_crosslink_positions`` and ``beta_proteins_crosslink_positions`` fields
    are set for crosslinks and crosslink-spectrum-matches.

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslinks or crosslink-spectrum-matches.
    filename : str, or None, default = None
        If not None, the exported data will be written to a file with the specified filename.
        The filename should end in ".xlsx" as the file is exported to Microsoft Excel file format.
    targets_only : bool, default = True
        Whether or not only target crosslinks or crosslink-spectrum-matches should be exported. For
        benchmarking purposes this is usually the case. If the crosslinks or crosslink-spectrum-matches
        do not contain target-decoy labels this should be set to False.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing crosslinks or crosslink-spectrum-matches in IMP-X-FDR format.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If data contains elements of mixed data type.
    ValueError
        If the provided data contains no elements or if none of the data has target-decoy labels
        and parameter 'targets_only' is set to True.
    RuntimeError
        If not all of the required information is present in the input data.

    Examples
    --------
    >>> from pyXLMS.exporter import to_impxfdr
    >>> from pyXLMS.parser import read
    >>> pr = read("data/xi/1perc_xl_boost_Links_xiFDR2.2.1.csv", engine="xiSearch/xiFDR", crosslinker="DSS")
    >>> crosslinks = pr["crosslinks"]
    >>> to_impxfdr(crosslinks, filename="crosslinks.xlsx")
        Crosslink Type             Sequence A  Position A Accession A In protein A  ... Position B  Accession B In protein B Best CSM Score  Decoy
    0            Intra          VVDELV[K]VMGR           7        Cas9          753  ...          7         Cas9          753         40.679  False
    1            Intra  MLASAGELQ[K]GNELALPSK          10        Cas9          753  ...          7         Cas9         1226         40.231  False
    2            Intra        MDGTEELLV[K]LNR          10        Cas9          396  ...         10         Cas9          396         39.582  False
    3            Intra         MTNFD[K]NLPNEK           6        Cas9          965  ...          2         Cas9          504         35.880  False
    4            Intra             DFQFY[K]VR           6        Cas9          978  ...          4         Cas9         1028         35.281  False
    ..             ...                    ...         ...         ...          ...  ...        ...          ...          ...            ...    ...
    220          Intra        LP[K]YSLFELENGR           3        Cas9          866  ...          3         Cas9         1204          9.877  False
    221          Intra               D[K]QSGK           2        Cas9          677  ...          2         Cas9          677          9.702  False
    222          Intra               AGFI[K]R           5        Cas9          922  ...         11         Cas9          881          9.666  False
    223          Intra                E[K]IEK           2        Cas9          443  ...          1         Cas9          562          9.656  False
    224          Intra                LS[K]SR           3        Cas9          222  ...          3         Cas9          222          9.619  False
    [225 rows x 11 columns]

    >>> from pyXLMS.exporter import to_impxfdr
    >>> from pyXLMS.parser import read
    >>> pr = read("data/xi/1perc_xl_boost_CSM_xiFDR2.2.1.csv", engine="xiSearch/xiFDR", crosslinker="DSS")
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> to_impxfdr(csms, filename="csms.xlsx")
        Crosslink Type          Sequence A  Position A Accession A In protein A  ... Position B  Accession B In protein B Best CSM Score  Decoy
    0            Intra  [K]IECFDSVEISGVEDR           1        Cas9          575  ...          1         Cas9          575         27.268  False
    1            Intra       LVDSTD[K]ADLR           7        Cas9          152  ...         11         Cas9          881         26.437  False
    2            Intra     GGLSELD[K]AGFIK           8        Cas9          917  ...          8         Cas9          917         26.134  False
    3            Intra       LVDSTD[K]ADLR           7        Cas9          152  ...          7         Cas9          152         25.804  False
    4            Intra       VVDELV[K]VMGR           7        Cas9          753  ...          7         Cas9          753         24.861  False
    ..             ...                 ...         ...         ...          ...  ...        ...          ...          ...            ...    ...
    406          Intra          [K]GILQTVK           1        Cas9          739  ...          3         Cas9          222          6.977  False
    407          Intra          QQLPE[K]YK           6        Cas9          350  ...          6         Cas9          350          6.919  False
    408          Intra           ESILP[K]R           6        Cas9         1117  ...          7         Cas9         1035          6.853  False
    409          Intra             LS[K]SR           3        Cas9          222  ...          2         Cas9          884          6.809  False
    410          Intra     QIT[K]HVAQILDSR           4        Cas9          933  ...          6         Cas9          350          6.808  False
    [411 rows x 11 columns]
    """
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    _ok = check_input(targets_only, "targets_only", bool)
    if targets_only:
        data = filter_target_decoy(data)["Target-Target"]
    if len(data) == 0:
        if targets_only:
            raise ValueError(
                "Provided data does not contain any crosslinks or crosslink-spectrum-matches after filtering for targets only!"
            )
        else:
            raise ValueError(
                "Provided data does not contain any crosslinks or crosslink-spectrum-matches!"
            )
    if "data_type" not in data[0] or data[0]["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslink or crosslink-spectrum-match!"
        )
    available_keys = get_available_keys(data)
    if (
        not available_keys["alpha_proteins"]
        or not available_keys["alpha_proteins_crosslink_positions"]
        or not available_keys["beta_proteins"]
        or not available_keys["beta_proteins_crosslink_positions"]
    ):
        raise RuntimeError(
            "Can't export to IMP-X-FDR because not all necessary information is available!"
        )
    if data[0]["data_type"] == "crosslink":
        return to_msannika(data, filename, format="xlsx")
    return to_msannika(
        [create_crosslink_from_csm(csm) for csm in data], filename, format="xlsx"
    )
