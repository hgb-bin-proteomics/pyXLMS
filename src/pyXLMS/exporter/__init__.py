#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

r"""
Export crosslink-spectrum-matches and crosslinks to different down-stream analysis
tools and formats.

Examples
--------
>>> from pyXLMS.exporter import to_xiview
>>> from pyXLMS.parser import read
>>> from pyXLMS.transform import targets_only
>>> from pyXLMS.transform import filter_proteins
>>> pr = read(
...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
...     engine="MS Annika",
...     crosslinker="DSS",
... )
>>> crosslinks = targets_only(pr)["crosslinks"]
>>> cas9 = filter_proteins(crosslinks, proteins=["Cas9"])["Both"]
>>> to_xiview(cas9, filename="crosslinks_xiVIEW.csv")
    AbsPos1 AbsPos2 Protein1 Protein2 Decoy1 Decoy2   Score
0       779     779     Cas9     Cas9  FALSE  FALSE  119.83
1       866     866     Cas9     Cas9  FALSE  FALSE  114.43
2       677     677     Cas9     Cas9  FALSE  FALSE  200.98
3       677      48     Cas9     Cas9  FALSE  FALSE   94.47
4        34      34     Cas9     Cas9  FALSE  FALSE  110.48
..      ...     ...      ...      ...    ...    ...     ...
248     396     396     Cas9     Cas9  FALSE  FALSE  305.63
249     688     952     Cas9     Cas9  FALSE  FALSE  110.46
250     793    1180     Cas9     Cas9  FALSE  FALSE  288.36
251     575     688     Cas9     Cas9  FALSE  FALSE  376.15
252    1180    1180     Cas9     Cas9  FALSE  FALSE  437.10
[253 rows x 7 columns]
"""

__all__ = [
    "to_xmas",
    "to_xlinkdb",
    "to_impxfdr",
    "to_msannika",
    "get_msannika_crosslink_sequence",
    "to_pyxlinkviewer",
    "to_xlmstools",
    "to_xinet",
    "to_xiview",
    "to_xifdr",
    "to_alphalink2",
    "to_proxl",
]

from ._to_xmas import to_xmas
from ._to_xlinkdb import to_xlinkdb
from ._to_impxfdr import to_impxfdr
from ._to_msannika import to_msannika
from ._to_msannika import get_msannika_crosslink_sequence
from ._to_pyxlinkviewer import to_pyxlinkviewer
from ._to_xlmstools import to_xlmstools
from ._to_xinet import to_xinet
from ._to_xiview import to_xiview
from ._to_xifdr import to_xifdr
from ._to_alphalink2 import to_alphalink2
from ._to_proxl import to_proxl
