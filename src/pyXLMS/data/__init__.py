#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

r"""
Core data structures and data type validation functions.

Examples
--------
>>> from pyXLMS.data import CrosslinkSpectrumMatch as CSM
>>> csm = CSM(
...     alpha_peptide="PEKP",
...     alpha_peptide_crosslink_position=3,
...     beta_peptide="TKIDE",
...     beta_peptide_crosslink_position=2,
...     spectrum_file="dsso.mzML",
...     scan_nr=1,
... )

>>> from pyXLMS.data import Crosslink
>>> xl = Crosslink(
...     alpha_peptide="PEKP",
...     alpha_peptide_crosslink_position=3,
...     beta_peptide="TKIDE",
...     beta_peptide_crosslink_position=2,
... )

>>> from pyXLMS.data import Crosslink
>>> from pyXLMS.data import ParserResult
>>> xl = Crosslink(
...     alpha_peptide="PEKP",
...     alpha_peptide_crosslink_position=3,
...     beta_peptide="TKIDE",
...     beta_peptide_crosslink_position=2,
... )
>>> pr = ParserResult(search_engine="My Search Engine", crosslinks=[xl])
"""

__all__ = [
    "check_input",
    "check_input_multi",
    "check_indexing",
    "Crosslink",
    "create_crosslink",
    "create_crosslink_min",
    "create_crosslink_from_csm",
    "CrosslinkSpectrumMatch",
    "create_csm",
    "create_csm_min",
    "ParserResult",
    "create_parser_result",
]

from ._util import check_input
from ._util import check_input_multi
from ._util import check_indexing
from ._crosslink import Crosslink
from ._crosslink import create_crosslink
from ._crosslink import create_crosslink_min
from ._csm import CrosslinkSpectrumMatch
from ._csm import create_csm
from ._csm import create_csm_min
from ._csm import create_crosslink_from_csm
from ._parser_result import ParserResult
from ._parser_result import create_parser_result
