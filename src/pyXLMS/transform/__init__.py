#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

r"""
Data summarization, transformation, and quality control functions.

Examples
--------
>>> from pyXLMS.parser import read
>>> from pyXLMS.transform import summary
>>> pr = read(
...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx",
...     engine="MS Annika",
...     crosslinker="DSS",
... )
>>> stats = summary(pr)
Number of crosslinks: 300.0
Number of unique crosslinks by peptide: 300.0
Number of unique crosslinks by protein: 298.0
Number of intra crosslinks: 279.0
Number of inter crosslinks: 21.0
Number of target-target crosslinks: 265.0
Number of target-decoy crosslinks: 0.0
Number of decoy-decoy crosslinks: 35.0
Minimum crosslink score: 1.11
Maximum crosslink score: 452.99
"""

__all__ = [
    "modifications_to_str",
    "assert_csms",
    "assert_xls",
    "assert_csms_or_xls",
    "assert_data_type_same",
    "get_available_keys",
    "check_available_keys",
    "filter_target_decoy",
    "filter_proteins",
    "filter_protein_distribution",
    "filter_crosslink_type",
    "filter_peptide_pair_distribution",
    "summary",
    "unique",
    "aggregate",
    "validate",
    "to_proforma",
    "to_dataframe",
    "from_dataframe",
    "targets_only",
    "fasta_title_to_accession",
    "reannotate_positions",
    "intersection",
    "annotate_fdr",
    "reannotate_decoy_labels",
    "filter_residue_pair_distribution",
    "get_string_ids",
    "get_string_network",
    "annotate_string_scores",
    "display",
    "to_json",
    "from_json",
]

from ._util import modifications_to_str
from ._util import assert_csms
from ._util import assert_xls
from ._util import assert_csms_or_xls
from ._util import assert_data_type_same
from ._util import get_available_keys
from ._util import check_available_keys
from ._filter import filter_target_decoy
from ._filter import filter_proteins
from ._filter import filter_protein_distribution
from ._filter import filter_crosslink_type
from ._filter import filter_peptide_pair_distribution
from ._summary import summary
from ._aggregate import unique
from ._aggregate import aggregate
from ._validate import validate
from ._to_proforma import to_proforma
from ._to_dataframe import to_dataframe
from ._from_dataframe import from_dataframe
from ._targets_only import targets_only
from ._reannotate_positions import fasta_title_to_accession
from ._intersection import intersection
from ._reannotate_positions import reannotate_positions
from ._annotate_fdr import annotate_fdr
from ._reannotate_decoy_labels import reannotate_decoy_labels
from ._filter import filter_residue_pair_distribution
from ._annotate_string_scores import get_string_ids
from ._annotate_string_scores import get_string_network
from ._annotate_string_scores import annotate_string_scores
from ._util import display
from ._json import to_json
from ._json import from_json
