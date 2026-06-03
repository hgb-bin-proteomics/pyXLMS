#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

r"""
Plot different visualizations of crosslink-spectrum-match and crosslink data.

Examples
--------
>>> from pyXLMS import parser
>>> from pyXLMS import plotting
>>> pr = parser.read_msannika(
...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
... )
>>> csms = pr["crosslink-spectrum-matches"]
>>> fig, ax = plotting.plot_crosslink_type_distribution(csms)
"""

__all__ = [
    "venn",
    "plot_venn_diagram",
    "plot_score_distribution",
    "plot_target_decoy_distribution",
    "plot_protein_distribution",
    "plot_peptide_pair_distribution",
    "plot_crosslink_type_distribution",
    "plot_residue_pair_distribution",
    "plot_string_score_distribution",
]

from ._plot_venn_diagram import venn
from ._plot_venn_diagram import plot_venn_diagram
from ._plot_score_distribution import plot_score_distribution
from ._plot_target_decoy_distribution import plot_target_decoy_distribution
from ._plot_protein_distribution import plot_protein_distribution
from ._plot_peptide_pair_distribution import plot_peptide_pair_distribution
from ._plot_crosslink_type_distribution import plot_crosslink_type_distribution
from ._plot_residue_pair_distribution import plot_residue_pair_distribution
from ._plot_string_score_distribution import plot_string_score_distribution
