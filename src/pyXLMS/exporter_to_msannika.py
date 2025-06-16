#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from .data import check_input
from .transform_util import assert_data_type_same

from typing import Optional
from typing import Dict
from typing import Any
from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def get_msannika_crosslink_sequence(peptide: str, crosslink_position: int) -> str:
    r"""
    Examples
    --------
    >>> from pyXLMS.exporter import get_msannika_crosslink_sequence
    >>> get_msannika_crosslink_sequence("PEPKTIDE", 4)
    'PEP[K]TIDE'

    >>> from pyXLMS.exporter import get_msannika_crosslink_sequence
    >>> get_msannika_crosslink_sequence("KPEPTIDE", 1)
    '[K]PEPTIDE'

    >>> from pyXLMS.exporter import get_msannika_crosslink_sequence
    >>> get_msannika_crosslink_sequence("PEPTIDEK", 8)
    'PEPTIDE[K]'
    """
    if crosslink_position < 1 or crosslink_position > len(peptide):
        raise ValueError(
            f"Crosslink position outside of range! Must be in range [1, {len(peptide)}]."
        )
    return f"{peptide[: crosslink_position - 1]}[{peptide[crosslink_position - 1]}]{peptide[crosslink_position:]}"


def __get_csm_td(value: Optional[bool]) -> str | None:
    _ok = check_input(value, "value", bool) if value is not None else True
    if value is None:
        return None
    if value:
        return "D"
    return "T"


def __get_xl_isdecoy(
    alpha_decoy: Optional[bool], beta_decoy: Optional[bool]
) -> bool | None:
    _ok = (
        check_input(alpha_decoy, "alpha_decoy", bool)
        if alpha_decoy is not None
        else None
    )
    _ok = (
        check_input(beta_decoy, "beta_decoy", bool) if beta_decoy is not None else None
    )
    if alpha_decoy is None or beta_decoy is None:
        return None
    return alpha_decoy or beta_decoy


def __csms_to_msannika(
    csms: List[Dict[str, Any]],
    filename: Optional[str],
    format: Literal["csv", "tsv", "xlsx"],
) -> pd.DataFrame:
    sequence = list()
    crosslink_type = list()
    sequence_a = list()
    crosslinker_position_a = list()
    accession_a = list()
    a_in_protein = list()
    score_alpha = list()
    alpha_td = list()
    sequence_b = list()
    crosslinker_position_b = list()
    accession_b = list()
    b_in_protein = list()
    score_beta = list()
    beta_td = list()
    combined_score = list()
    spectrum_file = list()
    first_scan = list()
    charge = list()
    rt_min = list()
    compensation_voltage = list()
    for csm in csms:
        sequence.append(f"{csm['alpha_peptide']}-{csm['beta_peptide']}")
        crosslink_type.append("Intra" if csm["crosslink_type"] == "intra" else "Inter")
        sequence_a.append(csm["alpha_peptide"])
        crosslinker_position_a.append(csm["alpha_peptide_crosslink_position"])
        accession_a.append(
            ";".join(csm["alpha_proteins"])
            if csm["alpha_proteins"] is not None
            else None
        )
        a_in_protein.append(
            ";".join([str(pos - 1) for pos in csm["alpha_proteins_peptide_positions"]])
            if csm["alpha_proteins_peptide_positions"] is not None
            else None
        )
        score_alpha.append(csm["alpha_score"])
        alpha_td.append(__get_csm_td(csm["alpha_decoy"]))
        sequence_b.append(csm["beta_peptide"])
        crosslinker_position_b.append(csm["beta_peptide_crosslink_position"])
        accession_b.append(
            ";".join(csm["beta_proteins"]) if csm["beta_proteins"] is not None else None
        )
        b_in_protein.append(
            ";".join([str(pos - 1) for pos in csm["beta_proteins_peptide_positions"]])
            if csm["beta_proteins_peptide_positions"] is not None
            else None
        )
        score_beta.append(csm["beta_score"])
        beta_td.append(__get_csm_td(csm["beta_decoy"]))
        combined_score.append(csm["score"])
        spectrum_file.append(csm["spectrum_file"])
        first_scan.append(csm["scan_nr"])
        charge.append(csm["charge"])
        rt_min.append(
            csm["retention_time"] / 60.0 if csm["retention_time"] is not None else None
        )
        compensation_voltage.append(csm["ion_mobility"])
    msannika_df = pd.DataFrame(
        {
            "Sequence": sequence,
            "Crosslink Type": crosslink_type,
            "Sequence A": sequence_a,
            "Crosslinker Position A": crosslinker_position_a,
            "Accession A": accession_a,
            "A in protein": a_in_protein,
            "Score Alpha": score_alpha,
            "Alpha T/D": alpha_td,
            "Sequence B": sequence_b,
            "Crosslinker Position B": crosslinker_position_b,
            "Accession B": accession_b,
            "B in protein": b_in_protein,
            "Score Beta": score_beta,
            "Beta T/D": beta_td,
            "Combined Score": combined_score,
            "Spectrum File": spectrum_file,
            "First Scan": first_scan,
            "Charge": charge,
            "RT [min]": rt_min,
            "Compensation Voltage": compensation_voltage,
        }
    )
    if filename is not None:
        if format == "csv":
            msannika_df.to_csv(filename, index=False)
        elif format == "tsv":
            msannika_df.to_csv(filename, sep="\t", index=False)
        else:
            msannika_df.to_excel(filename, engine="openpyxl", index=False)
    return msannika_df


def __xls_to_msannika(
    xls: List[Dict[str, Any]],
    filename: Optional[str],
    format: Literal["csv", "tsv", "xlsx"],
) -> pd.DataFrame:
    sequence_a = list()
    position_a = list()
    accession_a = list()
    in_protein_a = list()
    sequence_b = list()
    position_b = list()
    accession_b = list()
    in_protein_b = list()
    best_csm_score = list()
    decoy = list()
    for xl in xls:
        sequence_a.append(
            get_msannika_crosslink_sequence(
                xl["alpha_peptide"], xl["alpha_peptide_crosslink_position"]
            )
        )
        position_a.append(xl["alpha_peptide_crosslink_position"])
        accession_a.append(
            ";".join(xl["alpha_proteins"]) if xl["alpha_proteins"] is not None else None
        )
        in_protein_a.append(
            ";".join([str(pos) for pos in xl["alpha_proteins_crosslink_positions"]])
            if xl["alpha_proteins_crosslink_positions"] is not None
            else None
        )
        sequence_b.append(
            get_msannika_crosslink_sequence(
                xl["beta_peptide"], xl["beta_peptide_crosslink_position"]
            )
        )
        position_b.append(xl["beta_peptide_crosslink_position"])
        accession_b.append(
            ";".join(xl["beta_proteins"]) if xl["beta_proteins"] is not None else None
        )
        in_protein_b.append(
            ";".join([str(pos) for pos in xl["beta_proteins_crosslink_positions"]])
            if xl["beta_proteins_crosslink_positions"] is not None
            else None
        )
        best_csm_score.append(xl["score"])
        decoy.append(__get_xl_isdecoy(xl["alpha_decoy"], xl["beta_decoy"]))
    msannika_df = pd.DataFrame(
        {
            "Sequence A": sequence_a,
            "Position A": position_a,
            "Accession A": accession_a,
            "In protein A": in_protein_a,
            "Sequence B": sequence_b,
            "Position B": position_b,
            "Accession B": accession_b,
            "In protein B": in_protein_b,
            "Best CSM Score": best_csm_score,
            "Decoy": decoy,
        }
    )
    if filename is not None:
        if format == "csv":
            msannika_df.to_csv(filename, index=False)
        elif format == "tsv":
            msannika_df.to_csv(filename, sep="\t", index=False)
        else:
            msannika_df.to_excel(filename, engine="openpyxl", index=False)
    return msannika_df


def to_msannika(
    data: List[Dict[str, Any]],
    filename: Optional[str] = None,
    format: Literal["csv", "tsv", "xlsx"] = "csv",
) -> pd.DataFrame:
    _ok = check_input(data, "data", list, dict)
    _ok = check_input(filename, "filename", str) if filename is not None else True
    _ok = check_input(format, "format", str)
    if format not in ["csv", "tsv", "xlsx"]:
        raise TypeError("Parameter 'format' has to be one of 'csv', 'tsv', or 'xlsx'!")
    if len(data) == 0:
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
    if not assert_data_type_same(data):
        raise TypeError("Not all elements in data have the same data type!")
    if data[0]["data_type"] == "crosslink":
        return __xls_to_msannika(data, filename, format)
    return __csms_to_msannika(data, filename, format)
