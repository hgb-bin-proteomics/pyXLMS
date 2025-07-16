#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from ..data import check_input
from ..transform.util import assert_data_type_same
from .util import __get_filename

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
    r"""Returns the crosslinked peptide sequence in MS Annika format.

    Returns the crosslinked peptide sequence in MS Annika format, which is the peptide amino
    acid sequence with the crosslinked residue in square brackets (see examples).

    Parameters
    ----------
    peptide : str
        The (unmodified) amino acid sequence of the peptide.
    crosslink_position : int
        Position of the crosslinker in the peptide sequence (1-based).

    Returns
    -------
    str
        The crosslinked peptide sequence in MS Annika format.

    Raises
    ------
    ValueError
        If the crosslink position is outside the peptide's length.

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
    r"""Helper function to get the [Alpha|Beta] T/D value.

    Parameters
    ----------
    value : bool, or None
        Decoy value of the crosslink-spectrum-match, should be either "alpha_decoy" or "beta_decoy" attribute.

    Returns
    -------
    str, or None
        If None was provided, None is returned. If a boolean is provided, returns "D" if True or "T" if False.

    Notes
    -----
    This function should not be called directly, it is called from ``__csms_to_msannika()``.
    """
    _ok = check_input(value, "value", bool) if value is not None else True
    if value is None:
        return None
    if value:
        return "D"
    return "T"


def __get_xl_isdecoy(
    alpha_decoy: Optional[bool], beta_decoy: Optional[bool]
) -> bool | None:
    r"""Helper function to get the Decoy value.

    Parameters
    ----------
    alpha_decoy : bool, or None
        Decoy value for the alpha peptide of the crosslink, should be "alpha_decoy" attribute.
    beta_decoy : bool, or None
        Decoy value for the beta peptide of the crosslink, should be "beta_decoy" attribute.

    Returns
    -------
    bool, or None
        If None was provided for any of the inputs, None is returned. Otherwise returns True if any of the inputs
        is True (= a decoy).

    Notes
    -----
    This function should not be called directly, it is called from ``__xls_to_msannika()``.
    """
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
    r"""Exports crosslink-spectrum-matches to MS Annika format.

    Parameters
    ----------
    csms : list of dict of str, any
        A list of crosslink-spectrum-matches.
    filename : str, or None
        If not None, the data will be written to a file with the specified filename.
    format : str, one of "csv", "tsv", or "xlsx"
        Format of the output file if filename is not None.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame in MS Annika format.

    Notes
    -----
    This function should not be called directly, it is called from ``to_msannika()``.
    """
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
            msannika_df.to_csv(__get_filename(filename, format), index=False)
        elif format == "tsv":
            msannika_df.to_csv(__get_filename(filename, format), sep="\t", index=False)
        else:
            msannika_df.to_excel(
                __get_filename(filename, format), engine="openpyxl", index=False
            )
    return msannika_df


def __xls_to_msannika(
    xls: List[Dict[str, Any]],
    filename: Optional[str],
    format: Literal["csv", "tsv", "xlsx"],
) -> pd.DataFrame:
    r"""Exports crosslinks to MS Annika format.

    Parameters
    ----------
    xls : list of dict of str, any
        A list of crosslinks.
    filename : str, or None
        If not None, the data will be written to a file with the specified filename.
    format : str, one of "csv", "tsv", or "xlsx"
        Format of the output file if filename is not None.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame in MS Annika format.

    Notes
    -----
    This function should not be called directly, it is called from ``to_msannika()``.
    """
    crosslink_type = list()
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
        crosslink_type.append("Intra" if xl["crosslink_type"] == "intra" else "Inter")
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
            "Crosslink Type": crosslink_type,
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
            msannika_df.to_csv(__get_filename(filename, format), index=False)
        elif format == "tsv":
            msannika_df.to_csv(__get_filename(filename, format), sep="\t", index=False)
        else:
            msannika_df.to_excel(
                __get_filename(filename, format), engine="openpyxl", index=False
            )
    return msannika_df


def to_msannika(
    data: List[Dict[str, Any]],
    filename: Optional[str] = None,
    format: Literal["csv", "tsv", "xlsx"] = "csv",
) -> pd.DataFrame:
    r"""Exports a list of crosslinks or crosslink-spectrum-matches to MS Annika format.

    Exports a list of crosslinks or crosslink-spectrum-matches to MS Annika format. This might be useful
    for tools that support MS Annika input but are not supported by pyXLMS (yet).

    Parameters
    ----------
    data : list of dict of str, any
        A list of crosslinks or crosslink-spectrum-matches.
    filename : str, or None, default = None
        If not None, the exported data will be written to a file with the specified filename.
    format : str, one of "csv", "tsv", or "xlsx", default = "csv"
        File format of the exported file if filename is not None.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing crosslinks or crosslink-spectrum-matches in MS Annika format.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If data contains elements of mixed data type.
    TypeError
        If parameter format is not one of 'csv', 'tsv' or 'xlsx'.
    ValueError
        If the provided data contains no elements.

    Warnings
    --------
    The MS Annika exporter will not check if all necessary information is available for the exported
    crosslinks or crosslink-spectrum-matches. If a value is not available it will be denoted as a missing
    value in the dataframe and exported file. Please make sure all necessary information is available
    before using the exported file with another tool! Please also note that modifications are not exported,
    for modification down-stream analysis please refer to ``transform.to_proforma()`` or
    ``transform.to_dataframe()``!

    Examples
    --------
    >>> from pyXLMS.exporter import to_msannika
    >>> from pyXLMS.data import create_crosslink_min
    >>> xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    >>> xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    >>> crosslinks = [xl1, xl2]
    >>> to_msannika(crosslinks)
      Crosslink Type  Sequence A  Position A Accession A In protein A  Sequence B  Position B Accession B In protein B Best CSM Score Decoy
    0          Inter  [K]PEPTIDE           1        None         None  P[K]EPTIDE           2        None         None           None  None
    1          Inter  PE[K]PTIDE           3        None         None  PEP[K]TIDE           4        None         None           None  None

    >>> from pyXLMS.exporter import to_msannika
    >>> from pyXLMS.data import create_crosslink_min
    >>> xl1 = create_crosslink_min("KPEPTIDE", 1, "PKEPTIDE", 2)
    >>> xl2 = create_crosslink_min("PEKPTIDE", 3, "PEPKTIDE", 4)
    >>> crosslinks = [xl1, xl2]
    >>> df = to_msannika(crosslinks, filename = "crosslinks.csv", format = "csv")

    >>> from pyXLMS.exporter import to_msannika
    >>> from pyXLMS.data import create_csm_min
    >>> csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    >>> csm2 = create_csm_min("PEKPTIDE", 3, "PEPKTIDE", 4, "RUN_1", 2)
    >>> csms = [csm1, csm2]
    >>> to_msannika(csms)
                Sequence Crosslink Type Sequence A  Crosslinker Position A  ... First Scan Charge RT [min] Compensation Voltage
    0  KPEPTIDE-PKEPTIDE          Inter   KPEPTIDE                       1  ...          1   None     None                 None
    1  PEKPTIDE-PEPKTIDE          Inter   PEKPTIDE                       3  ...          2   None     None                 None
    [2 rows x 20 columns]

    >>> from pyXLMS.exporter import to_msannika
    >>> from pyXLMS.data import create_csm_min
    >>> csm1 = create_csm_min("KPEPTIDE", 1, "PKEPTIDE", 2, "RUN_1", 1)
    >>> csm2 = create_csm_min("PEKPTIDE", 3, "PEPKTIDE", 4, "RUN_1", 2)
    >>> csms = [csm1, csm2]
    >>> df = to_msannika(csms, filename = "csms.csv", format = "csv")
    """
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
