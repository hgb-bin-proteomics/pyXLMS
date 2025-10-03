#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from tqdm import tqdm

from ..data import check_input
from ..data import create_crosslink
from ..data import create_csm
from ..parser.parser_xldbse_custom import __get_value
from ..parser.parser_xldbse_custom import pyxlms_modification_str_parser
from ..parser.util import format_sequence
from ..parser.util import get_bool_from_value
from ..parser.util import __serialize_pandas_series

from typing import Optional
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List
from typing import Callable


def from_dataframe(
    df: pd.DataFrame,
    column_mapping: Optional[Dict[str, str]] = None,
    parse_modifications: bool = True,
    modification_parser: Optional[Callable[[str], Dict[int, Tuple[str, float]]]] = None,
    decoy_prefix: str = "REV_",
) -> List[Dict[str, Any]]:
    r"""Read a pandas DataFrame in custom or pyXLMS format.

    Reads a pandas DataFrame in custom or pyXLMS format and returns a list of crosslink-spectrum-matches or crosslinks.

    The minimum required columns for a crosslink-spectrum-matches pandas DataFrame are:

    - "Alpha Peptide": The unmodified amino acid sequence of the first peptide.
    - "Alpha Peptide Crosslink Position": The position of the crosslinker in the sequence of the first peptide (1-based).
    - "Beta Peptide": The unmodified amino acid sequence of the second peptide.
    - "Beta Peptide Crosslink Position": The position of the crosslinker in the sequence of the second peptide (1-based).
    - "Spectrum File": Name of the spectrum file the crosslink-spectrum-match was identified in.
    - "Scan Nr": The corresponding scan number of the crosslink-spectrum-match.

    The minimum required columns for crosslink pandas DataFrame are:

    - "Alpha Peptide": The unmodified amino acid sequence of the first peptide.
    - "Alpha Peptide Crosslink Position": The position of the crosslinker in the sequence of the first peptide (1-based).
    - "Beta Peptide": The unmodified amino acid sequence of the second peptide.
    - "Beta Peptide Crosslink Position": The position of the crosslinker in the sequence of the second peptide (1-based).

    A full specification of columns that can be parsed can be found in the
    `docs <https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md>`_.

    Parameters
    ----------
    df : pandas.DataFrame
        The pandas.DataFrame containing crosslink-spectrum-matches or crosslinks.
    column_mapping : dict of str, str
        A dictionary that maps the ``df`` columns to the required pyXLMS column names.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modification_parser' parameter.
    modification_parser : callable, or None
        A function that parses modification strings and returns the pyXLMS specific modifications object.
        If None, the function ``pyxlms_modification_str_parser()`` is used. If no modification columns are
        given this parameter is ignored.
    decoy_prefix : str, default = "REV\_"
        The prefix that indicates that a protein is from the decoy database.

    Returns
    -------
    list of dict of str, any
        If a crosslink-spectrum-matches DataFrame was given, a list of crosslink-spectrum-matches is returned.
        If a crosslinks DataFrame was given, a list of crosslinks is returned.

    Raises
    ------
    TypeError
        If one of the values could not be parsed.

    Examples
    --------
    >>> from pyXLMS import parser, transform
    >>> csms_from_pyxlms = parser.read_custom("data/pyxlms/csm.txt")
    >>> csms_df = transform.to_dataframe(csms_from_pyxlms["crosslink-spectrum-matches"])
    >>> csms_from_pyxlms = transform.from_dataframe(csms_df)

    >>> from pyXLMS import parser, transform
    >>> crosslinks_from_pyxlms = parser.read_custom("data/pyxlms/xl.txt")
    >>> crosslinks_df = transform.to_dataframe(crosslinks_from_pyxlms["crosslinks"])
    >>> crosslinks_from_pyxlms = transform.from_dataframe(crosslinks_df)
    """
    ## check input
    _ok = check_input(df, "df", pd.DataFrame)
    _ok = (
        check_input(column_mapping, "column_mapping", dict, str)
        if column_mapping is not None
        else True
    )
    _ok = check_input(parse_modifications, "parse_modifications", bool)
    _ok = (
        check_input(modification_parser, "modification_parser", Callable)
        if modification_parser is not None
        else True
    )
    _ok = check_input(decoy_prefix, "decoy_prefix", str)
    ## helper functions

    def get_is_decoy_value(
        row: pd.Series, decoy_prefix: str, alpha: bool
    ) -> bool | None:
        if alpha:
            if __get_value(row, "Alpha Decoy") is not None:
                return get_bool_from_value(__get_value(row, "Alpha Decoy"))
            if __get_value(row, "Alpha Proteins") is not None:
                return decoy_prefix in str(__get_value(row, "Alpha Proteins"))
            return None
        if __get_value(row, "Beta Decoy") is not None:
            return get_bool_from_value(__get_value(row, "Beta Decoy"))
        if __get_value(row, "Beta Proteins") is not None:
            return decoy_prefix in str(__get_value(row, "Beta Proteins"))
        return None

    def get_int(value: Any) -> int | None:
        if value is None:
            return None
        try:
            return int(value)
        except Exception as _e:
            pass
        raise TypeError(f"Could not parse int from value {value}!")
        return None

    def get_float(value: Any) -> float | None:
        if value is None:
            return None
        try:
            return float(value)
        except Exception as _e:
            pass
        raise TypeError(f"Could not parse float from value {value}!")
        return None

    ## set default parser
    if modification_parser is None:
        modification_parser = pyxlms_modification_str_parser

    ## data structures
    crosslinks = list()
    csms = list()

    ## detect input file type
    data = df.copy(deep=True)
    if column_mapping is not None:
        data.rename(columns=column_mapping, inplace=True)
    col_names = data.columns.values.tolist()
    is_crosslink_dataframe = "Scan Nr" not in col_names
    ## process data
    if is_crosslink_dataframe:
        for i, row in tqdm(
            data.iterrows(),
            total=data.shape[0],
            desc="Reading crosslinks...",
        ):
            # create crosslink
            crosslink = create_crosslink(
                peptide_a=format_sequence(str(row["Alpha Peptide"])),
                xl_position_peptide_a=int(row["Alpha Peptide Crosslink Position"]),
                proteins_a=[
                    protein.strip()
                    if protein.strip()[: len(decoy_prefix)] != decoy_prefix
                    else protein.strip()[len(decoy_prefix) :]
                    for protein in str(__get_value(row, "Alpha Proteins")).split(";")
                ]
                if __get_value(row, "Alpha Proteins") is not None
                else None,
                xl_position_proteins_a=[
                    int(position)
                    for position in str(
                        __get_value(row, "Alpha Proteins Crosslink Positions")
                    ).split(";")
                ]
                if __get_value(row, "Alpha Proteins Crosslink Positions") is not None
                else None,
                decoy_a=get_is_decoy_value(row, decoy_prefix, True),
                peptide_b=format_sequence(str(row["Beta Peptide"])),
                xl_position_peptide_b=int(row["Beta Peptide Crosslink Position"]),
                proteins_b=[
                    protein.strip()
                    if protein.strip()[: len(decoy_prefix)] != decoy_prefix
                    else protein.strip()[len(decoy_prefix) :]
                    for protein in str(__get_value(row, "Beta Proteins")).split(";")
                ]
                if __get_value(row, "Beta Proteins") is not None
                else None,
                xl_position_proteins_b=[
                    int(position)
                    for position in str(
                        __get_value(row, "Beta Proteins Crosslink Positions")
                    ).split(";")
                ]
                if __get_value(row, "Beta Proteins Crosslink Positions") is not None
                else None,
                decoy_b=get_is_decoy_value(row, decoy_prefix, False),
                score=get_float(__get_value(row, "Crosslink Score")),
                additional_information={"source": __serialize_pandas_series(row)},
            )
            crosslinks.append(crosslink)
        return crosslinks
    for i, row in tqdm(
        data.iterrows(),
        total=data.shape[0],
        desc="Reading CSMs...",
    ):
        # create csm
        csm = create_csm(
            peptide_a=format_sequence(str(row["Alpha Peptide"])),
            modifications_a=modification_parser(
                str(__get_value(row, "Alpha Peptide Modifications"))
            )
            if parse_modifications
            and __get_value(row, "Alpha Peptide Modifications") is not None
            else None,
            xl_position_peptide_a=int(row["Alpha Peptide Crosslink Position"]),
            proteins_a=[
                protein.strip()
                if protein.strip()[: len(decoy_prefix)] != decoy_prefix
                else protein.strip()[len(decoy_prefix) :]
                for protein in str(__get_value(row, "Alpha Proteins")).split(";")
            ]
            if __get_value(row, "Alpha Proteins") is not None
            else None,
            xl_position_proteins_a=[
                int(position)
                for position in str(
                    __get_value(row, "Alpha Proteins Crosslink Positions")
                ).split(";")
            ]
            if __get_value(row, "Alpha Proteins Crosslink Positions") is not None
            else None,
            pep_position_proteins_a=[
                int(position)
                for position in str(
                    __get_value(row, "Alpha Proteins Peptide Positions")
                ).split(";")
            ]
            if __get_value(row, "Alpha Proteins Peptide Positions") is not None
            else None,
            score_a=get_float(__get_value(row, "Alpha Score")),
            decoy_a=get_is_decoy_value(row, decoy_prefix, True),
            peptide_b=format_sequence(str(row["Beta Peptide"])),
            modifications_b=modification_parser(
                str(__get_value(row, "Beta Peptide Modifications"))
            )
            if parse_modifications
            and __get_value(row, "Beta Peptide Modifications") is not None
            else None,
            xl_position_peptide_b=int(row["Beta Peptide Crosslink Position"]),
            proteins_b=[
                protein.strip()
                if protein.strip()[: len(decoy_prefix)] != decoy_prefix
                else protein.strip()[len(decoy_prefix) :]
                for protein in str(__get_value(row, "Beta Proteins")).split(";")
            ]
            if __get_value(row, "Beta Proteins") is not None
            else None,
            xl_position_proteins_b=[
                int(position)
                for position in str(
                    __get_value(row, "Beta Proteins Crosslink Positions")
                ).split(";")
            ]
            if __get_value(row, "Beta Proteins Crosslink Positions") is not None
            else None,
            pep_position_proteins_b=[
                int(position)
                for position in str(
                    __get_value(row, "Beta Proteins Peptide Positions")
                ).split(";")
            ]
            if __get_value(row, "Beta Proteins Peptide Positions") is not None
            else None,
            score_b=get_float(__get_value(row, "Beta Score")),
            decoy_b=get_is_decoy_value(row, decoy_prefix, False),
            score=get_float(__get_value(row, "CSM Score")),
            spectrum_file=str(row["Spectrum File"]).strip(),
            scan_nr=int(row["Scan Nr"]),
            charge=get_int(__get_value(row, "Precursor Charge")),
            rt=get_float(__get_value(row, "Retention Time")),
            im_cv=get_float(__get_value(row, "Ion Mobility")),
            additional_information={"source": __serialize_pandas_series(row)},
        )
        csms.append(csm)
    return csms
