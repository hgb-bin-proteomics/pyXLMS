#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data._csm import CrosslinkSpectrumMatch
from ..data._crosslink import Crosslink
from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._util import check_input_multi

from typing import Optional
from typing import Any
from typing import Dict
from typing import Tuple
from typing import List

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def modifications_to_str(
    modifications: Optional[Dict[int, Tuple[str, float]]],
) -> str | None:
    r"""Returns the string representation of a modifications dictionary.

    Parameters
    ----------
    modifications : dict of [str, tuple], or None
        The modifications of a peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.
        ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications should be denoted with position ``len(peptide) + 1``.

    Returns
    -------
    str, or None
        The string representation of the modifications (or ``None`` if no modification was provided).

    Examples
    --------
    >>> from pyXLMS.transform import modifications_to_str
    >>> modifications_to_str(
    ...     {1: ("Oxidation", 15.994915), 5: ("Carbamidomethyl", 57.021464)}
    ... )
    '(1:[Oxidation|15.994915]);(5:[Carbamidomethyl|57.021464])'
    """
    ## check input
    _ok = (
        check_input(modifications, "modifcations", dict, tuple)
        if modifications is not None
        else True
    )

    modifications_str = ""
    if modifications is None:
        return None
    for modification_pos in sorted(modifications.keys()):
        modifications_str += f"({modification_pos}:[{modifications[modification_pos][0]}|{modifications[modification_pos][1]}]);"
    return modifications_str.rstrip(";")


def assert_csms(maybe_csms: Any) -> List[CrosslinkSpectrumMatch]:
    r"""Checks that the provided input is a list of type CrosslinkSpectrumMatch.

    Parameters
    ----------
    maybe_csms : any
        The input data to be checked.

    Returns
    -------
    list of CrosslinkSpectrumMatch
        Returns a list of type CrosslinkSpectrumMatch if the provided data was one.

    Raises
    ------
    TypeError
        If the provided data was not a list of CrosslinkSpectrumMatch.
    """
    csms: List[CrosslinkSpectrumMatch] = list()
    if isinstance(maybe_csms, list):
        for item in maybe_csms:
            if isinstance(item, CrosslinkSpectrumMatch):
                csms.append(item)
            else:
                raise TypeError(
                    "Provided input is not a valid list of type CrosslinkSpectrumMatch!"
                )
        return csms
    raise TypeError(
        "Provided input is not a valid list of type CrosslinkSpectrumMatch!"
    )
    return csms


def assert_xls(maybe_xls: Any) -> List[Crosslink]:
    r"""Checks that the provided input is a list of type Crosslink.

    Parameters
    ----------
    maybe_xls : any
        The input data to be checked.

    Returns
    -------
    list of Crosslink
        Returns a list of type Crosslink if the provided data was one.

    Raises
    ------
    TypeError
        If the provided data was not a list of Crosslink.
    """
    xls: List[Crosslink] = list()
    if isinstance(maybe_xls, list):
        for item in maybe_xls:
            if isinstance(item, Crosslink):
                xls.append(item)
            else:
                raise TypeError("Provided input is not a valid list of type Crosslink!")
        return xls
    raise TypeError("Provided input is not a valid list of type Crosslink!")
    return xls


def assert_csms_or_xls(
    maybe_csms_or_xls: Any,
) -> List[CrosslinkSpectrumMatch] | List[Crosslink]:
    r"""Checks that the provided input is a list of type CrosslinkSpectrumMatch or Crosslink.

    Parameters
    ----------
    maybe_csms_or_xls : any
        The input data to be checked.

    Returns
    -------
    list of CrosslinkSpectrumMatch, or list of Crosslink
        Returns a list of type CrosslinkSpectrumMatch, or a list of type Crosslink
        if the provided data was either.

    Raises
    ------
    TypeError
        If the provided data was neither a list of CrosslinkSpectrumMatch nor a list
        of Crosslink.
    """
    if isinstance(maybe_csms_or_xls, list):
        if len(maybe_csms_or_xls) == 0:
            return []
        if all(isinstance(item, CrosslinkSpectrumMatch) for item in maybe_csms_or_xls):
            return assert_csms(maybe_csms_or_xls)
        if all(isinstance(item, Crosslink) for item in maybe_csms_or_xls):
            return assert_xls(maybe_csms_or_xls)
    raise TypeError(
        "Provided input is not a valid list of type CrosslinkSpectrumMatch or Crosslink!"
    )
    return []


def assert_data_type_same(
    data_list: List[CrosslinkSpectrumMatch] | List[Crosslink] | List[ParserResult],
) -> bool:
    r"""Checks that all data is of the same data type.

    Verifies that all elements in the provided list are of the same data type.

    Parameters
    ----------
    data_list : list of dict of str, any
        A list of dictionaries with the ``data_type`` key.

    Returns
    -------
    bool
        If all elements are of the same data type.

    Raises
    ------
    TypeError
        If the item in the data list are not of type CrosslinkSpectrumMatch, Crosslink, or ParserResult.

    Examples
    --------
    >>> from pyXLMS.transform import assert_data_type_same
    >>> from pyXLMS import data
    >>> data_list = [
    ...     data.create_crosslink_min("PEPK", 4, "PKEP", 2),
    ...     data.create_crosslink_min("KPEP", 1, "PEKP", 3),
    ... ]
    >>> assert_data_type_same(data_list)
    True

    >>> from pyXLMS.transform import assert_data_type_same
    >>> from pyXLMS import data
    >>> data_list = [
    ...     data.create_crosslink_min("PEPK", 4, "PKEP", 2),
    ...     data.create_csm_min("KPEP", 1, "PEKP", 3, "RUN_1", 1),
    ... ]
    >>> assert_data_type_same(data_list)
    False
    """
    _ok = check_input(data_list, "data_list", list)
    if len(data_list) == 0:
        return True
    data_type = type(data_list[0])
    for item in data_list[1:]:
        if not isinstance(item, data_type):
            return False
    if isinstance(data_list[0], CrosslinkSpectrumMatch):
        return True
    if isinstance(data_list[0], Crosslink):
        return True
    if isinstance(data_list[0], ParserResult):
        return True
    raise TypeError(
        "Input list contains elements that are not of type CrosslinkSpectrumMatch, Crosslink, or ParserResult!"
    )
    return False


def get_available_keys(
    data_list: List[CrosslinkSpectrumMatch] | List[Crosslink],
    always_revalidate: bool = True,
) -> Dict[str, bool]:
    r"""Checks which data is available from a list of crosslinks or crosslink-spectrum-matches.

    Verifies which data fields have been set for all crosslinks or crosslink-spectrum-matches in the
    given list. Will return a dictionary structured the same as a crosslink or crosslink-spectrum-match,
    but instead of the data it will return either True or False, depending if the field was set or not.

    Parameters
    ----------
    data_list : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslinks or crosslink-spectrum-matches.
    always_revalidate : bool, default = True
        If ``True`` (default) the assigned ``completeness`` will be ignored and all data fields
        are re-checked. This is safer especially when data has been modified post reading.

    Returns
    -------
    dict of str, bool
        - If a list of crosslinks was provided, a dictionary with the following keys will be returned, where the value
          of each key denotes if the data field is available for all crosslinks in ``data_list``.
          Keys: ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_peptide_crosslink_position``,
          ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_decoy``, ``beta_peptide``, ``beta_peptide_crosslink_position``,
          ``beta_proteins``, ``beta_proteins_crosslink_positions``, ``beta_decoy``, ``crosslink_type``, ``score``, and ``additional_information``.
        - If a list of crosslink-spectrum-matches was provided, a dictionary with the following keys will be returned, where the value
          of each key denotes if the data field is available for all crosslink-spectrum-matches in ``data_list``.
          Keys: ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_modifications``,
          ``alpha_peptide_crosslink_position``, ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_proteins_peptide_positions``,
          ``alpha_score``, ``alpha_decoy``, ``beta_peptide``, ``beta_modifications``, ``beta_peptide_crosslink_position``, ``beta_proteins``,
          ``beta_proteins_crosslink_positions``, ``beta_proteins_peptide_positions``, ``beta_score``, ``beta_decoy``, ``crosslink_type``, ``score``,
          ``spectrum_file``, ``scan_nr``, ``retention_time``, ``ion_mobility``, and ``additional_information``.

    Raises
    ------
    TypeError
        If not all elements in ``data_list`` are of the same data type.
    TypeError
        If one or more elements in the list are of an unsupported data type.

    Examples
    --------
    >>> from pyXLMS.transform import get_available_keys
    >>> from pyXLMS import data
    >>> data_list = [
    ...     data.create_crosslink_min("PEPK", 4, "PKEP", 2),
    ...     data.create_crosslink_min("KPEP", 1, "PEKP", 3),
    ... ]
    >>> available_keys = get_available_keys(data_list)
    >>> available_keys["alpha_peptide"]
    True
    >>> available_keys["score"]
    False
    """
    if not assert_data_type_same(data_list):
        raise TypeError("Not all elements of the list have the same data type!")
    if len(data_list) == 0:
        raise ValueError("Provided data does not contain any elements!")
    # available keys
    modifications_a = True
    proteins_a = True
    xl_position_proteins_a = True
    pep_position_proteins_a = True
    score_a = True
    decoy_a = True
    modifications_b = True
    proteins_b = True
    xl_position_proteins_b = True
    pep_position_proteins_b = True
    score_b = True
    decoy_b = True
    score = True
    charge = True
    rt = True
    im_cv = True
    additional_information = True
    # parse available keys
    if isinstance(data_list[0], Crosslink):
        for data in data_list:
            if data["completeness"] != "full" or always_revalidate:
                if data["alpha_proteins"] is None:
                    proteins_a = False
                if data["alpha_proteins_crosslink_positions"] is None:
                    xl_position_proteins_a = False
                if data["alpha_decoy"] is None:
                    decoy_a = False
                if data["beta_proteins"] is None:
                    proteins_b = False
                if data["beta_proteins_crosslink_positions"] is None:
                    xl_position_proteins_b = False
                if data["beta_decoy"] is None:
                    decoy_b = False
                if data["score"] is None:
                    score = False
                if data["additional_information"] is None:
                    additional_information = False
        return {
            "data_type": True,
            "completeness": True,
            "alpha_peptide": True,
            "alpha_peptide_crosslink_position": True,
            "alpha_proteins": proteins_a,
            "alpha_proteins_crosslink_positions": xl_position_proteins_a,
            "alpha_decoy": decoy_a,
            "beta_peptide": True,
            "beta_peptide_crosslink_position": True,
            "beta_proteins": proteins_b,
            "beta_proteins_crosslink_positions": xl_position_proteins_b,
            "beta_decoy": decoy_b,
            "crosslink_type": True,
            "score": score,
            "additional_information": additional_information,
        }
    if isinstance(data_list[0], CrosslinkSpectrumMatch):
        for data in data_list:
            if data["completeness"] != "full" or always_revalidate:
                if data["alpha_modifications"] is None:
                    modifications_a = False
                if data["alpha_proteins"] is None:
                    proteins_a = False
                if data["alpha_proteins_crosslink_positions"] is None:
                    xl_position_proteins_a = False
                if data["alpha_proteins_peptide_positions"] is None:
                    pep_position_proteins_a = False
                if data["alpha_score"] is None:
                    score_a = False
                if data["alpha_decoy"] is None:
                    decoy_a = False
                if data["beta_modifications"] is None:
                    modifications_b = False
                if data["beta_proteins"] is None:
                    proteins_b = False
                if data["beta_proteins_crosslink_positions"] is None:
                    xl_position_proteins_b = False
                if data["beta_proteins_peptide_positions"] is None:
                    pep_position_proteins_b = False
                if data["beta_score"] is None:
                    score_b = False
                if data["beta_decoy"] is None:
                    decoy_b = False
                if data["score"] is None:
                    score = False
                if data["charge"] is None:
                    charge = False
                if data["retention_time"] is None:
                    rt = False
                if data["ion_mobility"] is None:
                    im_cv = False
                if data["additional_information"] is None:
                    additional_information = False
        return {
            "data_type": True,
            "completeness": True,
            "alpha_peptide": True,
            "alpha_modifications": modifications_a,
            "alpha_peptide_crosslink_position": True,
            "alpha_proteins": proteins_a,
            "alpha_proteins_crosslink_positions": xl_position_proteins_a,
            "alpha_proteins_peptide_positions": pep_position_proteins_a,
            "alpha_score": score_a,
            "alpha_decoy": decoy_a,
            "beta_peptide": True,
            "beta_modifications": modifications_b,
            "beta_peptide_crosslink_position": True,
            "beta_proteins": proteins_b,
            "beta_proteins_crosslink_positions": xl_position_proteins_b,
            "beta_proteins_peptide_positions": pep_position_proteins_b,
            "beta_score": score_b,
            "beta_decoy": decoy_b,
            "crosslink_type": True,
            "score": score,
            "spectrum_file": True,
            "scan_nr": True,
            "charge": charge,
            "retention_time": rt,
            "ion_mobility": im_cv,
            "additional_information": additional_information,
        }
    raise TypeError(
        f"Unknown data type {type(data_list[0])}. Data type must be Crosslink or CrosslinkSpectrumMatch!"
    )
    return {"err": True}


def check_available_keys(
    required_keys: List[
        Literal[
            "data_type",
            "completeness",
            "alpha_peptide",
            "alpha_modifications",
            "alpha_peptide_crosslink_position",
            "alpha_proteins",
            "alpha_proteins_crosslink_positions",
            "alpha_proteins_peptide_positions",
            "alpha_score",
            "alpha_decoy",
            "beta_peptide",
            "beta_modifications",
            "beta_peptide_crosslink_position",
            "beta_proteins",
            "beta_proteins_crosslink_positions",
            "beta_proteins_peptide_positions",
            "beta_score",
            "beta_decoy",
            "crosslink_type",
            "score",
            "spectrum_file",
            "scan_nr",
            "charge",
            "retention_time",
            "ion_mobility",
            "additional_information",
        ]
    ],
    data_list: List[CrosslinkSpectrumMatch] | List[Crosslink],
    always_revalidate: bool = True,
) -> bool:
    r"""Checks if all required keys are available in a list of crosslinks or crosslink-spectrum-matches.

    Parameters
    ----------
    required_keys : list of keys
        A list of valid Crosslink or CrosslinkSpectrumMatch keys/attributes to be checked.
    data_list : list of CrosslinkSpectrumMatch, or list of Crosslink
        A list of crosslinks or crosslink-spectrum-matches.
    always_revalidate : bool, default = True
        If ``True`` (default) the assigned ``completeness`` will be ignored and all data fields
        are re-checked. This is safer especially when data has been modified post reading.

    Returns
    -------
    bool
        True if all items in the data list have the required keys and the keys are not None.

    Raises
    ------
    ValueError
        If one of the keys is not available or None in any of the items in the data list.

    Examples
    --------
    >>> from pyXLMS.transform import check_available_keys
    >>> from pyXLMS import data
    >>> data_list = [
    ...     data.create_crosslink_min("PEPK", 4, "PKEP", 2),
    ...     data.create_crosslink_min("KPEP", 1, "PEKP", 3),
    ... ]
    >>> check_available_keys(["alpha_peptide"], data_list)
    True
    >>> check_available_keys(["score"], data_list)
    ValueError: Attribute 'score' is missing in at least one element but is required!
    """
    available_keys = get_available_keys(data_list, always_revalidate)
    for key in required_keys:
        if key not in available_keys or not available_keys[key]:
            raise ValueError(
                f"Attribute '{key}' is missing in at least one element but is required!"
            )
    return True


def display(
    data: CrosslinkSpectrumMatch | Crosslink | ParserResult,
    show_additional_information: bool = False,
    return_str: bool = False,
) -> None | str:
    r"""Pretty prints a crosslink-spectrum-match or crosslink or parser_result.

    Parameters
    ----------
    data : CrosslinkSpectrumMatch, Crosslink, or ParserResult
        A crosslink-spectrum-match or crosslink or parser_result to display.
    show_additional_information : bool, default = False
        Also display data in the ``additional_information``.
    return_str : bool, default = False
        If the display string should be returned.

    Returns
    -------
    None, or str
        The display string of the crosslink-spectrum-match, crosslink, or parser_result
        if ``return_str = True`` otherwise None.

    Raises
    ------
    TypeError
        If data is not a crosslink-spectrum-match, crosslink, or parser_result.

    Examples
    --------
    >>> from pyXLMS import parser
    >>> from pyXLMS import transform
    >>> pr = parser.read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> transform.display(pr)
    Data Type:                            parser_result
    Completeness:                         full
    Identifying Search Engine:            MS Annika
    Number of Crosslink-Spectrum-Matches: 826
    Number of Crosslinks:                 300

    >>> from pyXLMS import parser
    >>> from pyXLMS import transform
    >>> pr = parser.read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> csms = pr["crosslink-spectrum-matches"]
    >>> transform.display(csms[0])
    Data Type:                          crosslink-spectrum-match
    Completeness:                       full
    Alpha Peptide:                      GQKNSR
    Alpha Modifications:                {3: ('DSS', 138.06808)}
    Alpha Peptide Crosslink Position:   3
    Alpha Proteins:                     ['Cas9']
    Alpha Proteins Crosslink Positions: [779]
    Alpha Proteins Peptide Positions:   [777]
    Alpha Peptide Score:                119.82548987540834
    Alpha Decoy:                        False
    Beta Peptide:                       GQKNSR
    Beta Modifications:                 {3: ('DSS', 138.06808)}
    Beta Peptide Crosslink Position:    3
    Beta Proteins:                      ['Cas9']
    Beta Proteins Crosslink Positions:  [779]
    Beta Proteins Peptide Positions:    [777]
    Beta Peptide Score:                 119.82547820493929
    Beta Decoy:                         False
    Crosslink Type:                     intra
    CSM Score:                          119.82547820493929
    Spectrum File:                      XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw
    Scan Number:                        2257
    Precursor Charge:                   3
    Retention Time:                     733.1895599999999
    Ion Mobility/FAIMS CV:              0.0

    >>> from pyXLMS import parser
    >>> from pyXLMS import transform
    >>> pr = parser.read(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> xls = pr["crosslinks"]
    >>> transform.display(xls[0])
    Data Type:                          crosslink
    Completeness:                       full
    Alpha Peptide:                      GQKNSR
    Alpha Peptide Crosslink Position:   3
    Alpha Proteins:                     ['Cas9']
    Alpha Proteins Crosslink Positions: [779]
    Alpha Decoy:                        False
    Beta Peptide:                       GQKNSR
    Beta Peptide Crosslink Position:    3
    Beta Proteins:                      ['Cas9']
    Beta Proteins Crosslink Positions:  [779]
    Beta Decoy:                         False
    Crosslink Type:                     intra
    Crosslink Score:                    119.82547820493929
    """
    _ok = check_input_multi(
        data, "data", [CrosslinkSpectrumMatch, Crosslink, ParserResult]
    )
    _ok = check_input(show_additional_information, "show_additional_information", bool)
    _ok = check_input(return_str, "return_str", bool)
    display: str = ""
    if isinstance(data, Crosslink):
        display += f"Data Type:                          {data['data_type']}\n"
        display += f"Completeness:                       {data['completeness']}\n"
        display += f"Alpha Peptide:                      {data['alpha_peptide']}\n"
        display += f"Alpha Peptide Crosslink Position:   {data['alpha_peptide_crosslink_position']}\n"
        display += f"Alpha Proteins:                     {data['alpha_proteins']}\n"
        display += f"Alpha Proteins Crosslink Positions: {data['alpha_proteins_crosslink_positions']}\n"
        display += f"Alpha Decoy:                        {data['alpha_decoy']}\n"
        display += f"Beta Peptide:                       {data['beta_peptide']}\n"
        display += f"Beta Peptide Crosslink Position:    {data['beta_peptide_crosslink_position']}\n"
        display += f"Beta Proteins:                      {data['beta_proteins']}\n"
        display += f"Beta Proteins Crosslink Positions:  {data['beta_proteins_crosslink_positions']}\n"
        display += f"Beta Decoy:                         {data['beta_decoy']}\n"
        display += f"Crosslink Type:                     {data['crosslink_type']}\n"
        display += f"Crosslink Score:                    {data['score']}\n"
        if show_additional_information:
            display += f"Additional Information:             {data['additional_information']}\n"
        display = display.strip()
        print(display)
        if return_str:
            return display
        return
    if isinstance(data, CrosslinkSpectrumMatch):
        display += f"Data Type:                          {data['data_type']}\n"
        display += f"Completeness:                       {data['completeness']}\n"
        display += f"Alpha Peptide:                      {data['alpha_peptide']}\n"
        display += f"Alpha Modifications:                {data['alpha_modifications']}\n"  # fmt: skip
        display += f"Alpha Peptide Crosslink Position:   {data['alpha_peptide_crosslink_position']}\n"
        display += f"Alpha Proteins:                     {data['alpha_proteins']}\n"
        display += f"Alpha Proteins Crosslink Positions: {data['alpha_proteins_crosslink_positions']}\n"
        display += f"Alpha Proteins Peptide Positions:   {data['alpha_proteins_peptide_positions']}\n"
        display += f"Alpha Peptide Score:                {data['alpha_score']}\n"
        display += f"Alpha Decoy:                        {data['alpha_decoy']}\n"
        display += f"Beta Peptide:                       {data['beta_peptide']}\n"
        display += f"Beta Modifications:                 {data['beta_modifications']}\n"
        display += f"Beta Peptide Crosslink Position:    {data['beta_peptide_crosslink_position']}\n"
        display += f"Beta Proteins:                      {data['beta_proteins']}\n"
        display += f"Beta Proteins Crosslink Positions:  {data['beta_proteins_crosslink_positions']}\n"
        display += f"Beta Proteins Peptide Positions:    {data['beta_proteins_peptide_positions']}\n"
        display += f"Beta Peptide Score:                 {data['beta_score']}\n"
        display += f"Beta Decoy:                         {data['beta_decoy']}\n"
        display += f"Crosslink Type:                     {data['crosslink_type']}\n"
        display += f"CSM Score:                          {data['score']}\n"
        display += f"Spectrum File:                      {data['spectrum_file']}\n"
        display += f"Scan Number:                        {data['scan_nr']}\n"
        display += f"Precursor Charge:                   {data['charge']}\n"
        display += f"Retention Time:                     {data['retention_time']}\n"
        display += f"Ion Mobility/FAIMS CV:              {data['ion_mobility']}\n"
        if show_additional_information:
            display += f"Additional Information:             {data['additional_information']}\n"
        display = display.strip()
        print(display)
        if return_str:
            return display
        return
    if isinstance(data, ParserResult):
        csms = data["crosslink-spectrum-matches"]
        xls = data["crosslinks"]
        display += f"Data Type:                            {data['data_type']}\n"
        display += f"Completeness:                         {data['completeness']}\n"
        display += f"Identifying Search Engine:            {data['search_engine']}\n"
        display += f"Number of Crosslink-Spectrum-Matches: {len(csms) if csms is not None else None}\n"
        display += f"Number of Crosslinks:                 {len(xls) if xls is not None else None}\n"
        display = display.strip()
        print(display)
        if return_str:
            return display
        return
    raise TypeError(
        f"Unknown data type {type(data)}. Data type must be CrosslinkSpectrumMatch, Crosslink, or ParserResult!"
    )
    return
