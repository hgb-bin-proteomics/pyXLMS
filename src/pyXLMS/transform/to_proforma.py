#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd

from ..data import check_input
from ..data import check_input_multi

from typing import Optional
from typing import Dict
from typing import List
from typing import Tuple
from typing import Any


def __get_modified_peptide(
    sequence: str,
    modifications: Optional[Dict[int, Tuple[str, float]]],
    crosslink_position: int,
    crosslinker: Optional[str | float],
) -> str:
    r"""Returns the Proforma string for a single peptide.

    Parameters
    ----------
    sequence : str
        The unmodified peptide sequence.
    modifications : dict of int, tuple of str and float
        The pyXLMS specific modifications object. See ``data.create_csm()`` for reference.
    crosslink_position : int
        Crosslink position in the peptide sequence (1-based).
    crosslinker : str, or float, or None
        Optional name or mass of the crosslink reagent. If the name is given, it should be a valid
        name from XLMOD.

    Returns
    -------
    str
        The Proforma string of the peptidoform.

    Notes
    -----
    - This function should not be called directly, it is called from ``__to_proforma_csm()`` and ``__to_proforma_xl``.
    - Modifications with unknown mass are skipped.
    - If no modifications are given, only the crosslink modification will be encoded in the Proforma.
    - If no modifications are given and no crosslinker is given, the unmodified peptide Proforma will be returned.
    """
    if isinstance(crosslinker, float):
        crosslinker = f"+{crosslinker}" if crosslinker > 0.0 else f"{crosslinker}"
    pep_len = len(sequence)
    if modifications is not None:
        new_modifications = dict()
        for pos, mod in modifications.items():
            if not pd.isna(mod[1]):
                new_modifications[pos] = (
                    mod[0],
                    f"+{mod[1]}" if mod[1] > 0.0 else f"{mod[1]}",
                )
        if crosslink_position not in new_modifications and crosslinker is not None:
            new_modifications[crosslink_position] = ("", f"{crosslinker}")
        for pos in sorted(new_modifications.keys(), reverse=True):
            if pos == 0:
                sequence = f"[{new_modifications[pos][1]}]-" + sequence
            elif pos == pep_len + 1:
                sequence = sequence + f"-[{new_modifications[pos][1]}]"
            else:
                sequence = (
                    sequence[:pos] + f"[{new_modifications[pos][1]}]" + sequence[pos:]
                )
        return sequence
    if crosslinker is not None:
        sequence = (
            sequence[:crosslink_position]
            + f"[{crosslinker}]"
            + sequence[crosslink_position:]
        )
        return sequence
    return sequence


def __to_proforma_csm(csm: Dict[str, Any], crosslinker: Optional[str | float]) -> str:
    r"""Returns the Proforma string for a single crosslink-spectrum-match.

    Parameters
    ----------
    csm : dict of str, any
        A pyXLMS crosslink-spectrum-match object. See ``data.create_csm()``.
    crosslinker : str, or float, or None
        Optional name or mass of the crosslink reagent. If the name is given, it should be a valid
        name from XLMOD.

    Returns
    -------
    str
        The Proforma string of the crosslink-spectrum-match.

    Notes
    -----
    - This function should not be called directly, it is called from ``to_proforma()``.
    - Modifications with unknown mass are skipped.
    - If no modifications are given, only the crosslink modification will be encoded in the Proforma.
    - If no modifications are given and no crosslinker is given, the unmodified peptide Proforma will be returned.
    """
    peptide_a = __get_modified_peptide(
        csm["alpha_peptide"],
        csm["alpha_modifications"],
        csm["alpha_peptide_crosslink_position"],
        crosslinker,
    )
    peptide_b = __get_modified_peptide(
        csm["beta_peptide"],
        csm["beta_modifications"],
        csm["beta_peptide_crosslink_position"],
        crosslinker,
    )
    if csm["charge"] is not None:
        return f"{peptide_a}//{peptide_b}/{csm['charge']}"
    return f"{peptide_a}//{peptide_b}"


def __to_proforma_xl(xl: Dict[str, Any], crosslinker: Optional[str | float]) -> str:
    r"""Returns the Proforma string for a single crosslink.

    Parameters
    ----------
    xl : dict of str, any
        A pyXLMS crosslink object. See ``data.create_crosslink()``.
    crosslinker : str, or float, or None
        Optional name or mass of the crosslink reagent. If the name is given, it should be a valid
        name from XLMOD.

    Returns
    -------
    str
        The Proforma string of the crosslink.

    Notes
    -----
    - This function should not be called directly, it is called from ``to_proforma()``.
    - Modifications with unknown mass are skipped.
    - If no modifications are given, only the crosslink modification will be encoded in the Proforma.
    - If no modifications are given and no crosslinker is given, the unmodified peptide Proforma will be returned.
    """
    peptide_a = __get_modified_peptide(
        xl["alpha_peptide"], None, xl["alpha_peptide_crosslink_position"], crosslinker
    )
    peptide_b = __get_modified_peptide(
        xl["beta_peptide"], None, xl["beta_peptide_crosslink_position"], crosslinker
    )
    return f"{peptide_a}//{peptide_b}"


def to_proforma(
    data: Dict[str, Any] | List[Dict[str, Any]],
    crosslinker: Optional[str | float] = None,
) -> str | List[str]:
    r"""Returns the Proforma string for a single crosslink or crosslink-spectrum-match, or for
    a list of crosslinks or crosslink-spectrum-matches.

    Parameters
    ----------
    data : dict of str, any, or list of dict of str, any
        A pyXLMS crosslink object, e.g. see ``data.create_crosslink()``. Or a pyXLMS crosslink-spectrum-match
        object, e.g. see ``data.create_csm()``. Alternatively, a list of crosslinks or crosslink-spectrum-matches
        can also be provided.
    crosslinker : str, or float, or None, default = None
        Optional name or mass of the crosslink reagent. If the name is given, it should be a valid
        name from XLMOD. If the crosslink modification is contained in the crosslink-spectrum-match object
        this parameter has no effect.

    Returns
    -------
    str
        The Proforma string of the crosslink or crosslink-spectrum-match. If a list was provided
        a list containing all Proforma strings is returned.

    Raises
    ------
    TypeError
        If an unsupported data type is provided.

    Notes
    -----
    - Modifications with unknown mass are skipped.
    - If no modifications are given, only the crosslink modification will be encoded in the Proforma.
    - If no modifications are given and no crosslinker is given, the unmodified peptide Proforma will be returned.

    Examples
    --------
    >>> from pyXLMS.data import create_crosslink_min
    >>> from pyXLMS.transform import to_proforma
    >>> xl = create_crosslink_min("PEPKTIDE", 4, "KPEPTIDE", 1)
    >>> to_proforma(xl)
    'KPEPTIDE//PEPKTIDE'

    >>> from pyXLMS.data import create_crosslink_min
    >>> from pyXLMS.transform import to_proforma
    >>> xl = create_crosslink_min("PEPKTIDE", 4, "KPEPTIDE", 1)
    >>> to_proforma(xl, crosslinker="Xlink:DSSO")
    'K[Xlink:DSSO]PEPTIDE//PEPK[Xlink:DSSO]TIDE'

    >>> from pyXLMS.data import create_csm_min
    >>> from pyXLMS.transform import to_proforma
    >>> csm = create_csm_min("PEPKTIDE", 4, "KPEPTIDE", 1, "RUN_1", 1)
    >>> to_proforma(csm)
    'KPEPTIDE//PEPKTIDE'

    >>> from pyXLMS.data import create_csm_min
    >>> from pyXLMS.transform import to_proforma
    >>> csm = create_csm_min("PEPKTIDE", 4, "KPEPTIDE", 1, "RUN_1", 1)
    >>> to_proforma(csm, crosslinker="Xlink:DSSO")
    'K[Xlink:DSSO]PEPTIDE//PEPK[Xlink:DSSO]TIDE'

    >>> from pyXLMS.data import create_csm_min
    >>> from pyXLMS.transform import to_proforma
    >>> csm = create_csm_min("PEPKTIDE", 4, "KPMEPTIDE", 1, "RUN_1", 1, modifications_b={3:("Oxidation", 15.994915)})
    >>> to_proforma(csm, crosslinker="Xlink:DSSO")
    'K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE'

    >>> from pyXLMS.data import create_csm_min
    >>> from pyXLMS.transform import to_proforma
    >>> csm = create_csm_min("PEPKTIDE", 4, "KPMEPTIDE", 1, "RUN_1", 1, modifications_b={3:("Oxidation", 15.994915)}, charge=3)
    >>> to_proforma(csm, crosslinker="Xlink:DSSO")
    'K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE/3'

    >>> from pyXLMS.data import create_csm_min
    >>> from pyXLMS.transform import to_proforma
    >>> csm = create_csm_min("PEPKTIDE", 4, "KPMEPTIDE", 1, "RUN_1", 1, modifications_a={4:("DSSO", 158.00376)}, modifications_b={1:("DSSO", 158.00376), 3:("Oxidation", 15.994915)}, charge=3)
    >>> to_proforma(csm)
    'K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3'

    >>> from pyXLMS.data import create_csm_min
    >>> from pyXLMS.transform import to_proforma
    >>> csm = create_csm_min("PEPKTIDE", 4, "KPMEPTIDE", 1, "RUN_1", 1, modifications_a={4:("DSSO", 158.00376)}, modifications_b={1:("DSSO", 158.00376), 3:("Oxidation", 15.994915)}, charge=3)
    >>> to_proforma(csm, crosslinker="Xlink:DSSO")
    'K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3'
    """
    _ok = (
        check_input_multi(crosslinker, "crosslinker", [str, float])
        if crosslinker is not None
        else True
    )
    if isinstance(data, list):
        _ok = check_input(data, "data", list, dict)
        return [to_proforma(item, crosslinker) for item in data]  # pyright: ignore[reportReturnType]
    _ok = check_input(data, "data", dict)
    if "data_type" not in data or data["data_type"] not in [
        "crosslink",
        "crosslink-spectrum-match",
    ]:
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a (list of) crosslink or crosslink-spectrum-match!"
        )
    if data["data_type"] == "crosslink":
        return __to_proforma_xl(data, crosslinker)
    return __to_proforma_csm(data, crosslinker)
