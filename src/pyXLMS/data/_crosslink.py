#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import copy
import numpy as np
from pydantic import BaseModel
from pydantic import Field
from pydantic import ConfigDict
from pydantic import computed_field

from ._util import check_input
from ._util import check_indexing
from ._util import __get_modified_peptide as get_modified_peptide

from typing import override
from typing import Annotated
from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


class Crosslink(BaseModel):
    r"""Core data structure representing a single crosslink.

    Crosslinks represent two crosslinked peptides. Crosslinks can be unique peptide
    pairs or unique residue pairs, depending on their grouping.

    Attributes Summary
    ------------------
    Here is a short summary about the crosslink attributes, for more details
    on the specific Pydantic validation requirements please refer to the corresponding attributes
    themselves.

    Required
    ^^^^^^^^
    The following attributes are required:

    alpha_peptide : str
        The unmodified amino acid sequence of the first peptide. Amino acids should be
        in upper case. Modifications should not be included in the sequence.
    alpha_peptide_crosslink_position : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    beta_peptide : str
        The unmodified amino acid sequence of the second peptide. Amino acids should be
        in upper case. Modifications should not be included in the sequence.
    beta_peptide_crosslink_position : int
        The position of the crosslinker in the sequence of the second peptide (1-based).

    Optional
    ^^^^^^^^
    The following attributes are optional:

    alpha_proteins : list of str, or None, default = None
        The accessions of proteins that the first peptide is associated with.
    alpha_proteins_crosslink_positions : list of int, or None, default = None
        Positions of the crosslink in the proteins of the first peptide (1-based). If given the list
        should be of the same length as ``alpha_proteins`` and crosslink position at list index ``i``
        should correspond to the protein at list index ``i`` in ``alpha_proteins``.
    alpha_decoy : bool, or None, default = None
        Whether the first peptide is from the decoy database (``True``) or not (``False``).
    beta_proteins : list of str, or None, default = None
        The accessions of proteins that the second peptide is associated with.
    beta_proteins_crosslink_positions : list of int, or None, default = None
        Positions of the crosslink in the proteins of the second peptide (1-based). If given the list
        should be of the same length as ``beta_proteins`` and crosslink position at list index ``i``
        should correspond to the protein at list index ``i`` in ``beta_proteins``.
    beta_decoy : bool, or None, default = None
        Whether the second peptide is from the decoy database (``True``) or not (``False``).
    score : float, or None, default = None
        Score of the crosslink.
    additional_information : dict of str, any, or None, default = None
        A dictionary with additional information associated with the crosslink.

    Notes
    -----
    Alpha and beta assignment is internally decided by whichever peptide's sequence
    is alphabetically first. If the ``beta_peptide``'s sequence comes alphabetically
    first it will be assigned to ``alpha_peptide`` and the original ``alpha_peptide``
    will be assigned to ``beta_peptide`` (and the same happens for all other corresponding
    alpha and beta values).

    Examples
    --------
    >>> from pyXLMS.data import Crosslink
    >>> xl = Crosslink(
    ...     alpha_peptide="PEKP",
    ...     alpha_peptide_crosslink_position=3,
    ...     beta_peptide="TKIDE",
    ...     beta_peptide_crosslink_position=2,
    ... )
    """

    alpha_peptide: Annotated[
        str,
        Field(
            frozen=True,
            description="The unmodified amino acid sequence of the first peptide.",
        ),
    ]
    r"""
    The unmodified amino acid sequence of the first peptide. Amino acids should be
    in upper case. Modifications should not be included in the sequence.
    """
    alpha_peptide_crosslink_position: Annotated[
        int,
        Field(
            frozen=True,
            description="The position of the crosslinker in the sequence of the first peptide (1-based).",
        ),
    ]
    r"""
    The position of the crosslinker in the sequence of the first peptide (1-based).
    """
    beta_peptide: Annotated[
        str,
        Field(
            frozen=True,
            description="The unmodified amino acid sequence of the second peptide.",
        ),
    ]
    r"""
    The unmodified amino acid sequence of the second peptide. Amino acids should be
    in upper case. Modifications should not be included in the sequence.
    """
    beta_peptide_crosslink_position: Annotated[
        int,
        Field(
            frozen=True,
            description="The position of the crosslinker in the sequence of the second peptide (1-based).",
        ),
    ]
    r"""
    The position of the crosslinker in the sequence of the second peptide (1-based).
    """
    alpha_proteins: Annotated[
        Optional[List[str]],
        Field(
            frozen=True,
            description="The accessions of proteins that the first peptide is associated with.",
        ),
    ] = None
    r"""
    The accessions of proteins that the first peptide is associated with.
    """
    alpha_proteins_crosslink_positions: Annotated[
        Optional[List[int]],
        Field(
            frozen=True,
            description="Positions of the crosslink in the proteins of the first peptide (1-based).",
        ),
    ] = None
    r"""
    Positions of the crosslink in the proteins of the first peptide (1-based). If given the list
    should be of the same length as ``alpha_proteins`` and crosslink position at list index ``i``
    should correspond to the protein at list index ``i`` in ``alpha_proteins``.
    """
    alpha_decoy: Annotated[
        Optional[bool],
        Field(
            frozen=True,
            description="Whether the alpha peptide is from the decoy database or not.",
        ),
    ] = None
    r"""
    Whether the first peptide is from the decoy database (``True``) or not (``False``).
    """
    beta_proteins: Annotated[
        Optional[List[str]],
        Field(
            frozen=True,
            description="The accessions of proteins that the second peptide is associated with.",
        ),
    ] = None
    r"""
    The accessions of proteins that the second peptide is associated with.
    """
    beta_proteins_crosslink_positions: Annotated[
        Optional[List[int]],
        Field(
            frozen=True,
            description="Positions of the crosslink in the proteins of the second peptide (1-based).",
        ),
    ] = None
    r"""
    Positions of the crosslink in the proteins of the second peptide (1-based). If given the list
    should be of the same length as ``beta_proteins`` and crosslink position at list index ``i``
    should correspond to the protein at list index ``i`` in ``beta_proteins``.
    """
    beta_decoy: Annotated[
        Optional[bool],
        Field(
            frozen=True,
            description="Whether the beta peptide is from the decoy database or not.",
        ),
    ] = None
    r"""
    Whether the second peptide is from the decoy database (``True``) or not (``False``).
    """
    score: Annotated[
        Optional[float], Field(frozen=True, description="Score of the crosslink.")
    ] = None
    r"""
    Score of the crosslink.
    """
    additional_information: Annotated[
        Optional[Dict[str, Any]],
        Field(
            frozen=False,
            description="A dictionary with additional information associated with the crosslink.",
        ),
    ] = None
    r"""
    A dictionary with additional information associated with the crosslink.
    """
    model_config = ConfigDict(
        validate_assignment=True, strict=True, str_strip_whitespace=True
    )
    r"""
    Pydantic configuration for the underlying validation model.
    """

    @computed_field(description="Data type of the object.")
    @property
    def data_type(self) -> Literal["crosslink"]:
        r"""
        Data type of the object.
        """
        return "crosslink"

    @computed_field(description="Completeness of the crosslink.")
    @property
    def completeness(self) -> Literal["full", "partial"]:
        r"""
        Completeness of the crosslink, e.g. ``"full"`` if all attributes
        are not ``None`` and else ``"partial"``.
        """
        full = all(
            [
                self.alpha_proteins is not None,
                self.alpha_proteins_crosslink_positions is not None,
                self.alpha_decoy is not None,
                self.beta_proteins is not None,
                self.beta_proteins_crosslink_positions is not None,
                self.beta_decoy is not None,
                self.score is not None,
            ]
        )
        return "full" if full else "partial"

    @computed_field(description="Link type of the crosslink.")
    @property
    def crosslink_type(self) -> Literal["intra", "inter"]:
        r"""
        Link type of the crosslink, e.g. ``"intra"`` if the proteins in
        ``alpha_proteins`` and ``beta_proteins`` overlap, otherwise ``"inter"``.
        """
        a_prot = set(
            [str(protein).strip() for protein in self.alpha_proteins]
            if self.alpha_proteins is not None
            else []
        )
        b_prot = set(
            [str(protein).strip() for protein in self.beta_proteins]
            if self.beta_proteins is not None
            else []
        )
        return "intra" if len(a_prot.intersection(b_prot)) > 0 else "inter"

    @override
    def model_post_init(self, context: Any = None) -> None:
        r"""
        Performs extra validation and post init functions.

        Notes
        -----
        Alpha and beta assignment is internally decided by whichever peptide's sequence
        is alphabetically first. If the ``beta_peptide``'s sequence comes alphabetically
        first it will be assigned to ``alpha_peptide`` and the original ``alpha_peptide``
        will be assigned to ``beta_peptide`` (and the same happens for all other corresponding
        alpha and beta values).

        Warnings
        --------
        This method should not be called manually!
        """
        # extra validation
        if (
            self.alpha_proteins is not None
            and self.alpha_proteins_crosslink_positions is not None
        ):
            if len(self.alpha_proteins) != len(self.alpha_proteins_crosslink_positions):
                raise ValueError(
                    "Crosslink position has to be given for every protein! Length of alpha_proteins and alpha_proteins_crosslink_positions has to match!"
                )
        if (
            self.beta_proteins is not None
            and self.beta_proteins_crosslink_positions is not None
        ):
            if len(self.beta_proteins) != len(self.beta_proteins_crosslink_positions):
                raise ValueError(
                    "Crosslink position has to be given for every protein! Length of beta_proteins and beta_proteins_crosslink_positions has to match!"
                )
        _ok = check_indexing(self.alpha_peptide_crosslink_position)
        _ok = check_indexing(self.beta_peptide_crosslink_position)
        _ok = (
            check_indexing(self.alpha_proteins_crosslink_positions)
            if self.alpha_proteins_crosslink_positions is not None
            else True
        )
        _ok = (
            check_indexing(self.beta_proteins_crosslink_positions)
            if self.beta_proteins_crosslink_positions is not None
            else True
        )
        ## processing
        key_a = f"{self.alpha_peptide.strip()}{self.alpha_peptide_crosslink_position}"
        key_b = f"{self.beta_peptide.strip()}{self.beta_peptide_crosslink_position}"
        # if homomeric crosslink
        if key_a == key_b:
            key_a += "_0"
            key_b += "_1"
        crosslink = {
            key_a: {
                "peptide": self.alpha_peptide.strip(),
                "xl_position_peptide": self.alpha_peptide_crosslink_position,
                "proteins": copy.deepcopy(self.alpha_proteins),
                "xl_position_proteins": copy.deepcopy(
                    self.alpha_proteins_crosslink_positions
                ),
                "decoy": self.alpha_decoy,
            },
            key_b: {
                "peptide": self.beta_peptide.strip(),
                "xl_position_peptide": self.beta_peptide_crosslink_position,
                "proteins": copy.deepcopy(self.beta_proteins),
                "xl_position_proteins": copy.deepcopy(
                    self.beta_proteins_crosslink_positions
                ),
                "decoy": self.beta_decoy,
            },
        }
        keys = sorted(list(crosslink.keys()))
        alpha_proteins_clean = (
            [str(protein).strip() for protein in crosslink[keys[0]]["proteins"]]  # ty: ignore[not-iterable]
            if crosslink[keys[0]]["proteins"] is not None
            else None
        )
        beta_proteins_clean = (
            [str(protein).strip() for protein in crosslink[keys[1]]["proteins"]]  # ty: ignore[not-iterable]
            if crosslink[keys[1]]["proteins"] is not None
            else None
        )
        # re-assign
        self.__dict__["alpha_peptide"] = crosslink[keys[0]]["peptide"]
        self.__dict__["alpha_peptide_crosslink_position"] = crosslink[keys[0]]["xl_position_peptide"]  # fmt: skip
        self.__dict__["alpha_proteins"] = alpha_proteins_clean
        self.__dict__["alpha_proteins_crosslink_positions"] = crosslink[keys[0]]["xl_position_proteins"]  # fmt: skip
        self.__dict__["alpha_decoy"] = crosslink[keys[0]]["decoy"]
        self.__dict__["beta_peptide"] = crosslink[keys[1]]["peptide"]
        self.__dict__["beta_peptide_crosslink_position"] = crosslink[keys[1]]["xl_position_peptide"]  # fmt: skip
        self.__dict__["beta_proteins"] = beta_proteins_clean
        self.__dict__["beta_proteins_crosslink_positions"] = crosslink[keys[1]]["xl_position_proteins"]  # fmt: skip
        self.__dict__["beta_decoy"] = crosslink[keys[1]]["decoy"]
        if self.score is not None:
            if np.isnan(self.score):
                self.__dict__["score"] = None
        return

    def __getitem__(self, key: str) -> Any:
        r"""
        Support for dict-like access.
        """
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(f"'{key}' is not a valid field!")

    def __contains__(self, key: str) -> bool:
        r"""
        Support for ``in`` operator.
        """
        return hasattr(self, key)

    def items(self) -> List[Tuple[str, Any]]:
        r"""
        Support for dict-like read access for backward compatibility.

        Returns
        -------
        list of tuple of str, any
            Returns a list of tuples of attribute name, attribute value.

        Notes
        -----
        This internally just calls ``self.model_dump(mode="python").items()``.
        See `model_dump <https://pydantic.dev/docs/validation/latest/api/pydantic/base_model/#pydantic.BaseModel.model_dump>`_.
        """
        return self.model_dump(mode="python").items()

    def keys(self) -> List[str]:
        r"""
        Support for dict-like read access for backward compatibility.

        Returns
        -------
        list of str
            Returns a list of attribute names.

        Notes
        -----
        This internally just calls ``self.model_dump(mode="python").keys()``.
        See `model_dump <https://pydantic.dev/docs/validation/latest/api/pydantic/base_model/#pydantic.BaseModel.model_dump>`_.
        """
        return self.model_dump(mode="python").keys()

    def values(self) -> List[Any]:
        r"""
        Support for dict-like read access for backward compatibility.

        Returns
        -------
        list of any
            Returns a list of attribute values.

        Notes
        -----
        This internally just calls ``self.model_dump(mode="python").values()``.
        See `model_dump <https://pydantic.dev/docs/validation/latest/api/pydantic/base_model/#pydantic.BaseModel.model_dump>`_.
        """
        return self.model_dump(mode="python").values()

    def copy_with_update(self, update: Dict[str, Any] = {}) -> Crosslink:
        r"""Creates a deep copy of the crosslink with optional attribute updates.

        Parameters
        ----------
        update : dict of str, any, default = empty dict
            Dictionary mapping attribute names (str) to their updated values.
            The default (empty dict) will create a deep copy with the original
            attribute values.

        Returns
        -------
        Crosslink
            New crosslink with optionally updated attributes.

        Examples
        --------
        >>> from pyXLMS.data import Crosslink
        >>> xl = Crosslink(
        ...     alpha_peptide="PEKP",
        ...     alpha_peptide_crosslink_position=3,
        ...     alpha_proteins=["PROT"],
        ...     beta_peptide="PEKP",
        ...     beta_peptide_crosslink_position=3,
        ...     beta_proteins=["PROT"],
        ... )
        >>> xl_copy = xl.copy_with_update(
        ...     update={"additional_information": {"homomeric": True}}
        ... )
        """
        _ok = check_input(update, "update", dict)
        return Crosslink(
            alpha_peptide=self.alpha_peptide
            if "alpha_peptide" not in update
            else update["alpha_peptide"],
            alpha_peptide_crosslink_position=self.alpha_peptide_crosslink_position
            if "alpha_peptide_crosslink_position" not in update
            else update["alpha_peptide_crosslink_position"],
            beta_peptide=self.beta_peptide
            if "beta_peptide" not in update
            else update["beta_peptide"],
            beta_peptide_crosslink_position=self.beta_peptide_crosslink_position
            if "beta_peptide_crosslink_position" not in update
            else update["beta_peptide_crosslink_position"],
            alpha_proteins=copy.deepcopy(self.alpha_proteins)
            if "alpha_proteins" not in update
            else update["alpha_proteins"],
            alpha_proteins_crosslink_positions=copy.deepcopy(
                self.alpha_proteins_crosslink_positions
            )
            if "alpha_proteins_crosslink_positions" not in update
            else update["alpha_proteins_crosslink_positions"],
            alpha_decoy=self.alpha_decoy
            if "alpha_decoy" not in update
            else update["alpha_decoy"],
            beta_proteins=copy.deepcopy(self.beta_proteins)
            if "beta_proteins" not in update
            else update["beta_proteins"],
            beta_proteins_crosslink_positions=copy.deepcopy(
                self.beta_proteins_crosslink_positions
            )
            if "beta_proteins_crosslink_positions" not in update
            else update["beta_proteins_crosslink_positions"],
            beta_decoy=self.beta_decoy
            if "beta_decoy" not in update
            else update["beta_decoy"],
            score=self.score if "score" not in update else update["score"],
            additional_information=copy.deepcopy(self.additional_information)
            if "additional_information" not in update
            else update["additional_information"],
        )

    def display(
        self,
        show_additional_information: bool = False,
        return_str: bool = False,
    ) -> None | str:
        r"""Pretty prints the crosslink.

        Parameters
        ----------
        show_additional_information : bool, default = False
            Also display data in the ``additional_information``.
        return_str : bool, default = False
            If the display string should be returned.

        Returns
        -------
        None, or str
            The display string of the crosslink if ``return_str = True`` otherwise None.

        Examples
        --------
        >>> from pyXLMS import parser
        >>> pr = parser.read(
        ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        ...     engine="MS Annika",
        ...     crosslinker="DSS",
        ... )
        >>> xls = pr["crosslinks"]
        >>> xls[0].display()
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
        _ok = check_input(
            show_additional_information, "show_additional_information", bool
        )
        _ok = check_input(return_str, "return_str", bool)
        display: str = ""
        display += f"Data Type:                          {self.data_type}\n"
        display += f"Completeness:                       {self.completeness}\n"
        display += f"Alpha Peptide:                      {self.alpha_peptide}\n"
        display += f"Alpha Peptide Crosslink Position:   {self.alpha_peptide_crosslink_position}\n"
        display += f"Alpha Proteins:                     {self.alpha_proteins}\n"
        display += f"Alpha Proteins Crosslink Positions: {self.alpha_proteins_crosslink_positions}\n"
        display += f"Alpha Decoy:                        {self.alpha_decoy}\n"
        display += f"Beta Peptide:                       {self.beta_peptide}\n"
        display += f"Beta Peptide Crosslink Position:    {self.beta_peptide_crosslink_position}\n"
        display += f"Beta Proteins:                      {self.beta_proteins}\n"
        display += f"Beta Proteins Crosslink Positions:  {self.beta_proteins_crosslink_positions}\n"
        display += f"Beta Decoy:                         {self.beta_decoy}\n"
        display += f"Crosslink Type:                     {self.crosslink_type}\n"
        display += f"Crosslink Score:                    {self.score}\n"
        if show_additional_information:
            display += (
                f"Additional Information:             {self.additional_information}\n"
            )
        display = display.strip()
        print(display)
        if return_str:
            return display
        return

    def to_proforma(self, crosslinker: Optional[str | float] = None) -> str:
        r"""Returns the Proforma string for the crosslink.

        Parameters
        ----------
        crosslinker : str, or float, or None, default = None
            Optional name or mass of the crosslink reagent. If the name is given, it should be a valid
            name from XLMOD.

        Returns
        -------
        str
            The Proforma string of the crosslink.

        Notes
        -----
        - If no crosslinker is given, the unmodified peptide Proforma will be returned.

        Examples
        --------
        >>> from pyXLMS.data import create_crosslink_min
        >>> xl = create_crosslink_min("PEPKTIDE", 4, "KPEPTIDE", 1)
        >>> xl.to_proforma()
        'KPEPTIDE//PEPKTIDE'

        >>> from pyXLMS.data import create_crosslink_min
        >>> xl = create_crosslink_min("PEPKTIDE", 4, "KPEPTIDE", 1)
        >>> xl.to_proforma(crosslinker="Xlink:DSSO")
        'K[Xlink:DSSO]PEPTIDE//PEPK[Xlink:DSSO]TIDE'
        """
        peptide_a = get_modified_peptide(
            self.alpha_peptide,
            None,
            self.alpha_peptide_crosslink_position,
            crosslinker,
        )
        peptide_b = get_modified_peptide(
            self.beta_peptide, None, self.beta_peptide_crosslink_position, crosslinker
        )
        return f"{peptide_a}//{peptide_b}"


def create_crosslink(
    peptide_a: str,
    xl_position_peptide_a: int,
    proteins_a: Optional[List[str]],
    xl_position_proteins_a: Optional[List[int]],
    decoy_a: Optional[bool],
    peptide_b: str,
    xl_position_peptide_b: int,
    proteins_b: Optional[List[str]],
    xl_position_proteins_b: Optional[List[int]],
    decoy_b: Optional[bool],
    score: Optional[float],
    additional_information: Optional[Dict[str, Any]] = None,
) -> Crosslink:
    r"""Creates a crosslink data structure.

    Contains minimal data necessary for representing a single crosslink. The returned crosslink data structure is a dictionary with keys
    as detailed in the return section.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    proteins_a : list of str, or None
        The accessions of proteins that the first peptide is associated with.
    xl_position_proteins_a : list of int, or None
        Positions of the crosslink in the proteins of the first peptide (1-based).
    decoy_a : bool, or None
        Whether the alpha peptide is from the decoy database or not.
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    proteins_b : list of str, or None
        The accessions of proteins that the second peptide is associated with.
    xl_position_proteins_b : list of int, or None
        Positions of the crosslink in the proteins of the second peptide (1-based).
    decoy_b : bool, or None
        Whether the beta peptide is from the decoy database or not.
    score: float, or None
        Score of the crosslink.
    additional_information: dict with str keys, or None, default = None
        A dictionary with additional information associated with the crosslink.

    Returns
    -------
    dict
        The dictionary representing the crosslink with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_peptide_crosslink_position``,
        ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_decoy``, ``beta_peptide``, ``beta_peptide_crosslink_position``,
        ``beta_proteins``, ``beta_proteins_crosslink_positions``, ``beta_decoy``, ``crosslink_type``, ``score``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.
    ValueError
        If the length of crosslink positions is not equal to the length of proteins.

    Notes
    -----
    The minimum required data for creating a crosslink is:

    - ``peptide_a``: The unmodified amino acid sequence of the first peptide.
    - ``peptide_b``: The unmodified amino acid sequence of the second peptide.
    - ``xl_position_peptide_a``: The position of the crosslinker in the sequence of the first peptide (1-based).
    - ``xl_position_peptide_b``: The position of the crosslinker in the sequence of the second peptide (1-based).

    Examples
    --------
    >>> from pyXLMS.data import create_crosslink
    >>> minimal_crosslink = create_crosslink(
    ...     peptide_a="PEPTIDEA",
    ...     xl_position_peptide_a=1,
    ...     proteins_a=None,
    ...     xl_position_proteins_a=None,
    ...     decoy_a=None,
    ...     peptide_b="PEPTIDEB",
    ...     xl_position_peptide_b=5,
    ...     proteins_b=None,
    ...     xl_position_proteins_b=None,
    ...     decoy_b=None,
    ...     score=None,
    ... )

    >>> from pyXLMS.data import create_crosslink
    >>> crosslink = create_crosslink(
    ...     peptide_a="PEPTIDEA",
    ...     xl_position_peptide_a=1,
    ...     proteins_a=["PROTEINA"],
    ...     xl_position_proteins_a=[1],
    ...     decoy_a=False,
    ...     peptide_b="PEPTIDEB",
    ...     xl_position_peptide_b=5,
    ...     proteins_b=["PROTEINB"],
    ...     xl_position_proteins_b=[3],
    ...     decoy_b=False,
    ...     score=34.5,
    ... )
    """
    return Crosslink(
        alpha_peptide=peptide_a,
        alpha_peptide_crosslink_position=xl_position_peptide_a,
        beta_peptide=peptide_b,
        beta_peptide_crosslink_position=xl_position_peptide_b,
        alpha_proteins=proteins_a,
        alpha_proteins_crosslink_positions=xl_position_proteins_a,
        alpha_decoy=decoy_a,
        beta_proteins=proteins_b,
        beta_proteins_crosslink_positions=xl_position_proteins_b,
        beta_decoy=decoy_b,
        score=score,
        additional_information=additional_information,
    )


def create_crosslink_min(
    peptide_a: str,
    xl_position_peptide_a: int,
    peptide_b: str,
    xl_position_peptide_b: int,
    **kwargs,
) -> Crosslink:
    r"""Creates a crosslink data structure from minimal input.

    Contains minimal data necessary for representing a single crosslink. This is an alias for
    ``data.create_crosslink()``that sets all optional parameters to ``None`` for convenience.
    The returned crosslink data structure is a dictionary with keys as detailed in the return
    section.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    **kwargs
        Any additional parameters will be passed to ``data.create_crosslink()``.

    Returns
    -------
    dict
        The dictionary representing the crosslink with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_peptide_crosslink_position``,
        ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_decoy``, ``beta_peptide``, ``beta_peptide_crosslink_position``,
        ``beta_proteins``, ``beta_proteins_crosslink_positions``, ``beta_decoy``, ``crosslink_type``, ``score``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Notes
    -----
    See also ``data.create_crosslink()``.

    Examples
    --------
    >>> from pyXLMS.data import create_crosslink_min
    >>> minimal_crosslink = create_crosslink_min("PEPTIDEA", 1, "PEPTIDEB", 5)
    """
    return create_crosslink(
        peptide_a=peptide_a,
        xl_position_peptide_a=xl_position_peptide_a,
        proteins_a=kwargs["proteins_a"] if "proteins_a" in kwargs else None,
        xl_position_proteins_a=kwargs["xl_position_proteins_a"]
        if "xl_position_proteins_a" in kwargs
        else None,
        decoy_a=kwargs["decoy_a"] if "decoy_a" in kwargs else None,
        peptide_b=peptide_b,
        xl_position_peptide_b=xl_position_peptide_b,
        proteins_b=kwargs["proteins_b"] if "proteins_b" in kwargs else None,
        xl_position_proteins_b=kwargs["xl_position_proteins_b"]
        if "xl_position_proteins_b" in kwargs
        else None,
        decoy_b=kwargs["decoy_b"] if "decoy_b" in kwargs else None,
        score=kwargs["score"] if "score" in kwargs else None,
        additional_information=kwargs["additional_information"]
        if "additional_information" in kwargs
        else None,
    )
