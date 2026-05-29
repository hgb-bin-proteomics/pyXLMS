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

from ._crosslink import Crosslink
from ._crosslink import create_crosslink
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


class CrosslinkSpectrumMatch(BaseModel):
    r"""Core data structure representing a single crosslink-spectrum-match.

    Crosslink-spectrum-matches associate two crosslinked peptides with a specific
    mass spectrum. They contain spectrum level information additionally to crosslink
    information.

    Attributes Summary
    ------------------
    Here is a short summary about the crosslink-spectrum-match attributes, for more details
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
    spectrum_file : str
        Name of the spectrum file the crosslink-spectrum-match was identified in.
    scan_nr : int
        The corresponding scan number of the crosslink-spectrum-match. If the scan number
        is not available the spectrum index should be provided.

    Optional
    ^^^^^^^^
    The following attributes are optional:

    alpha_modifications : dict of int, tuple of str, float, or None, default = None
        The modifications of the first peptide given as a dictionary that maps peptide position
        (1-based) to modification given as a tuple of modification name and modification delta mass.
        ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications
        should be denoted with position ``len(peptide) + 1``. If the peptide is not modified an empty
        dictionary should be given.
    alpha_proteins : list of str, or None, default = None
        The accessions of proteins that the first peptide is associated with.
    alpha_proteins_crosslink_positions : list of int, or None, default = None
        Positions of the crosslink in the proteins of the first peptide (1-based). If given the list
        should be of the same length as ``alpha_proteins`` and crosslink position at list index ``i``
        should correspond to the protein at list index ``i`` in ``alpha_proteins``.
    alpha_proteins_peptide_positions : list of int, or None, default = None
        Positions of the first peptide in the corresponding proteins (1-based). If given the list
        should be of the same length as ``alpha_proteins`` and peptide position at list index ``i``
        should correspond to the protein at list index ``i`` in ``alpha_proteins``.
    alpha_score : float, or None, default = None
        Identification score of the first peptide.
    alpha_decoy : bool, or None, default = None
        Whether the first peptide is from the decoy database (``True``) or not (``False``).
    beta_modifications : dict of int, tuple of str, float, or None, default = None
        The modifications of the second peptide given as a dictionary that maps peptide position
        (1-based) to modification given as a tuple of modification name and modification delta mass.
        ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications
        should be denoted with position ``len(peptide) + 1``. If the peptide is not modified an empty
        dictionary should be given.
    beta_proteins : list of str, or None, default = None
        The accessions of proteins that the second peptide is associated with.
    beta_proteins_crosslink_positions : list of int, or None, default = None
        Positions of the crosslink in the proteins of the second peptide (1-based). If given the list
        should be of the same length as ``beta_proteins`` and crosslink position at list index ``i``
        should correspond to the protein at list index ``i`` in ``beta_proteins``.
    beta_proteins_peptide_positions : list of int, or None, default = None
        Positions of the second peptide in the corresponding proteins (1-based). If given the list
        should be of the same length as ``beta_proteins`` and peptide position at list index ``i``
        should correspond to the protein at list index ``i`` in ``beta_proteins``.
    beta_score : float, or None, default = None
        Identification score of the second peptide.
    beta_decoy : bool, or None, default = None
        Whether the second peptide is from the decoy database (``True``) or not (``False``).
    score : float, or None, default = None
        Score of the crosslink-spectrum-match.
    charge : int, or None, default = None
        The precursor charge of the corresponding mass spectrum of the crosslink-spectrum-match.
    retention_time : float, or None, default = None
        The retention time of the corresponding mass spectrum of the crosslink-spectrum-match in seconds.
    ion_mobility : float, or None, default = None
        The ion mobility or compensation voltage of the corresponding mass spectrum of the crosslink-spectrum-match.
    additional_information : dict of str, any, or None, default = None
        A dictionary with additional information associated with the crosslink-spectrum-match.

    Notes
    -----
    Alpha and beta assignment is internally decided by whichever peptide's sequence
    is alphabetically first. If the ``beta_peptide``'s sequence comes alphabetically
    first it will be assigned to ``alpha_peptide`` and the original ``alpha_peptide``
    will be assigned to ``beta_peptide`` (and the same happens for all other corresponding
    alpha and beta values).

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
    spectrum_file: Annotated[
        str,
        Field(
            frozen=True,
            description="Name of the spectrum file the crosslink-spectrum-match was identified in.",
        ),
    ]
    r"""
    Name of the spectrum file the crosslink-spectrum-match was identified in.
    """
    scan_nr: Annotated[
        int,
        Field(
            frozen=True,
            description="The corresponding scan number of the crosslink-spectrum-match.",
        ),
    ]
    r"""
    The corresponding scan number of the crosslink-spectrum-match. If the scan number
    is not available the spectrum index should be provided.
    """
    alpha_modifications: Annotated[
        Optional[Dict[int, Tuple[str, float]]],
        Field(frozen=True, description="The modifications of the first peptide."),
    ] = None
    r"""
    The modifications of the first peptide given as a dictionary that maps peptide position
    (1-based) to modification given as a tuple of modification name and modification delta mass.
    ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications
    should be denoted with position ``len(peptide) + 1``. If the peptide is not modified an empty
    dictionary should be given.
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
    alpha_proteins_peptide_positions: Annotated[
        Optional[List[int]],
        Field(
            frozen=True,
            description="Positions of the first peptide in the corresponding proteins (1-based).",
        ),
    ] = None
    r"""
    Positions of the first peptide in the corresponding proteins (1-based). If given the list
    should be of the same length as ``alpha_proteins`` and peptide position at list index ``i``
    should correspond to the protein at list index ``i`` in ``alpha_proteins``.
    """
    alpha_score: Annotated[
        Optional[float],
        Field(frozen=True, description="Identification score of the first peptide."),
    ] = None
    r"""
    Identification score of the first peptide.
    """
    alpha_decoy: Annotated[
        Optional[bool],
        Field(
            frozen=True,
            description="Whether the first peptide is from the decoy database or not.",
        ),
    ] = None
    r"""
    Whether the first peptide is from the decoy database (``True``) or not (``False``).
    """
    beta_modifications: Annotated[
        Optional[Dict[int, Tuple[str, float]]],
        Field(frozen=True, description="The modifications of the second peptide."),
    ] = None
    r"""
    The modifications of the second peptide given as a dictionary that maps peptide position
    (1-based) to modification given as a tuple of modification name and modification delta mass.
    ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications
    should be denoted with position ``len(peptide) + 1``. If the peptide is not modified an empty
    dictionary should be given.
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
    beta_proteins_peptide_positions: Annotated[
        Optional[List[int]],
        Field(
            frozen=True,
            description="Positions of the second peptide in the corresponding proteins (1-based).",
        ),
    ] = None
    r"""
    Positions of the second peptide in the corresponding proteins (1-based). If given the list
    should be of the same length as ``beta_proteins`` and peptide position at list index ``i``
    should correspond to the protein at list index ``i`` in ``beta_proteins``.
    """
    beta_score: Annotated[
        Optional[float],
        Field(frozen=True, description="Identification score of the second peptide."),
    ] = None
    r"""
    Identification score of the second peptide.
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
        Optional[float],
        Field(frozen=True, description="Score of the crosslink-spectrum-match."),
    ] = None
    r"""
    Score of the crosslink-spectrum-match.
    """
    charge: Annotated[
        Optional[int],
        Field(
            frozen=True,
            description="The precursor charge of the corresponding mass spectrum of the crosslink-spectrum-match.",
        ),
    ] = None
    r"""
    The precursor charge of the corresponding mass spectrum of the crosslink-spectrum-match.
    """
    retention_time: Annotated[
        Optional[float],
        Field(
            frozen=True,
            description="The retention time of the corresponding mass spectrum of the crosslink-spectrum-match in seconds.",
        ),
    ] = None
    r"""
    The retention time of the corresponding mass spectrum of the crosslink-spectrum-match in seconds.
    """
    ion_mobility: Annotated[
        Optional[float],
        Field(
            frozen=True,
            description="The ion mobility or compensation voltage of the corresponding mass spectrum of the crosslink-spectrum-match.",
        ),
    ] = None
    r"""
    The ion mobility or compensation voltage of the corresponding mass spectrum of the crosslink-spectrum-match.
    """
    additional_information: Annotated[
        Optional[Dict[str, Any]],
        Field(
            frozen=False,
            description="A dictionary with additional information associated with the crosslink-spectrum-match.",
        ),
    ] = None
    r"""
    A dictionary with additional information associated with the crosslink-spectrum-match.
    """
    model_config = ConfigDict(
        validate_assignment=True, strict=True, str_strip_whitespace=True
    )
    r"""
    Pydantic configuration for the underlying validation model.
    """

    @computed_field(description="Data type of the object.")
    @property
    def data_type(self) -> Literal["crosslink-spectrum-match"]:
        r"""
        Data type of the object.
        """
        return "crosslink-spectrum-match"

    @computed_field(description="Completeness of the crosslink-spectrum-match.")
    @property
    def completeness(self) -> Literal["full", "partial"]:
        r"""
        Completeness of the crosslink-spectrum-match, e.g. ``"full"`` if all attributes
        are not ``None`` and else ``"partial"``.
        """
        full = all(
            [
                self.alpha_modifications is not None,
                self.alpha_proteins is not None,
                self.alpha_proteins_crosslink_positions is not None,
                self.alpha_proteins_peptide_positions is not None,
                self.alpha_score is not None,
                self.alpha_decoy is not None,
                self.beta_modifications is not None,
                self.beta_proteins is not None,
                self.beta_proteins_crosslink_positions is not None,
                self.beta_proteins_peptide_positions is not None,
                self.beta_score is not None,
                self.beta_decoy is not None,
                self.score is not None,
                self.charge is not None,
                self.retention_time is not None,
                self.ion_mobility is not None,
            ]
        )
        return "full" if full else "partial"

    @computed_field(description="Link type of the crosslink-spectrum-match.")
    @property
    def crosslink_type(self) -> Literal["intra", "inter"]:
        r"""
        Link type of the crosslink-spectrum-match, e.g. ``"intra"`` if the proteins in
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
        if (
            self.alpha_proteins is not None
            and self.alpha_proteins_peptide_positions is not None
        ):
            if len(self.alpha_proteins) != len(self.alpha_proteins_peptide_positions):
                raise ValueError(
                    "Peptide position has to be given for every protein! Length of alpha_proteins and alpha_proteins_peptide_positions has to match!"
                )
        if (
            self.beta_proteins is not None
            and self.beta_proteins_peptide_positions is not None
        ):
            if len(self.beta_proteins) != len(self.beta_proteins_peptide_positions):
                raise ValueError(
                    "Peptide position has to be given for every protein! Length of beta_proteins and beta_proteins_peptide_positions has to match!"
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
        _ok = (
            check_indexing(self.alpha_proteins_peptide_positions)
            if self.alpha_proteins_peptide_positions is not None
            else True
        )
        _ok = (
            check_indexing(self.beta_proteins_peptide_positions)
            if self.beta_proteins_peptide_positions is not None
            else True
        )
        ## validity
        if (
            self.alpha_proteins_crosslink_positions is not None
            and self.alpha_proteins_peptide_positions is not None
        ):
            for i in range(len(self.alpha_proteins_crosslink_positions)):
                if (
                    self.alpha_proteins_crosslink_positions[i]
                    - self.alpha_proteins_peptide_positions[i]
                    + 1
                    != self.alpha_peptide_crosslink_position
                ):
                    _ok = check_indexing(0)
        if (
            self.beta_proteins_crosslink_positions is not None
            and self.beta_proteins_peptide_positions is not None
        ):
            for i in range(len(self.beta_proteins_crosslink_positions)):
                if (
                    self.beta_proteins_crosslink_positions[i]
                    - self.beta_proteins_peptide_positions[i]
                    + 1
                    != self.beta_peptide_crosslink_position
                ):
                    _ok = check_indexing(0)
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
                "modifications": copy.deepcopy(
                    {
                        int(key): (
                            self.alpha_modifications[key][0].strip(),
                            float(self.alpha_modifications[key][1]),
                        )
                        for key in self.alpha_modifications.keys()
                    }
                )
                if self.alpha_modifications is not None
                else None,
                "xl_position_peptide": self.alpha_peptide_crosslink_position,
                "proteins": copy.deepcopy(self.alpha_proteins),
                "xl_position_proteins": copy.deepcopy(
                    self.alpha_proteins_crosslink_positions
                ),
                "pep_position_proteins": copy.deepcopy(
                    self.alpha_proteins_peptide_positions
                ),
                "score": self.alpha_score,
                "decoy": self.alpha_decoy,
            },
            key_b: {
                "peptide": self.beta_peptide.strip(),
                "modifications": copy.deepcopy(
                    {
                        int(key): (
                            self.beta_modifications[key][0].strip(),
                            float(self.beta_modifications[key][1]),
                        )
                        for key in self.beta_modifications.keys()
                    }
                )
                if self.beta_modifications is not None
                else None,
                "xl_position_peptide": self.beta_peptide_crosslink_position,
                "proteins": copy.deepcopy(self.beta_proteins),
                "xl_position_proteins": copy.deepcopy(
                    self.beta_proteins_crosslink_positions
                ),
                "pep_position_proteins": copy.deepcopy(
                    self.beta_proteins_peptide_positions
                ),
                "score": self.beta_score,
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
        self.__dict__["alpha_modifications"] = crosslink[keys[0]]["modifications"]
        self.__dict__["alpha_peptide_crosslink_position"] = crosslink[keys[0]]["xl_position_peptide"]  # fmt: skip
        self.__dict__["alpha_proteins"] = alpha_proteins_clean
        self.__dict__["alpha_proteins_crosslink_positions"] = crosslink[keys[0]]["xl_position_proteins"]  # fmt: skip
        self.__dict__["alpha_proteins_peptide_positions"] = crosslink[keys[0]]["pep_position_proteins"]  # fmt: skip
        self.__dict__["alpha_score"] = crosslink[keys[0]]["score"]
        self.__dict__["alpha_decoy"] = crosslink[keys[0]]["decoy"]
        self.__dict__["beta_peptide"] = crosslink[keys[1]]["peptide"]
        self.__dict__["beta_modifications"] = crosslink[keys[1]]["modifications"]
        self.__dict__["beta_peptide_crosslink_position"] = crosslink[keys[1]]["xl_position_peptide"]  # fmt: skip
        self.__dict__["beta_proteins"] = beta_proteins_clean
        self.__dict__["beta_proteins_crosslink_positions"] = crosslink[keys[1]]["xl_position_proteins"]  # fmt: skip
        self.__dict__["beta_proteins_peptide_positions"] = crosslink[keys[1]]["pep_position_proteins"]  # fmt: skip
        self.__dict__["beta_score"] = crosslink[keys[1]]["score"]
        self.__dict__["beta_decoy"] = crosslink[keys[1]]["decoy"]
        if self.alpha_score is not None:
            if np.isnan(self.alpha_score):
                self.__dict__["alpha_score"] = None
        if self.beta_score is not None:
            if np.isnan(self.beta_score):
                self.__dict__["beta_score"] = None
        if self.score is not None:
            if np.isnan(self.score):
                self.__dict__["score"] = None
        if self.retention_time is not None:
            if np.isnan(self.retention_time):
                self.__dict__["retention_time"] = None
        if self.ion_mobility is not None:
            if np.isnan(self.ion_mobility):
                self.__dict__["ion_mobility"] = None
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

    def copy_with_update(self, update: Dict[str, Any] = {}) -> CrosslinkSpectrumMatch:
        r"""Creates a deep copy of the crosslink-spectrum-match with optional attribute updates.

        Parameters
        ----------
        update : dict of str, any, default = empty dict
            Dictionary mapping attribute names (str) to their updated values.
            The default (empty dict) will create a deep copy with the original
            attribute values.

        Returns
        -------
        CrosslinkSpectrumMatch
            New crosslink-spectrum-match with optionally updated attributes.

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
        >>> csm_copy = csm.copy_with_update(update={"scan_nr": 2})
        """
        _ok = check_input(update, "update", dict)
        return CrosslinkSpectrumMatch(
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
            spectrum_file=self.spectrum_file
            if "spectrum_file" not in update
            else update["spectrum_file"],
            scan_nr=self.scan_nr if "scan_nr" not in update else update["scan_nr"],
            alpha_modifications=copy.deepcopy(self.alpha_modifications)
            if "alpha_modifications" not in update
            else update["alpha_modifications"],
            alpha_proteins=copy.deepcopy(self.alpha_proteins)
            if "alpha_proteins" not in update
            else update["alpha_proteins"],
            alpha_proteins_crosslink_positions=copy.deepcopy(
                self.alpha_proteins_crosslink_positions
            )
            if "alpha_proteins_crosslink_positions" not in update
            else update["alpha_proteins_crosslink_positions"],
            alpha_proteins_peptide_positions=copy.deepcopy(
                self.alpha_proteins_peptide_positions
            )
            if "alpha_proteins_peptide_positions" not in update
            else update["alpha_proteins_peptide_positions"],
            alpha_score=self.alpha_score
            if "alpha_score" not in update
            else update["alpha_score"],
            alpha_decoy=self.alpha_decoy
            if "alpha_decoy" not in update
            else update["alpha_decoy"],
            beta_modifications=copy.deepcopy(self.beta_modifications)
            if "beta_modifications" not in update
            else update["beta_modifications"],
            beta_proteins=copy.deepcopy(self.beta_proteins)
            if "beta_proteins" not in update
            else update["beta_proteins"],
            beta_proteins_crosslink_positions=copy.deepcopy(
                self.beta_proteins_crosslink_positions
            )
            if "beta_proteins_crosslink_positions" not in update
            else update["beta_proteins_crosslink_positions"],
            beta_proteins_peptide_positions=copy.deepcopy(
                self.beta_proteins_peptide_positions
            )
            if "beta_proteins_peptide_positions" not in update
            else update["beta_proteins_peptide_positions"],
            beta_score=self.beta_score
            if "beta_score" not in update
            else update["beta_score"],
            beta_decoy=self.beta_decoy
            if "beta_decoy" not in update
            else update["beta_decoy"],
            score=self.score if "score" not in update else update["score"],
            charge=self.charge if "charge" not in update else update["charge"],
            retention_time=self.retention_time
            if "retention_time" not in update
            else update["retention_time"],
            ion_mobility=self.ion_mobility
            if "ion_mobility" not in update
            else update["ion_mobility"],
            additional_information=copy.deepcopy(self.additional_information)
            if "additional_information" not in update
            else update["additional_information"],
        )

    def to_crosslink(self) -> Crosslink:
        r"""Creates a crosslink from the crosslink-spectrum-match.

        Returns
        -------
        Crosslink
            The corresponding crosslink created from the crosslink-spectrum-match.
        """
        return create_crosslink(
            peptide_a=self.alpha_peptide,
            xl_position_peptide_a=self.alpha_peptide_crosslink_position,
            proteins_a=copy.deepcopy(self.alpha_proteins),
            xl_position_proteins_a=copy.deepcopy(
                self.alpha_proteins_crosslink_positions
            ),
            decoy_a=self.alpha_decoy,
            peptide_b=self.beta_peptide,
            xl_position_peptide_b=self.beta_peptide_crosslink_position,
            proteins_b=copy.deepcopy(self.beta_proteins),
            xl_position_proteins_b=copy.deepcopy(
                self.beta_proteins_crosslink_positions
            ),
            decoy_b=self.beta_decoy,
            score=self.score,
            additional_information=copy.deepcopy(self.additional_information),
        )

    def display(
        self,
        show_additional_information: bool = False,
        return_str: bool = False,
    ) -> None | str:
        r"""Pretty prints the crosslink-spectrum-match.

        Parameters
        ----------
        show_additional_information : bool, default = False
            Also display data in the ``additional_information``.
        return_str : bool, default = False
            If the display string should be returned.

        Returns
        -------
        None, or str
            The display string of the crosslink-spectrum-match if ``return_str = True``
            otherwise None.

        Examples
        --------
        >>> from pyXLMS import parser
        >>> pr = parser.read(
        ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        ...     engine="MS Annika",
        ...     crosslinker="DSS",
        ... )
        >>> csms = pr["crosslink-spectrum-matches"]
        >>> csms[0].display()
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
        """
        _ok = check_input(
            show_additional_information, "show_additional_information", bool
        )
        _ok = check_input(return_str, "return_str", bool)
        display: str = ""
        display += f"Data Type:                          {self.data_type}\n"
        display += f"Completeness:                       {self.completeness}\n"
        display += f"Alpha Peptide:                      {self.alpha_peptide}\n"
        display += f"Alpha Modifications:                {self.alpha_modifications}\n"
        display += f"Alpha Peptide Crosslink Position:   {self.alpha_peptide_crosslink_position}\n"
        display += f"Alpha Proteins:                     {self.alpha_proteins}\n"
        display += f"Alpha Proteins Crosslink Positions: {self.alpha_proteins_crosslink_positions}\n"
        display += f"Alpha Proteins Peptide Positions:   {self.alpha_proteins_peptide_positions}\n"
        display += f"Alpha Peptide Score:                {self.alpha_score}\n"
        display += f"Alpha Decoy:                        {self.alpha_decoy}\n"
        display += f"Beta Peptide:                       {self.beta_peptide}\n"
        display += f"Beta Modifications:                 {self.beta_modifications}\n"
        display += f"Beta Peptide Crosslink Position:    {self.beta_peptide_crosslink_position}\n"
        display += f"Beta Proteins:                      {self.beta_proteins}\n"
        display += f"Beta Proteins Crosslink Positions:  {self.beta_proteins_crosslink_positions}\n"
        display += f"Beta Proteins Peptide Positions:    {self.beta_proteins_peptide_positions}\n"
        display += f"Beta Peptide Score:                 {self.beta_score}\n"
        display += f"Beta Decoy:                         {self.beta_decoy}\n"
        display += f"Crosslink Type:                     {self.crosslink_type}\n"
        display += f"CSM Score:                          {self.score}\n"
        display += f"Spectrum File:                      {self.spectrum_file}\n"
        display += f"Scan Number:                        {self.scan_nr}\n"
        display += f"Precursor Charge:                   {self.charge}\n"
        display += f"Retention Time:                     {self.retention_time}\n"
        display += f"Ion Mobility/FAIMS CV:              {self.ion_mobility}\n"
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
        r"""Returns the Proforma string for the crosslink-spectrum-match.

        Parameters
        ----------
        crosslinker : str, or float, or None, default = None
            Optional name or mass of the crosslink reagent. If the name is given, it should be a valid
            name from XLMOD.

        Returns
        -------
        str
            The Proforma string of the crosslink-spectrum-match.

        Notes
        -----
        - Modifications with unknown mass are skipped.
        - If no modifications are given, only the crosslink modification will be encoded in the Proforma.
        - If no modifications are given and no crosslinker is given, the unmodified peptide Proforma will be returned.

        Examples
        --------
        >>> from pyXLMS.data import create_csm_min
        >>> csm = create_csm_min("PEPKTIDE", 4, "KPEPTIDE", 1, "RUN_1", 1)
        >>> csm.to_proforma()
        'KPEPTIDE//PEPKTIDE'

        >>> from pyXLMS.data import create_csm_min
        >>> csm = create_csm_min("PEPKTIDE", 4, "KPEPTIDE", 1, "RUN_1", 1)
        >>> csm.to_proforma(crosslinker="Xlink:DSSO")
        'K[Xlink:DSSO]PEPTIDE//PEPK[Xlink:DSSO]TIDE'

        >>> from pyXLMS.data import create_csm_min
        >>> csm = create_csm_min(
        ...     "PEPKTIDE",
        ...     4,
        ...     "KPMEPTIDE",
        ...     1,
        ...     "RUN_1",
        ...     1,
        ...     modifications_b={3: ("Oxidation", 15.994915)},
        ... )
        >>> csm.to_proforma(crosslinker="Xlink:DSSO")
        'K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE'

        >>> from pyXLMS.data import create_csm_min
        >>> csm = create_csm_min(
        ...     "PEPKTIDE",
        ...     4,
        ...     "KPMEPTIDE",
        ...     1,
        ...     "RUN_1",
        ...     1,
        ...     modifications_b={3: ("Oxidation", 15.994915)},
        ...     charge=3,
        ... )
        >>> csm.to_proforma(crosslinker="Xlink:DSSO")
        'K[Xlink:DSSO]PM[+15.994915]EPTIDE//PEPK[Xlink:DSSO]TIDE/3'

        >>> from pyXLMS.data import create_csm_min
        >>> csm = create_csm_min(
        ...     "PEPKTIDE",
        ...     4,
        ...     "KPMEPTIDE",
        ...     1,
        ...     "RUN_1",
        ...     1,
        ...     modifications_a={4: ("DSSO", 158.00376)},
        ...     modifications_b={1: ("DSSO", 158.00376), 3: ("Oxidation", 15.994915)},
        ...     charge=3,
        ... )
        >>> csm.to_proforma()
        'K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3'

        >>> from pyXLMS.data import create_csm_min
        >>> csm = create_csm_min(
        ...     "PEPKTIDE",
        ...     4,
        ...     "KPMEPTIDE",
        ...     1,
        ...     "RUN_1",
        ...     1,
        ...     modifications_a={4: ("DSSO", 158.00376)},
        ...     modifications_b={1: ("DSSO", 158.00376), 3: ("Oxidation", 15.994915)},
        ...     charge=3,
        ... )
        >>> csm.to_proforma(crosslinker="Xlink:DSSO")
        'K[+158.00376]PM[+15.994915]EPTIDE//PEPK[+158.00376]TIDE/3'
        """
        peptide_a = get_modified_peptide(
            self.alpha_peptide,
            self.alpha_modifications,
            self.alpha_peptide_crosslink_position,
            crosslinker,
        )
        peptide_b = get_modified_peptide(
            self.beta_peptide,
            self.beta_modifications,
            self.beta_peptide_crosslink_position,
            crosslinker,
        )
        if self.charge is not None:
            return f"{peptide_a}//{peptide_b}/{self.charge}"
        return f"{peptide_a}//{peptide_b}"


def create_csm(
    peptide_a: str,
    modifications_a: Optional[Dict[int, Tuple[str, float]]],
    xl_position_peptide_a: int,
    proteins_a: Optional[List[str]],
    xl_position_proteins_a: Optional[List[int]],
    pep_position_proteins_a: Optional[List[int]],
    score_a: Optional[float],
    decoy_a: Optional[bool],
    peptide_b: str,
    modifications_b: Optional[Dict[int, Tuple[str, float]]],
    xl_position_peptide_b: int,
    proteins_b: Optional[List[str]],
    xl_position_proteins_b: Optional[List[int]],
    pep_position_proteins_b: Optional[List[int]],
    score_b: Optional[float],
    decoy_b: Optional[bool],
    score: Optional[float],
    spectrum_file: str,
    scan_nr: int,
    charge: Optional[int],
    rt: Optional[float],
    im_cv: Optional[float],
    additional_information: Optional[Dict[str, Any]] = None,
) -> CrosslinkSpectrumMatch:
    r"""Creates a crosslink-spectrum-match data structure.

    Contains minimal data necessary for representing a single crosslink-spectrum-match. The returned crosslink-spectrum-match data structure
    is a dictionary with keys as detailed in the return section.

    Parameters
    ----------
    peptide_a : str
        The unmodified amino acid sequence of the first peptide.
    modifications_a : dict of [int, tuple], or None
        The modifications of the first peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.
        ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications should be denoted with position ``len(peptide) + 1``.
        If the peptide is not modified an empty dictionary should be given.
    xl_position_peptide_a : int
        The position of the crosslinker in the sequence of the first peptide (1-based).
    proteins_a : list of str, or None
        The accessions of proteins that the first peptide is associated with.
    xl_position_proteins_a : list of int, or None
        Positions of the crosslink in the proteins of the first peptide (1-based).
    pep_position_proteins_a : list of int, or None
        Positions of the first peptide in the corresponding proteins (1-based).
    score_a : float, or None
        Identification score of the first peptide.
    decoy_a : bool, or None
        Whether the alpha peptide is from the decoy database or not.
    peptide_b : str
        The unmodified amino acid sequence of the second peptide.
    modifications_b : dict of [int, tuple], or None
        The modifications of the second peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass.
        ``N-terminal`` modifications should be denoted with position ``0``. ``C-terminal`` modifications should be denoted with position ``len(peptide) + 1``.
        If the peptide is not modified an empty dictionary should be given.
    xl_position_peptide_b : int
        The position of the crosslinker in the sequence of the second peptide (1-based).
    proteins_b : list of str, or None
        The accessions of proteins that the second peptide is associated with.
    xl_position_proteins_b : list of int, or None
        Positions of the crosslink in the proteins of the second peptide (1-based).
    pep_position_proteins_b : list of int, or None
        Positions of the second peptide in the corresponding proteins (1-based).
    score_b : float, or None
        Identification score of the second peptide.
    decoy_b : bool, or None
        Whether the beta peptide is from the decoy database or not.
    score: float, or None
        Score of the crosslink-spectrum-match.
    spectrum_file : str
        Name of the spectrum file the crosslink-spectrum-match was identified in.
    scan_nr : int
        The corresponding scan number of the crosslink-spectrum-match.
    charge : int, or None
        The precursor charge of the corresponding mass spectrum of the crosslink-spectrum-match.
    rt : float, or None
        The retention time of the corresponding mass spectrum of the crosslink-spectrum-match in seconds.
    im_cv : float, or None
        The ion mobility or compensation voltage of the corresponding mass spectrum of the crosslink-spectrum-match.
    additional_information: dict with str keys, or None, default = None
        A dictionary with additional information associated with the crosslink-spectrum-match.

    Returns
    -------
    dict
        The dictionary representing the crosslink-spectrum-match with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_modifications``,
        ``alpha_peptide_crosslink_position``, ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_proteins_peptide_positions``,
        ``alpha_score``, ``alpha_decoy``, ``beta_peptide``, ``beta_modifications``, ``beta_peptide_crosslink_position``, ``beta_proteins``,
        ``beta_proteins_crosslink_positions``, ``beta_proteins_peptide_positions``, ``beta_score``, ``beta_decoy``, ``crosslink_type``, ``score``,
        ``spectrum_file``, ``scan_nr``, ``retention_time``, ``ion_mobility``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.
    ValueError
        If the length of crosslink positions or peptide positions is not equal to the length of proteins.

    Notes
    -----
    The minimum required data for creating a crosslink-spectrum-match is:

    - ``peptide_a``: The unmodified amino acid sequence of the first peptide.
    - ``peptide_b``: The unmodified amino acid sequence of the second peptide.
    - ``xl_position_peptide_a``: The position of the crosslinker in the sequence of the first peptide (1-based).
    - ``xl_position_peptide_b``: The position of the crosslinker in the sequence of the second peptide (1-based).
    - ``spectrum_file``: Name of the spectrum file the crosslink-spectrum-match was identified in.
    - ``scan_nr``: The corresponding scan number of the crosslink-spectrum-match.

    Examples
    --------
    >>> from pyXLMS.data import create_csm
    >>> minimal_csm = create_csm(
    ...     peptide_a="PEPTIDEA",
    ...     modifications_a={},
    ...     xl_position_peptide_a=1,
    ...     proteins_a=None,
    ...     xl_position_proteins_a=None,
    ...     pep_position_proteins_a=None,
    ...     score_a=None,
    ...     decoy_a=None,
    ...     peptide_b="PEPTIDEB",
    ...     modifications_b={},
    ...     xl_position_peptide_b=5,
    ...     proteins_b=None,
    ...     xl_position_proteins_b=None,
    ...     pep_position_proteins_b=None,
    ...     score_b=None,
    ...     decoy_b=None,
    ...     score=None,
    ...     spectrum_file="MS_EXP1",
    ...     scan_nr=1,
    ...     charge=None,
    ...     rt=None,
    ...     im_cv=None,
    ... )

    >>> from pyXLMS.data import create_csm
    >>> csm = create_csm(
    ...     peptide_a="PEPTIDEA",
    ...     modifications_a={1: ("Oxidation", 15.994915)},
    ...     xl_position_peptide_a=1,
    ...     proteins_a=["PROTEINA"],
    ...     xl_position_proteins_a=[1],
    ...     pep_position_proteins_a=[1],
    ...     score_a=20.1,
    ...     decoy_a=False,
    ...     peptide_b="PEPTIDEB",
    ...     modifications_b={},
    ...     xl_position_peptide_b=5,
    ...     proteins_b=["PROTEINB"],
    ...     xl_position_proteins_b=[3],
    ...     pep_position_proteins_b=[1],
    ...     score_b=33.7,
    ...     decoy_b=False,
    ...     score=20.1,
    ...     spectrum_file="MS_EXP1",
    ...     scan_nr=1,
    ...     charge=3,
    ...     rt=13.5,
    ...     im_cv=-50,
    ... )
    """
    return CrosslinkSpectrumMatch(
        alpha_peptide=peptide_a,
        alpha_peptide_crosslink_position=xl_position_peptide_a,
        beta_peptide=peptide_b,
        beta_peptide_crosslink_position=xl_position_peptide_b,
        spectrum_file=spectrum_file,
        scan_nr=scan_nr,
        alpha_modifications=modifications_a,
        alpha_proteins=proteins_a,
        alpha_proteins_crosslink_positions=xl_position_proteins_a,
        alpha_proteins_peptide_positions=pep_position_proteins_a,
        alpha_score=score_a,
        alpha_decoy=decoy_a,
        beta_modifications=modifications_b,
        beta_proteins=proteins_b,
        beta_proteins_crosslink_positions=xl_position_proteins_b,
        beta_proteins_peptide_positions=pep_position_proteins_b,
        beta_score=score_b,
        beta_decoy=decoy_b,
        score=score,
        charge=charge,
        retention_time=rt,
        ion_mobility=im_cv,
        additional_information=additional_information,
    )


def create_csm_min(
    peptide_a: str,
    xl_position_peptide_a: int,
    peptide_b: str,
    xl_position_peptide_b: int,
    spectrum_file: str,
    scan_nr: int,
    **kwargs,
) -> CrosslinkSpectrumMatch:
    r"""Creates a crosslink-spectrum-match data structure from minimal input.

    Contains minimal data necessary for representing a single crosslink-spectrum-match. This
    is an alias for ``data.create_csm()``that sets all optional parameters to ``None`` for convenience.
    The returned crosslink-spectrum-match data structure is a dictionary with keys as detailed in the
    return section.

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
    spectrum_file : str
        Name of the spectrum file the crosslink-spectrum-match was identified in.
    scan_nr : int
        The corresponding scan number of the crosslink-spectrum-match.
    **kwargs
        Any additional parameters will be passed to ``data.create_csm()``.

    Returns
    -------
    dict
        The dictionary representing the crosslink-spectrum-match with keys ``data_type``, ``completeness``, ``alpha_peptide``, ``alpha_modifications``,
        ``alpha_peptide_crosslink_position``, ``alpha_proteins``, ``alpha_proteins_crosslink_positions``, ``alpha_proteins_peptide_positions``,
        ``alpha_score``, ``alpha_decoy``, ``beta_peptide``, ``beta_modifications``, ``beta_peptide_crosslink_position``, ``beta_proteins``,
        ``beta_proteins_crosslink_positions``, ``beta_proteins_peptide_positions``, ``beta_score``, ``beta_decoy``, ``crosslink_type``, ``score``,
        ``spectrum_file``, ``scan_nr``, ``retention_time``, ``ion_mobility``, and ``additional_information``.
        Alpha and beta are assigned based on peptide sequence, the peptide that alphabetically comes first is assigned to alpha.

    Notes
    -----
    See also ``data.create_csm()``.

    Examples
    --------
    >>> from pyXLMS.data import create_csm_min
    >>> minimal_csm = create_csm("PEPTIDEA", 1, "PEPTIDEB", 5, "MS_EXP1", 1)
    """
    return create_csm(
        peptide_a=peptide_a,
        modifications_a=kwargs["modifications_a"]
        if "modifications_a" in kwargs
        else None,
        xl_position_peptide_a=xl_position_peptide_a,
        proteins_a=kwargs["proteins_a"] if "proteins_a" in kwargs else None,
        xl_position_proteins_a=kwargs["xl_position_proteins_a"]
        if "xl_position_proteins_a" in kwargs
        else None,
        pep_position_proteins_a=kwargs["pep_position_proteins_a"]
        if "pep_position_proteins_a" in kwargs
        else None,
        score_a=kwargs["score_a"] if "score_a" in kwargs else None,
        decoy_a=kwargs["decoy_a"] if "decoy_a" in kwargs else None,
        peptide_b=peptide_b,
        modifications_b=kwargs["modifications_b"]
        if "modifications_b" in kwargs
        else None,
        xl_position_peptide_b=xl_position_peptide_b,
        proteins_b=kwargs["proteins_b"] if "proteins_b" in kwargs else None,
        xl_position_proteins_b=kwargs["xl_position_proteins_b"]
        if "xl_position_proteins_b" in kwargs
        else None,
        pep_position_proteins_b=kwargs["pep_position_proteins_b"]
        if "pep_position_proteins_b" in kwargs
        else None,
        score_b=kwargs["score_b"] if "score_b" in kwargs else None,
        decoy_b=kwargs["decoy_b"] if "decoy_b" in kwargs else None,
        score=kwargs["score"] if "score" in kwargs else None,
        spectrum_file=spectrum_file,
        scan_nr=scan_nr,
        charge=kwargs["charge"] if "charge" in kwargs else None,
        rt=kwargs["rt"] if "rt" in kwargs else None,
        im_cv=kwargs["im_cv"] if "im_cv" in kwargs else None,
        additional_information=kwargs["additional_information"]
        if "additional_information" in kwargs
        else None,
    )


def create_crosslink_from_csm(csm: CrosslinkSpectrumMatch) -> Crosslink:
    r"""Creates a crosslink data structure from a crosslink-spectrum-match.

    Creates a crosslink data structure from a crosslink-spectrum-match. The returned crosslink data structure is a dictionary with keys
    as detailed in the return section.

    Parameters
    ----------
    csm : dict of str
        The crosslink-spectrum-match item to be converted to a crosslink item.

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
        If parameter ``csm`` is not a valid crosslink-spectrum-match.

    Notes
    -----
    See also ``data.create_crosslink()``.

    Examples
    --------
    >>> from pyXLMS.data import create_csm_min, create_crosslink_from_csm
    >>> csm = create_csm_min("PEPTIDEA", 1, "PEPTIDEB", 5, "RUN_1", 1)
    >>> crosslink = create_crosslink_from_csm(csm)
    """
    _ok = check_input(csm, "csm", CrosslinkSpectrumMatch)
    return create_crosslink(
        peptide_a=csm.alpha_peptide,
        xl_position_peptide_a=csm.alpha_peptide_crosslink_position,
        proteins_a=csm.alpha_proteins,
        xl_position_proteins_a=csm.alpha_proteins_crosslink_positions,
        decoy_a=csm.alpha_decoy,
        peptide_b=csm.beta_peptide,
        xl_position_peptide_b=csm.beta_peptide_crosslink_position,
        proteins_b=csm.beta_proteins,
        xl_position_proteins_b=csm.beta_proteins_crosslink_positions,
        decoy_b=csm.beta_decoy,
        score=csm.score,
        additional_information=csm.additional_information,
    )
