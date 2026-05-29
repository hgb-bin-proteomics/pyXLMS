#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import copy
from pydantic import BaseModel
from pydantic import Field
from pydantic import ConfigDict
from pydantic import computed_field

from ._csm import CrosslinkSpectrumMatch
from ._crosslink import Crosslink
from ._util import check_input

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


class ParserResult(BaseModel):
    r"""Core data structure for parser results.

    Data structure returned by any (parser) function that reads crosslink-spectrum-matches
    and/or crosslinks.

    Attributes Summary
    ------------------
    Here is a short summary about the parser result attributes, for more details
    on the specific Pydantic validation requirements please refer to the corresponding attributes
    themselves.

    Required
    ^^^^^^^^
    The following attributes are required:

    search_engine : str
        The name of the identifying crosslink search engine.

    Optional
    ^^^^^^^^
    The following attributes are optional:

    crosslink_spectrum_matches : list of CrosslinkSpectrumMatch, or None, default = None
        List of parsed crosslink-spectrum-matches.
    crosslinks : list of Crosslink, or None, default = None
        List of parsed crosslinks.

    Examples
    --------
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

    search_engine: Annotated[
        str,
        Field(
            frozen=True,
            description="The name of the identifying crosslink search engine.",
        ),
    ]
    r"""
    The name of the identifying crosslink search engine.
    """
    crosslink_spectrum_matches: Annotated[
        Optional[List[CrosslinkSpectrumMatch]],
        Field(frozen=True, description="List of parsed crosslink-spectrum-matches."),
    ] = None
    r"""
    List of parsed crosslink-spectrum-matches.
    """
    crosslinks: Annotated[
        Optional[List[Crosslink]],
        Field(frozen=True, description="List of parsed crosslinks."),
    ] = None
    r"""
    List of parsed crosslinks.
    """
    model_config = ConfigDict(
        validate_assignment=True, strict=True, str_strip_whitespace=True
    )
    r"""
    Pydantic configuration for the underlying validation model.
    """

    @computed_field(description="Data type of the object.")
    @property
    def data_type(self) -> Literal["parser_result"]:
        r"""
        Data type of the object.
        """
        return "parser_result"

    @computed_field(description="Completeness of the parser result.")
    @property
    def completeness(self) -> Literal["full", "partial", "empty"]:
        r"""
        Completeness of the parser result, e.g. ``"full"`` if all attributes
        are not ``None``, ``"empty"`` if crosslink-spectrum-matches and crosslinks
        are ``None``, and otherwise ``"partial"``.
        """
        if self.crosslink_spectrum_matches is not None and self.crosslinks is not None:
            return "full"
        if self.crosslink_spectrum_matches is None and self.crosslinks is None:
            return "empty"
        return "partial"

    def __getitem__(self, key: str) -> Any:
        r"""
        Support for dict-like access.
        """
        # this is for legacy support
        if key == "crosslink-spectrum-matches":
            return self.crosslink_spectrum_matches
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(f"'{key}' is not a valid field!")

    def __contains__(self, key: str) -> bool:
        r"""
        Support for ``in`` operator.
        """
        # this is for legacy support
        if key == "crosslink-spectrum-matches":
            return True
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

    def copy_with_update(self, update: Dict[str, Any] = {}) -> ParserResult:
        r"""Creates a deep copy of the parser result with optional attribute updates.

        Parameters
        ----------
        update : dict of str, any, default = empty dict
            Dictionary mapping attribute names (str) to their updated values.
            The default (empty dict) will create a deep copy with the original
            attribute values.

        Returns
        -------
        ParserResult
            New parser result with optionally updated attributes.

        Examples
        --------
        >>> from pyXLMS.data import Crosslink
        >>> from pyXLMS.data import ParserResult
        >>> pr = ParserResult(search_engine="My Search Engine")
        >>> xl = Crosslink(
        ...     alpha_peptide="PEKP",
        ...     alpha_peptide_crosslink_position=3,
        ...     beta_peptide="TKIDE",
        ...     beta_peptide_crosslink_position=2,
        ... )
        >>> pr_copy = pr.copy_with_update(update={"crosslinks": [xl]})
        """
        _ok = check_input(update, "update", dict)
        if (
            "crosslink_spectrum_matches" in update
            and "crosslink-spectrum-matches" in update
        ):
            raise ValueError(
                "Dict 'update' must only contain key 'crosslink_spectrum_matches' "
                "or key 'crosslink-spectrum-matches' but not both!"
            )
        new_csms = copy.deepcopy(self.crosslink_spectrum_matches)
        if "crosslink_spectrum_matches" in update:
            new_csms = update["crosslink_spectrum_matches"]
        if "crosslink-spectrum-matches" in update:
            new_csms = update["crosslink-spectrum-matches"]
        return ParserResult(
            search_engine=self.search_engine
            if "search_engine" not in update
            else update["search_engine"],
            crosslink_spectrum_matches=new_csms,
            crosslinks=copy.deepcopy(self.crosslinks)
            if "crosslinks" not in update
            else update["crosslinks"],
        )

    def csms(self, create_copy: bool = True) -> List[CrosslinkSpectrumMatch] | None:
        r"""Shorthand function to retrieve crosslink-spectrum-matches.

        Parameters
        ----------
        create_copy : bool, default = True
            Whether a deep copy of the crosslink-spectrum-matches should be returned
            (default) or ``self.crosslink_spectrum_matches`` directly.

        Returns
        -------
        list of CrosslinkSpectrumMatch, or None
            Returns (a deep copy of) ``self.crosslink_spectrum_matches``.

        Notes
        -----
        Please be aware that by default this explicitly creates a deep copy of the
        underlying data!
        """
        if create_copy:
            return copy.deepcopy(self.crosslink_spectrum_matches)
        return self.crosslink_spectrum_matches

    def xls(self, create_copy: bool = True) -> List[Crosslink] | None:
        r"""Shorthand function to retrieve crosslinks.

        Parameters
        ----------
        create_copy : bool, default = True
            Whether a deep copy of the crosslinks should be returned (default) or
            ``self.crosslinks`` directly.

        Returns
        -------
        list of Crosslink, or None
            Returns (a deep copy of) ``self.crosslinks``.

        Notes
        -----
        Please be aware that by default this explicitly creates a deep copy of the
        underlying data!
        """
        if create_copy:
            return copy.deepcopy(self.crosslinks)
        return self.crosslinks

    def display(
        self,
        show_additional_information: bool = False,
        return_str: bool = False,
    ) -> None | str:
        r"""Pretty prints the parser result.

        Parameters
        ----------
        show_additional_information : bool, default = False
            Also display data in the ``additional_information``.
        return_str : bool, default = False
            If the display string should be returned.

        Returns
        -------
        None, or str
            The display string of the parser result if ``return_str = True`` otherwise None.

        Examples
        --------
        >>> from pyXLMS import parser
        >>> pr = parser.read(
        ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        ...     engine="MS Annika",
        ...     crosslinker="DSS",
        ... )
        >>> pr.display()
        Data Type:                            parser_result
        Completeness:                         full
        Identifying Search Engine:            MS Annika
        Number of Crosslink-Spectrum-Matches: 826
        Number of Crosslinks:                 300
        """
        _ok = check_input(
            show_additional_information, "show_additional_information", bool
        )
        _ok = check_input(return_str, "return_str", bool)
        display: str = ""
        csms = self.crosslink_spectrum_matches
        xls = self.crosslinks
        display += f"Data Type:                            {self.data_type}\n"
        display += f"Completeness:                         {self.completeness}\n"
        display += f"Identifying Search Engine:            {self.search_engine}\n"
        display += f"Number of Crosslink-Spectrum-Matches: {len(csms) if csms is not None else None}\n"
        display += f"Number of Crosslinks:                 {len(xls) if xls is not None else None}\n"
        display = display.strip()
        print(display)
        if return_str:
            return display
        return


def create_parser_result(
    search_engine: str,
    csms: Optional[List[CrosslinkSpectrumMatch]] = None,
    crosslinks: Optional[List[Crosslink]] = None,
) -> ParserResult:
    r"""Creates a parser result data structure.

    Contains all necessary data elements that should be contained in a result returned by a crosslink search engine result parser.

    Parameters
    ----------
    search_engine : str
        Name of the identifying crosslink search engine.
    csms : list of dict, or None, default = None
        List of crosslink-spectrum-matches as created by ``data.create_csm()``.
    crosslinks : list of dict, or None, default = None
        List of crosslinks as created by ``data.create_crosslink()``.

    Returns
    -------
    dict
        The parser result data structure which is a dictionary with keys ``data_type``, ``completeness``, ``search_engine``, ``crosslink-spectrum-matches`` and
        ``crosslinks``.

    Examples
    --------
    >>> from pyXLMS.data import create_parser_result
    >>> result = create_parser_result("MS Annika", None, None)
    >>> result["data_type"]
    'parser_result'
    >>> result["completeness"]
    'empty'
    >>> result["search_engine"]
    'MS Annika'
    """
    return ParserResult(
        search_engine=search_engine,
        crosslink_spectrum_matches=csms,
        crosslinks=crosslinks,
    )
