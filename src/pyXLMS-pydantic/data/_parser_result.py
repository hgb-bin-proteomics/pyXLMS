#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import copy
import numpy as np
from pydantic import BaseModel
from pydantic import computed_field

from ._csm import CrosslinkSpectrumMatch
from ._crosslink import Crosslink
from ._util import check_input
from ._util import check_indexing

from typing import override
from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any

class ParserResult(BaseModel):
    search_engine: str,
    crosslink_spectrum_matches: Optional[List[CrosslinkSpectrumMatch]]
    crosslinks: Optional[List[Crosslink]]
    
    @computed_field
    @property
    def data_type(self) -> str:
        return "parser_result"
    
    @computed_field
    @property
    def completeness(self) -> str:
        if self.crosslink_spectrum_matches is not None and self.crosslinks is not None:
            return "full"
        if self.crosslink_spectrum_matches is None and self.crosslinks not None:
            return "empty"
        return "partial"
    
    def __getitem__(self, key: str) -> Any:
        if key == "crosslink-spectrum-matches":
            return self.crosslink_spectrum_matches
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(f"'{key}' is not a valid field!")
    
    
    def csms(self) -> List[CrosslinkSpectrumMatch] | None:
        return self.crosslink_spectrum_matches
    
    
    def xls(self) -> List[Crosslink] | None:
        return self.crosslinks
    
    
    def display(
        self,
        show_additional_information: bool = False,
        return_str: bool = False,
    ) -> None | str:
        r"""Pretty prints the parser_result.
    
        Parameters
        ----------
        show_additional_information : bool, default = False
            Also display data in the ``additional_information``.
        return_str : bool, default = False
            If the display string should be returned.
    
        Returns
        -------
        None, or str
            The display string of the crosslink-spectrum-match, crosslink, or parser_result
            if ``return_str = True`` otherwise None.
    
        Examples
        --------
        >>> from pyXLMS import parser
        >>> from pyXLMS import transform
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
        _ok = check_input(show_additional_information, "show_additional_information", bool)
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
    csms: Optional[CrosslinkSpectrumMatch]],
    crosslinks: Optional[List[Crosslink]]],
) -> ParserResult:
    r"""Creates a parser result data structure.

    Contains all necessary data elements that should be contained in a result returned by a crosslink search engine result parser.

    Parameters
    ----------
    search_engine : str
        Name of the identifying crosslink search engine.
    csms : list of dict, or None
        List of crosslink-spectrum-matches as created by ``data.create_csm()``.
    crosslinks : list of dict, or None
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
