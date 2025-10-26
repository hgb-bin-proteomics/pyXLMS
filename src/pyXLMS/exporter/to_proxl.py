#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import os
from lxml import etree  # pyright: ignore[reportAttributeAccessIssue]
import urllib.request as ur
from Bio.SeqIO.FastaIO import SimpleFastaParser

from ..data import check_input
from .to_proxl_util import __local_schema
from .util import __get_filename
from .to_alphalink2 import __protein_supported_by_crosslink
from ..transform.util import get_available_keys
from ..transform.util import modifications_to_str as mts
from ..constants import MODIFICATIONS

from typing import Optional
from typing import Callable
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

SCHEMA_URL = "https://github.com/yeastrc/proxl-import-api/raw/refs/heads/master/xsd/proxl-xml.xsd"


def __build_header(
    fasta_filename: str,
    search_engine: str,
    version: str,
    score: Literal["higher_better", "lower_better"],
    crosslinker: str,
    crosslinker_mass: float,
) -> List[str]:
    r"""Builds the header section of the ProXL XML.

    Parameters
    ----------
    fasta_filename : str
        The name/path of the fasta file for reading protein sequences.
    search_engine : str
        Name of the used crosslink search engine.
    version : str
        Version identifier of the used crosslink search engine.
    score : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float
        Monoisotopic delta mass of the crosslink modification.

    Returns
    -------
    list of str
        A list of lines of the header section of the ProXL XML.

    Notes
    -----
    This function should not be called directly, it is called from ``to_proxl()``.
    """
    filter_direction = "above" if score == "higher_better" else "below"
    lines = [
        r"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>""",
        f'<proxl_input fasta_filename="{fasta_filename}">',
        r"""<search_program_info>""",
        r"""<search_programs>""",
        f'<search_program name="{search_engine}" display_name="{search_engine}" version="{version}">',
        r"""<psm_annotation_types>""",
        r"""<filterable_psm_annotation_types>""",
        f'<filterable_psm_annotation_type name="score" description="Score of the search engine" filter_direction="{filter_direction}" default_filter="false"/>',
        r"""</filterable_psm_annotation_types>""",
        r"""<descriptive_psm_annotation_types>""",
        r"""<descriptive_psm_annotation_type name="spectrum filename" description="Name of the MS file"/>""",
        r"""<descriptive_psm_annotation_type name="scan number" description="Scan number"/>""",
        r"""</descriptive_psm_annotation_types>""",
        r"""</psm_annotation_types>""",
        r"""</search_program>""",
        r"""</search_programs>""",
        r"""<default_visible_annotations>""",
        r"""<visible_psm_annotations>""",
        f'<search_annotation search_program="{search_engine}" annotation_name="score"/>',
        f'<search_annotation search_program="{search_engine}" annotation_name="spectrum filename"/>',
        f'<search_annotation search_program="{search_engine}" annotation_name="scan number"/>',
        r"""</visible_psm_annotations>""",
        r"""</default_visible_annotations>""",
        r"""</search_program_info>""",
        r"""<linkers>""",
        f'<linker name="{crosslinker}">',
        r"""<crosslink_masses>""",
        f'<crosslink_mass mass="{crosslinker_mass}"/>',
        r"""</crosslink_masses>""",
        r"""</linker>""",
        r"""</linkers>""",
    ]
    return lines


def __get_reported_peptide_string(csm: Dict[str, Any]) -> str:
    r"""Creates a unique 'reported_peptide_string' for a crosslink-spectrum-match.

    Parameters
    ----------
    csm : dict of str, any
        A crosslink-spectrum-match.

    Returns
    -------
    str
        The unique 'reported_peptide_string' of the crosslink-spectrum-match.

    Notes
    -----
    This function should not be called directly, it is called from ``to_proxl()``.
    """
    return (
        f"{csm['alpha_peptide']}({csm['alpha_peptide_crosslink_position']})-{csm['beta_peptide']}({csm['beta_peptide_crosslink_position']})"
        f"_({mts(csm['alpha_modifications'])})-({mts(csm['beta_modifications'])})"
    )


def __get_reported_peptides(
    csms: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    r"""Groups crosslink-spectrum-matches by their unique 'reported_peptide_string'.

    Parameters
    ----------
    csms : list of dict of str, any
        A list of crosslink-spectrum-matches.

    Returns
    -------
    dict of str, list of dict of str, any
        A dictionary that maps unique 'reported_peptide_string' keys to lists of associated
        crosslink-spectrum-matches as values.

    Notes
    -----
    This function should not be called directly, it is called from ``to_proxl()``.
    """
    reported_peptides = dict()
    for csm in csms:
        reported_peptide_string = __get_reported_peptide_string(csm)
        if reported_peptide_string in reported_peptides:
            reported_peptides[reported_peptide_string].append(csm)
        else:
            reported_peptides[reported_peptide_string] = [csm]
    return reported_peptides


def __build_psm(
    csm: Dict[str, Any], crosslinker_mass: float, search_engine: str
) -> List[str]:
    r"""Builds the 'psm' section of the ProXL XML.

    Parameters
    ----------
    csms : dict of str, any
        A crosslink-spectrum-match.
    crosslinker_mass : float
        Monoisotopic delta mass of the crosslink modification.
    search_engine : str
        Name of the used crosslink search engine.

    Returns
    -------
    list of str
        A list of lines of one 'psm' section of the ProXL XML.

    Notes
    -----
    This function should not be called directly, it is called from ``to_proxl()``.
    """
    scan_file_name = os.path.splitext(csm["spectrum_file"])[0] + ".mzML"
    lines = [
        f'<psm scan_file_name="{scan_file_name}" scan_number="{csm["scan_nr"]}" precursor_charge="{csm["charge"]}" linker_mass="{crosslinker_mass}">',
        r"""<filterable_psm_annotations>""",
        f'<filterable_psm_annotation search_program="{search_engine}" annotation_name="score" value="{csm["score"]}"/>',
        r"""</filterable_psm_annotations>""",
        r"""<descriptive_psm_annotations>""",
        f'<descriptive_psm_annotation search_program="{search_engine}" annotation_name="spectrum filename" value="{scan_file_name}"/>',
        f'<descriptive_psm_annotation search_program="{search_engine}" annotation_name="scan number" value="{csm["scan_nr"]}"/>',
        r"""</descriptive_psm_annotations>""",
        r"""</psm>""",
    ]
    return lines


def __build_modifications(csm: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    r"""Builds the 'modifications' section of the ProXL XML.

    Parameters
    ----------
    csm : dict of str, any
        A crosslink-spectrum-match.

    Returns
    -------
    tuple of list of str, list of str
        A list of lines of one 'modifications' section of the ProXL XML for each crosslinked
        peptide.

    Notes
    -----
    This function should not be called directly, it is called from ``to_proxl()``.
    """
    modifications_a = list()
    if csm["alpha_modifications"] is not None:
        modifications_a.append(r"""<modifications>""")
        for pos, modification in csm["alpha_modifications"].items():
            if pos == 0:
                pos = 1
            if pos > len(csm["alpha_peptide"]):
                pos = len(csm["alpha_peptide"])
            if pos != csm["alpha_peptide_crosslink_position"]:
                modifications_a.append(
                    f'<modification mass="{modification[1]}" position="{pos}" isMonolink="false"/>'
                )
        modifications_a.append(r"""</modifications>""")
    if len(modifications_a) <= 2:
        modifications_a.clear()
    modifications_b = list()
    if csm["beta_modifications"] is not None:
        modifications_b.append(r"""<modifications>""")
        for pos, modification in csm["beta_modifications"].items():
            if pos == 0:
                pos = 1
            if pos > len(csm["beta_peptide"]):
                pos = len(csm["beta_peptide"])
            if pos != csm["beta_peptide_crosslink_position"]:
                modifications_b.append(
                    f'<modification mass="{modification[1]}" position="{pos}" isMonolink="false"/>'
                )
        modifications_b.append(r"""</modifications>""")
    if len(modifications_b) <= 2:
        modifications_b.clear()
    return (modifications_a, modifications_b)


def __build_reported_peptides(
    reported_peptides: Dict[str, List[Dict[str, Any]]],
    crosslinker_mass: float,
    search_engine: str,
) -> List[str]:
    r"""Builds the 'reported_peptides' section of the ProXL XML.

    Parameters
    ----------
    reported_peptides : dict of str, list of dict of str, any
        A dictionary that maps unique 'reported_peptide_string' keys to lists of associated
        crosslink-spectrum-matches as values.
    crosslinker_mass : float
        Monoisotopic delta mass of the crosslink modification.
    search_engine : str
        Name of the used crosslink search engine.

    Returns
    -------
    list of str
        A list of lines of the 'reported_peptides' section of the ProXL XML.

    Notes
    -----
    This function should not be called directly, it is called from ``to_proxl()``.
    """
    lines = [r"""<reported_peptides>"""]
    for reported_peptide in reported_peptides:
        example_csm = reported_peptides[reported_peptide][0]
        modifications = __build_modifications(example_csm)
        lines += [
            f'<reported_peptide reported_peptide_string="{reported_peptide}" type="crosslink">',
            r"""<peptides>""",
            f'<peptide sequence="{example_csm["alpha_peptide"]}">',
        ]
        lines += modifications[0]
        lines += [
            r"""<linked_positions>""",
            f'<linked_position position="{example_csm["alpha_peptide_crosslink_position"]}"/>',
            r"""</linked_positions>""",
            r"""</peptide>""",
            f'<peptide sequence="{example_csm["beta_peptide"]}">',
        ]
        lines += modifications[1]
        lines += [
            r"""<linked_positions>""",
            f'<linked_position position="{example_csm["beta_peptide_crosslink_position"]}"/>',
            r"""</linked_positions>""",
        ]
        lines.append(r"""</peptide>""")
        lines.append(r"""</peptides>""")
        lines.append(r"""<psms>""")
        for csm in reported_peptides[reported_peptide]:
            lines += __build_psm(csm, crosslinker_mass, search_engine)
        lines.append(r"""</psms>""")
        lines.append(r"""</reported_peptide>""")
    lines.append(r"""</reported_peptides>""")
    return lines


def __build_matched_proteins(
    csms: List[Dict[str, Any]],
    fasta_filename: str,
    title_to_accession: Optional[Callable[[str], str]],
) -> List[str]:
    r"""Builds the 'matched_proteins' section of the ProXL XML.

    Parameters
    ----------
    csms : list of dict of str, any
        A list of crosslink-spectrum-matches.
    fasta_filename : str
        The name/path of the fasta file for reading protein sequences.
    title_to_accession : callable, or None
        A function that parses the protein accession from the fasta title/header. If None (default)
        the full fasta headers are used. An example function would be ``transform.fasta_title_to_accession()``.

    Returns
    -------
    list of str
        A list of lines of the 'matched_proteins' section of the ProXL XML.

    Notes
    -----
    This function should not be called directly, it is called from ``to_proxl()``.
    """
    lines = [r"""<matched_proteins>"""]
    fasta_items = list()
    with open(fasta_filename, "r", encoding="utf-8") as f:
        for item in SimpleFastaParser(f):
            fasta_items.append(item)
    for item in fasta_items:
        header = item[0].replace('"', "")
        name = title_to_accession(header) if title_to_accession is not None else header
        sequence = item[1]
        if __protein_supported_by_crosslink(sequence, csms):
            lines.append(f'<protein sequence="{sequence}">')
            lines.append(f'<protein_annotation name="{name}"/>')
            lines.append(r"""</protein>""")
    lines.append(r"""</matched_proteins>""")
    return lines


def __validate_schema(
    xml_str: str, schema_validation: Literal["online", "offline"]
) -> bool:
    r"""Validates a ProXL XML string against the ProXL XML schema.

    Parameters
    ----------
    xml_str : str
        The ProXL XML string to validate.
    schema_validation : str, one of "online" or "offline"
        If XML schema validation should use the most recent online schema or the locally stored but
        possibly outdated offline schema.

    Returns
    -------
    bool
        True if the ProXL XML string validates successfully against the schema, otherwise False.

    Notes
    -----
    This function should not be called directly, it is called from ``to_proxl()``.
    """
    proxl_schema = __local_schema.encode("utf-8")
    if schema_validation == "online":
        proxl_schema = ur.urlopen(SCHEMA_URL).read()
    parser = etree.XMLParser(encoding="utf-8")
    schema_doc = etree.fromstring(proxl_schema, parser=parser)
    schema = etree.XMLSchema(schema_doc)
    xml_doc = etree.fromstring(xml_str.encode("utf-8"), parser=parser)
    return schema.validate(xml_doc)


def to_proxl(
    csms: List[Dict[str, Any]],
    fasta_filename: str,
    search_engine: str,
    search_engine_version: str,
    score: Literal["higher_better", "lower_better"],
    crosslinker: str,
    crosslinker_mass: Optional[float] = None,
    modifications: Dict[str, float] = MODIFICATIONS,
    fasta_filename_override: Optional[str] = None,
    fasta_title_to_accession: Optional[Callable[[str], str]] = None,
    filename: Optional[str] = None,
    schema_validation: Literal["online", "offline"] = "online",
) -> str:
    r"""Exports a list of crosslink-spectrum-matches to ProXL format.

    Exports a list of crosslink-spectrum-matches to ProXL format. The tool ProXL is accessible
    via the link `yeastrc.org/proxl_public <https://www.yeastrc.org/proxl_public/>`_.
    Requires that ``charge`` and ``score`` fields are set for all crosslink-spectrum-matches.

    Parameters
    ----------
    csms : list of dict of str, any
        A list of crosslink-spectrum-matches.
    fasta_filename : str
        The name/path of the fasta file for reading protein sequences.
    search_engine : str
        Name of the used crosslink search engine.
    search_engine_version : str
        Version identifier of the used crosslink search engine.
    score : str, one of "higher_better" or "lower_better"
        If a higher score is considered better, or a lower score is considered better.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float, or None, default = None
        Monoisotopic delta mass of the crosslink modification. If the crosslinker is
        defined in parameter "modifications" this can be omitted.
    modifications: dict of str, float, default = ``constants.MODIFICATIONS``
        Mapping of modification names to modification masses. By default uses ``constants.MODIFICATIONS``.
    fasta_filename_override : str, or None, default = None
        Name that should be used in the ProXL XML for the fasta file. If None (default)
        uses the filename of parameter 'fasta_filename' (preceding directories are pruned).
    fasta_title_to_accession : callable, or None, default = None
        A function that parses the protein accession from the fasta title/header. If None (default)
        the full fasta headers are used. An example function would be ``transform.fasta_title_to_accession()``.
    filename : str, or None, default = None
        If not None, the exported data will be written to a file with the specified filename.
    schema_validation : str, one of "online" or "offline", default = "online"
        If XML schema validation should use the most recent online schema or the locally stored but
        possibly outdated offline schema.

    Returns
    -------
    str
        The ProXL XML as a string.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If 'csms' parameter contains elements of mixed data type.
    TypeError
        If parameter score is not one of 'higher_better' or 'lower_better'.
    TypeError
        If parameter schema_validation is not one of 'online' or 'offline'.
    ValueError
        If the provided 'csms' parameter contains no elements.
    KeyError
        If the specified crosslinker could not be found/mapped.
    RuntimeError
        If not all of the required information is present in the input data.
    RuntimeError
        If the created ProXL XML file fails validation against the ProXL XML schema.

    Examples
    --------
    >>> from pyXLMS.pipelines import pipeline
    >>> from pyXLMS.exporter import to_proxl
    >>> pr = pipeline(
    ...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
    ...     engine="MS Annika",
    ...     crosslinker="DSS",
    ... )
    >>> xml = to_proxl(
    ...     pr["crosslink-spectrum-matches"],
    ...     fasta_filename="data/_fasta/Cas9_plus10.fasta",
    ...     search_engine="MS Annika",
    ...     search_engine_version="3.0.1",
    ...     score="higher_better",
    ...     crosslinker="DSS",
    ...     filename="DSS_Cas9_ProXL.xml",
    ... )
    """
    _ok = check_input(csms, "csms", list, dict)
    _ok = check_input(fasta_filename, "fasta_filename", str)
    _ok = check_input(search_engine, "search_engine", str)
    _ok = check_input(search_engine_version, "search_engine_version", str)
    _ok = check_input(score, "score", str)
    _ok = check_input(crosslinker, "crosslinker", str)
    _ok = (
        check_input(crosslinker_mass, "crosslink_mass", float)
        if crosslinker_mass is not None
        else True
    )
    _ok = check_input(modifications, "modifications", dict, float)
    _ok = (
        check_input(fasta_filename_override, "fasta_filename_override", str)
        if fasta_filename_override is not None
        else True
    )
    _ok = (
        check_input(fasta_title_to_accession, "title_to_accession", Callable)
        if fasta_title_to_accession is not None
        else True
    )
    _ok = check_input(filename, "filename", str) if filename is not None else True
    _ok = check_input(schema_validation, "schema_validation", str)
    if score not in ["higher_better", "lower_better"]:
        raise TypeError(
            "Parameter 'score' has to be one of 'higher_better' or 'lower_better'!"
        )
    if schema_validation not in ["online", "offline"]:
        raise TypeError(
            "Parameter 'schema_validation' has to be one of 'online' or 'offline'!"
        )
    if crosslinker_mass is None:
        if crosslinker not in modifications:
            raise KeyError(
                "Cannot infer crosslinker mass because crosslinker is not defined in "
                "parameter 'modifications'. Please specify crosslinker mass manually!"
            )
        else:
            crosslinker_mass = modifications[crosslinker]
    if len(csms) == 0:
        raise ValueError("Provided crosslink-spectrum-matches contain no elements!")
    if "data_type" not in csms[0] or csms[0]["data_type"] != "crosslink-spectrum-match":
        raise TypeError(
            "Unsupported data type for input csms! Parameter csms has to be a list of crosslink-spectrum-matches!"
        )
    available_keys = get_available_keys(csms)
    if not available_keys["score"] or not available_keys["charge"]:
        raise RuntimeError(
            "Can't export to ProXL because not all necessary information is available!"
        )
    fasta_name = (
        os.path.basename(fasta_filename)
        if fasta_filename_override is None
        else fasta_filename_override
    )
    reported_peptides = __get_reported_peptides(csms)
    lines = (
        __build_header(
            fasta_name,
            search_engine,
            search_engine_version,
            score,
            crosslinker,
            crosslinker_mass,
        )
        + __build_reported_peptides(reported_peptides, crosslinker_mass, search_engine)
        + __build_matched_proteins(csms, fasta_filename, fasta_title_to_accession)
        + [r"""</proxl_input>"""]
    )
    xml_str = "\n".join(lines)
    if __validate_schema(xml_str, schema_validation):
        print(
            f"Successfully created ProXL XML and validated it against {schema_validation} XML schema!"
        )
    else:
        raise RuntimeError(
            f"Created ProXL XML but validation against {schema_validation} XML schema failed!"
        )
    if filename is not None:
        with open(__get_filename(filename, "xml"), "w", encoding="utf-8") as f:
            f.write(xml_str)
            f.close()
    return xml_str
