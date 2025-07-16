#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from ..data import check_input
from ..data import check_input_multi
from ..transform.util import assert_data_type_same
from .to_pyxlinkviewer import to_pyxlinkviewer

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import List


def to_xlmstools(
    crosslinks: List[Dict[str, Any]],
    pdb_file: str | BinaryIO,
    gap_open: int | float = -10.0,
    gap_extension: int | float = -1.0,
    min_sequence_identity: float = 0.8,
    allow_site_mismatch: bool = False,
    ignore_chains: List[str] = [],
    filename_prefix: Optional[str] = None,
) -> Dict[str, Any]:
    r"""Exports a list of crosslinks to xlms-tools format.

    Exports a list of crosslinks to xlms-tools format for protein structure analysis. The python package
    xlms-tools is available from
    `gitlab.com/topf-lab/xlms-tools <https://gitlab.com/topf-lab/xlms-tools>`_.
    This exporter performs basical local sequence alignment to align crosslinked peptides to a protein
    structure in PDB format. Gap open and gap extension penalties can be chosen as well as a threshold
    for sequence identity that must be satisfied in order for a match to be reported. Additionally the
    alignment is checked if the supposedly crosslinked residue can be modified with a crosslinker in
    the protein structure. Due to the alignment shift amino acids might change and a crosslink is
    reported at a position that is not able to react with the crosslinker. Optionally, these positions
    can still be reported.

    Parameters
    ----------
    crosslinks : list of dict of str, any
        A list of crosslinks.
    pdb_file : str, or file stream
        The name/path of the PDB file or a file-like object/stream. If a string is
        provided but no file is found locally, it's assumed to be an identifier and
        the file is fetched from the PDB.
    gap_open : int, or float, default = -10.0
        Gap open penalty for sequence alignment.
    gap_extension : int, or float, default = -1.0,
        Gap extension penalty for sequence alignment.
    min_sequence_identity : float, default = 0.8
        Minimum sequence identity to consider an aligned crosslinked peptide a match with
        its corresponding position in the protein structure. Should be given as a fraction
        between 0 and 1, e.g. the default of 0.8 corresponds to a minimum of 80% sequence
        identity.
    allow_site_mismatch : bool, default = False
        If the crosslink position after alignment is not a reactive amino acid in the protein
        structure, should the position still be reported. By default such cases are not reported.
    ignore_chains : list of str, default = empty list
        A list of chains to ignore in the protein structure.
    filename_prefix : str, or None, default = None
        If not None, the exported data will be written to files with the specified filename prefix.
        The full list of written files can be accessed via the returned dictionary.

    Returns
    -------
    dict of str, any
        Returns a dictionary with key ``xlms-tools`` containing the formatted text for xlms-tools,
        with key ``xlms-tools DataFrame`` containing the information from ``xlms-tools`` but as a
        pandas DataFrame, with key ``Number of mapped crosslinks`` containing the total number of mapped
        crosslinks, with key ``Mapping`` containing a string that logs how crosslinks were mapped to the
        protein structure, with key ``Parsed PDB sequence`` containing the protein sequence that was
        parsed from the PDB file, with key ``Parsed PDB chains`` containing the parsed chains from the
        PDB file, with key ``Parsed PDB residue numbers`` containing the parsed residue numbers from the
        PDB file, and with key ``Exported files`` containing a list of filenames of all files that were
        written to disk.

    Raises
    ------
    TypeError
        If a wrong data type is provided.
    TypeError
        If data contains elements of mixed data type.
    ValueError
        If parameter min_sequence_identity is out of bounds.
    ValueError
        If the provided data contains no elements.

    Notes
    -----
    Internally this exporter just calls ``exporter.to_pyxlinkviewer()`` and re-writes some of the files
    since the two tools share the same input file structure.

    Examples
    --------
    >>> from pyXLMS.exporter import to_xlmstools
    >>> from pyXLMS.parser import read_custom
    >>> pr = read_custom("data/_test/exporter/xlms-tools/unique_links_all_pyxlms.csv")
    >>> crosslinks = pr["crosslinks"]
    >>> xlmstools_result = to_xlmstools(crosslinks, pdb_file="6YHU", filename_prefix="6YHU")
    >>> xlmstools_output_file_str = xlmstools_result["xlms-tools"]
    >>> xlmstools_dataframe = xlmstools_result["xlms-tools DataFrame"]
    >>> nr_mapped_crosslinks = xlmstools_result["Number of mapped crosslinks"]
    >>> crosslink_mapping = xlmstools_result["Mapping"]
    >>> parsed_pdb_sequenece = xlmstools_result["Parsed PDB sequence"]
    >>> parsed_pdb_chains = xlmstools_result["Parsed PDB chains"]
    >>> parsed_pdb_residue_numbers = xlmstools_result["Parsed PDB residue numbers"]
    >>> exported_files = xlmstools_result["Exported files"]
    """
    _ok = check_input(crosslinks, "crosslinks", list, dict)
    _ok = check_input_multi(gap_open, "gap_open", [int, float])
    _ok = check_input_multi(gap_extension, "gap_extension", [int, float])
    _ok = check_input(min_sequence_identity, "min_sequence_identity", float)
    _ok = check_input(allow_site_mismatch, "allow_site_mismatch", bool)
    _ok = check_input(ignore_chains, "ignore_chains", list, str)
    _ok = (
        check_input(filename_prefix, "filename_prefix", str)
        if filename_prefix is not None
        else True
    )
    if min_sequence_identity < 0.0 or min_sequence_identity > 1.0:
        raise ValueError(
            "Minimum sequence identity should be given as a fraction, e.g. 0.8 for 80% minimum sequence identity!"
        )
    if len(crosslinks) == 0:
        raise ValueError("Provided data does not contain any crosslinks!")
    if "data_type" not in crosslinks[0] or crosslinks[0]["data_type"] != "crosslink":
        raise TypeError(
            "Unsupported data type for input data! Parameter data has to be a list of crosslinks!"
        )
    if not assert_data_type_same(crosslinks):
        raise TypeError("Not all elements in data have the same data type!")
    pyxlinkviewer = to_pyxlinkviewer(
        crosslinks=crosslinks,
        pdb_file=pdb_file,
        gap_open=float(gap_open),
        gap_extension=float(gap_extension),
        min_sequence_identity=min_sequence_identity,
        allow_site_mismatch=allow_site_mismatch,
        ignore_chains=ignore_chains,
        filename_prefix=None,
    )
    exported_files = list()
    parsed_pdb_str = ""
    pdb_sequence = pyxlinkviewer["Parsed PDB sequence"]
    pdb_chains = pyxlinkviewer["Parsed PDB chains"]
    pdb_residue_numbers = pyxlinkviewer["Parsed PDB residue numbers"]
    for i, r in enumerate(pdb_residue_numbers):
        parsed_pdb_str = (
            parsed_pdb_str + pdb_sequence[i] + " " + pdb_chains[i] + " " + str(r) + "\n"
        )
    fasta = f">db|PARSEDPDB|sequence parsed from PDB file\n{pdb_sequence}"
    if filename_prefix is not None:
        with open(filename_prefix + "_xlms-tools.txt", "w", encoding="utf-8") as f:
            f.write(pyxlinkviewer["PyXlinkViewer"])
            f.close()
        exported_files.append(filename_prefix + "_xlms-tools.txt")
        with open(filename_prefix + "_mapping.txt", "w", encoding="utf-8") as f:
            f.write(pyxlinkviewer["Mapping"])
            f.close()
        exported_files.append(filename_prefix + "_mapping.txt")
        with open(filename_prefix + "_parsedPDB.txt", "w", encoding="utf-8") as f:
            f.write(parsed_pdb_str)
            f.close()
        exported_files.append(filename_prefix + "_parsedPDB.txt")
        with open(filename_prefix + "_sequence.fasta", "w", encoding="utf-8") as f:
            f.write(fasta)
            f.close()
        exported_files.append(filename_prefix + "_sequence.fasta")
    return {
        "xlms-tools": pyxlinkviewer["PyXlinkViewer"],
        "xlms-tools DataFrame": pyxlinkviewer["PyXlinkViewer DataFrame"],
        "Number of mapped crosslinks": pyxlinkviewer["Number of mapped crosslinks"],
        "Mapping": pyxlinkviewer["Mapping"],
        "Parsed PDB sequence": pdb_sequence,
        "Parsed PDB chains": pdb_chains,
        "Parsed PDB residue numbers": pdb_residue_numbers,
        "Exported files": exported_files,
    }
