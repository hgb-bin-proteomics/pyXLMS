#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import re
import os
import warnings
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from biopandas.pdb import PandasPdb

from ..constants import AMINO_ACIDS_3TO1
from ..data import check_input
from ..data import check_input_multi
from ..transform.util import assert_data_type_same

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import List

# BLOSUM62
try:
    BLOSUM62 = substitution_matrices.load("BLOSUM62")
except Exception as _e:
    warnings.warn(
        RuntimeWarning("Unable to load BLOSUM62 from biopython. Using local version...")
    )
    # fmt: off
    bl62_matrix = [[ 4., -1., -2., -2.,  0., -1., -1.,  0., -2., -1., -1., -1., -1., -2., -1.,  1.,  0., -3., -2.,  0., -2., -1.,  0., -4.],
                   [-1.,  5.,  0., -2., -3.,  1.,  0., -2.,  0., -3., -2.,  2., -1., -3., -2., -1., -1., -3., -2., -3., -1.,  0., -1., -4.],
                   [-2.,  0.,  6.,  1., -3.,  0.,  0.,  0.,  1., -3., -3.,  0., -2., -3., -2.,  1.,  0., -4., -2., -3.,  3.,  0., -1., -4.],
                   [-2., -2.,  1.,  6., -3.,  0.,  2., -1., -1., -3., -4., -1., -3., -3., -1.,  0., -1., -4., -3., -3.,  4.,  1., -1., -4.],
                   [ 0., -3., -3., -3.,  9., -3., -4., -3., -3., -1., -1., -3., -1., -2., -3., -1., -1., -2., -2., -1., -3., -3., -2., -4.],
                   [-1.,  1.,  0.,  0., -3.,  5.,  2., -2.,  0., -3., -2.,  1.,  0., -3., -1.,  0., -1., -2., -1., -2.,  0.,  3., -1., -4.],
                   [-1.,  0.,  0.,  2., -4.,  2.,  5., -2.,  0., -3., -3.,  1., -2., -3., -1.,  0., -1., -3., -2., -2.,  1.,  4., -1., -4.],
                   [ 0., -2.,  0., -1., -3., -2., -2.,  6., -2., -4., -4., -2., -3., -3., -2.,  0., -2., -2., -3., -3., -1., -2., -1., -4.],
                   [-2.,  0.,  1., -1., -3.,  0.,  0., -2.,  8., -3., -3., -1., -2., -1., -2., -1., -2., -2.,  2., -3.,  0.,  0., -1., -4.],
                   [-1., -3., -3., -3., -1., -3., -3., -4., -3.,  4.,  2., -3.,  1.,  0., -3., -2., -1., -3., -1.,  3., -3., -3., -1., -4.],
                   [-1., -2., -3., -4., -1., -2., -3., -4., -3.,  2.,  4., -2.,  2.,  0., -3., -2., -1., -2., -1.,  1., -4., -3., -1., -4.],
                   [-1.,  2.,  0., -1., -3.,  1.,  1., -2., -1., -3., -2.,  5., -1., -3., -1.,  0., -1., -3., -2., -2.,  0.,  1., -1., -4.],
                   [-1., -1., -2., -3., -1.,  0., -2., -3., -2.,  1.,  2., -1.,  5.,  0., -2., -1., -1., -1., -1.,  1., -3., -1., -1., -4.],
                   [-2., -3., -3., -3., -2., -3., -3., -3., -1.,  0.,  0., -3.,  0.,  6., -4., -2., -2.,  1.,  3., -1., -3., -3., -1., -4.],
                   [-1., -2., -2., -1., -3., -1., -1., -2., -2., -3., -3., -1., -2., -4.,  7., -1., -1., -4., -3., -2., -2., -1., -2., -4.],
                   [ 1., -1.,  1.,  0., -1.,  0.,  0.,  0., -1., -2., -2.,  0., -1., -2., -1.,  4.,  1., -3., -2., -2.,  0.,  0.,  0., -4.],
                   [ 0., -1.,  0., -1., -1., -1., -1., -2., -2., -1., -1., -1., -1., -2., -1.,  1.,  5., -2., -2.,  0., -1., -1.,  0., -4.],
                   [-3., -3., -4., -4., -2., -2., -3., -2., -2., -3., -2., -3., -1.,  1., -4., -3., -2., 11.,  2., -3., -4., -3., -2., -4.],
                   [-2., -2., -2., -3., -2., -1., -2., -3.,  2., -1., -1., -2., -1.,  3., -3., -2., -2.,  2.,  7., -1., -3., -2., -1., -4.],
                   [ 0., -3., -3., -3., -1., -2., -2., -3., -3.,  3.,  1., -2.,  1., -1., -2., -2.,  0., -3., -1.,  4., -3., -2., -1., -4.],
                   [-2., -1.,  3.,  4., -3.,  0.,  1., -1.,  0., -3., -4.,  0., -3., -3., -2.,  0., -1., -4., -3., -3.,  4.,  1., -1., -4.],
                   [-1.,  0.,  0.,  1., -3.,  3.,  4., -2.,  0., -3., -3.,  1., -1., -3., -1.,  0., -1., -3., -2., -2.,  1.,  4., -1., -4.],
                   [ 0., -1., -1., -1., -2., -1., -1., -1., -1., -1., -1., -1., -1., -1., -2.,  0.,  0., -2., -1., -1., -1., -1., -1., -4.],
                   [-4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4.,  1.]
                  ]
    # fmt: on
    bl62_alphabet = "ARNDCQEGHILKMFPSTWYVBZX*"
    BLOSUM62 = substitution_matrices.Array(
        alphabet=bl62_alphabet, data=np.array(bl62_matrix)
    )


## ------------------------------ INFO ------------------------------ ##
##                                                                    ##
## This code is based on pyXlinkViewerExporter_msannika.py            ##
## retrieved: 18. June 2025                                           ##
## base version: 1.1.0                                                ##
## accessible via:                                                    ##
## https://github.com/hgb-bin-proteomics/MSAnnika_exporters/          ##
##     blob/master/pyXlinkViewerExporter_msannika.py                  ##
##                                                                    ##
## Because this code has been written a very long time ago and is     ##
## very complicated it is unfortunately largely undocumented. For me  ##
## documentation of the private functions is not a priority, so this  ##
## will probably largely stay undocumented.                           ##
##                                                                    ##
## signed -MJB                                                        ##
##                                                                    ##
## ------------------------------------------------------------------ ##


def __get_pdb(pdb_file: str | BinaryIO) -> PandasPdb:
    if isinstance(pdb_file, str):
        if os.path.isfile(pdb_file):
            return PandasPdb().read_pdb(pdb_file)
        return PandasPdb().fetch_pdb(pdb_file.split(".pdb")[0])
    lines = pdb_file.readlines()
    return PandasPdb().read_pdb_from_list(lines)


def __get_pdb_data(
    pdb_file: str | BinaryIO, ignore_chains: List[str]
) -> Dict[str, Any]:
    pdb_df = __get_pdb(pdb_file)

    sequence = list()
    chains = list()
    residue_numbers = dict()
    residue_numbers_lst = list()

    for i, row in pdb_df.df["ATOM"].iterrows():
        three_letter_aa = row["residue_name"]
        residue_number = row["residue_number"]
        chain = row["chain_id"]

        if three_letter_aa.strip() in AMINO_ACIDS_3TO1:
            residue = AMINO_ACIDS_3TO1[three_letter_aa.strip()]
            if chain not in ignore_chains:
                if chain in residue_numbers:
                    if residue_number not in residue_numbers[chain]:
                        residue_numbers[chain].append(residue_number)
                        sequence.append(residue)
                        chains.append(chain)
                else:
                    residue_numbers[chain] = [residue_number]
                    sequence.append(residue)
                    chains.append(chain)
        else:
            warnings.warn(
                RuntimeWarning(
                    f"Found amino acid: {three_letter_aa.strip()}, which is not a supported amino acid!"
                )
            )

    for chain_id in sorted(residue_numbers.keys()):
        residue_numbers_lst = residue_numbers_lst + residue_numbers[chain_id]

    return {
        "sequence": "".join(sequence),
        "chains": "".join(chains),
        "residue_numbers": residue_numbers_lst,
    }


def __calculate_sequence_identity(sequence_a: str, sequence_b: str) -> float:
    ident = 0
    for i in range(len(sequence_a)):
        if sequence_a[i] == sequence_b[i]:
            ident += 1
    return float(ident / len(sequence_a))


def __calculate_shifted_xl_pos(alignment: str, xl_pos_in_pep: int) -> int:
    new_xl_pos = xl_pos_in_pep
    if len(alignment) <= xl_pos_in_pep:
        return len(alignment) - 1
    if "-" not in alignment[: xl_pos_in_pep + 1]:
        return xl_pos_in_pep
    else:
        gaps = [m.start() for m in re.finditer("-", alignment)]
        curr_limit = xl_pos_in_pep + 1
        for gap in gaps:
            if gap < curr_limit:
                curr_limit += 1
                new_xl_pos += 1
        return new_xl_pos


def __get_pep_pos(candidates: List[int], alignment_position: int) -> int:
    distances = dict()
    for candidate in candidates:
        distances[abs(candidate - alignment_position)] = candidate
    return distances[sorted(distances.keys())[0]]


def __get_xl_position_and_chain_in_protein(
    pdb_sequence: str,
    pdb_chains: str,
    pdb_residue_numbers: str,
    peptide_sequence: str,
    crosslink_position_in_peptide: int,
    gap_open: float,
    gap_extension: float,
    min_sequence_identity: float,
    allow_site_mismatch: bool,
) -> List[str]:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extension
    aligner.substitution_matrix = BLOSUM62
    pep_seq = peptide_sequence
    pep_pos_in_proteins = [m.start() for m in re.finditer(pep_seq, pdb_sequence)]
    xl_pos_in_pep = crosslink_position_in_peptide - 1
    if len(pep_pos_in_proteins) == 0:
        alignments = sorted(
            aligner.align(Seq(pdb_sequence), Seq(pep_seq)),
            key=lambda alignment: alignment.score,  # pyright: ignore[reportAttributeAccessIssue]
            reverse=True,
        )
        if len(alignments) == 0:
            return []
        else:
            top_alignment = alignments[0]
            if top_alignment.coordinates is not None:
                a_start = top_alignment.coordinates[0][0]
                a_end = top_alignment.coordinates[0][1]
                b_start = top_alignment.coordinates[1][0]
                b_end = top_alignment.coordinates[1][1]
            else:
                raise RuntimeError("Could not extract positions of alignment!")
            seqA = top_alignment.target[a_start:a_end]
            seqB = top_alignment.query[b_start:b_end]
            sequence_identity = __calculate_sequence_identity(seqA, seqB)
            if sequence_identity > min_sequence_identity:
                xl_pos_in_alignment = xl_pos_in_pep
                if len(pep_seq) != len(seqB):
                    xl_pos_in_alignment = __calculate_shifted_xl_pos(
                        seqB, xl_pos_in_pep
                    )
                    ## not sure what this block does anymore, also I am pretty sure
                    ## that seqA should not have gaps, because we align seqB to seqA?
                    ## but maybe it should be xl_pos_in_alignment -= len()?
                    # if "-" in seqA[:xl_pos_in_alignment]:
                    #    xl_pos_in_alignment - len(
                    #        [m for m in re.finditer("-", seqA[:xl_pos_in_alignment])]
                    #    )
                    ## xl_pos_in_alignment -= len() should make sense since if there
                    ## is gaps in seqA we need to subtract those to get the correct
                    ## position?
                    ## reminder: be cautious
                    if "-" in seqA[:xl_pos_in_alignment]:
                        xl_pos_in_alignment -= len(
                            [m for m in re.finditer("-", seqA[:xl_pos_in_alignment])]
                        )
                if not allow_site_mismatch:
                    if xl_pos_in_alignment < xl_pos_in_pep:
                        return []
                    if seqA[xl_pos_in_alignment] == seqB[xl_pos_in_alignment]:
                        pep_pos_in_protein = __get_pep_pos(
                            [
                                m.start()
                                for m in re.finditer(
                                    seqA.replace("-", ""), pdb_sequence
                                )
                            ],
                            a_start,
                        )
                        xl_position = pep_pos_in_protein + xl_pos_in_alignment
                        xl_chain = pdb_chains[xl_position]
                        xl_residue = pdb_residue_numbers[xl_position]
                        return [str(xl_residue) + "|" + str(xl_chain) + "|"]
                    else:
                        return []
                else:
                    pep_pos_in_protein = __get_pep_pos(
                        [
                            m.start()
                            for m in re.finditer(seqA.replace("-", ""), pdb_sequence)
                        ],
                        a_start,
                    )
                    xl_position = pep_pos_in_protein + xl_pos_in_alignment
                    xl_chain = pdb_chains[xl_position]
                    xl_residue = pdb_residue_numbers[xl_position]
                    return [str(xl_residue) + "|" + str(xl_chain) + "|"]
            else:
                return []
    else:
        chain_residues = []
        for pep_pos_in_protein in pep_pos_in_proteins:
            xl_position = pep_pos_in_protein + xl_pos_in_pep
            xl_chain = pdb_chains[xl_position]
            xl_residue = pdb_residue_numbers[xl_position]
            chain_residues.append(str(xl_residue) + "|" + str(xl_chain) + "|")
        return chain_residues


def __to_dataframe(pyxlinkviewer: str) -> pd.DataFrame:
    pos_a = list()
    chain_a = list()
    pos_b = list()
    chain_b = list()
    for line in pyxlinkviewer.split("\n"):
        if len(line.strip()) > 0:
            xl = line.split("|")
            pos_a.append(int(xl[0]))
            chain_a.append(xl[1].strip())
            pos_b.append(int(xl[2]))
            chain_b.append(xl[3].strip())
    return pd.DataFrame(
        {"residue 1": pos_a, "chain 1": chain_a, "residue 2": pos_b, "chain 2": chain_b}
    )


def __to_pyxlinkviewer(
    crosslinks: List[Dict[str, Any]],
    pdb_file: str | BinaryIO,
    gap_open: float,
    gap_extension: float,
    min_sequence_identity: float,
    allow_site_mismatch: bool,
    ignore_chains: List[str],
    filename_prefix: Optional[str],
) -> Dict[str, Any]:
    pdb_data = __get_pdb_data(pdb_file, ignore_chains)
    pdb_sequence = pdb_data["sequence"]
    pdb_chains = pdb_data["chains"]
    pdb_residue_numbers = pdb_data["residue_numbers"]

    if len(pdb_sequence) == len(pdb_chains) == len(pdb_residue_numbers):
        # all good, do nothing
        pass
    else:
        raise RuntimeError(
            "Parsed PDB sequence, chain and residue numbers are not matching!"
        )

    output_string = ""
    mapping_string = ""
    nr_of_mapped_xl = 0
    for i, crosslink in enumerate(crosslinks):
        links_a = __get_xl_position_and_chain_in_protein(
            pdb_sequence=pdb_sequence,
            pdb_chains=pdb_chains,
            pdb_residue_numbers=pdb_residue_numbers,
            peptide_sequence=crosslink["alpha_peptide"],
            crosslink_position_in_peptide=crosslink["alpha_peptide_crosslink_position"],
            gap_open=gap_open,
            gap_extension=gap_extension,
            min_sequence_identity=min_sequence_identity,
            allow_site_mismatch=allow_site_mismatch,
        )
        links_b = __get_xl_position_and_chain_in_protein(
            pdb_sequence=pdb_sequence,
            pdb_chains=pdb_chains,
            pdb_residue_numbers=pdb_residue_numbers,
            peptide_sequence=crosslink["beta_peptide"],
            crosslink_position_in_peptide=crosslink["beta_peptide_crosslink_position"],
            gap_open=gap_open,
            gap_extension=gap_extension,
            min_sequence_identity=min_sequence_identity,
            allow_site_mismatch=allow_site_mismatch,
        )
        if len(links_a) != 0 and len(links_b) != 0:
            for link_a in links_a:
                for link_b in links_b:
                    output_string = output_string + link_a + link_b + "\n"
                    mapping_string = (
                        mapping_string
                        + link_a
                        + link_b
                        + "\n"
                        + crosslink["alpha_peptide"]
                        + " - "
                        + crosslink["beta_peptide"]
                        + "\n"
                    )
                    nr_of_mapped_xl += 1
    exported_files = list()
    parsed_pdb_str = ""
    for i, r in enumerate(pdb_residue_numbers):
        parsed_pdb_str = (
            parsed_pdb_str + pdb_sequence[i] + " " + pdb_chains[i] + " " + str(r) + "\n"
        )
    fasta = f">db|PARSEDPDB|sequence parsed from PDB file\n{pdb_sequence}"
    if filename_prefix is not None:
        with open(filename_prefix + "_PyXlinkViewer.txt", "w", encoding="utf-8") as f:
            f.write(output_string)
            f.close()
        exported_files.append(filename_prefix + "_PyXlinkViewer.txt")
        with open(filename_prefix + "_mapping.txt", "w", encoding="utf-8") as f:
            f.write(mapping_string)
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
        "PyXlinkViewer": output_string,
        "PyXlinkViewer DataFrame": __to_dataframe(output_string),
        "Number of mapped crosslinks": nr_of_mapped_xl,
        "Mapping": mapping_string,
        "Parsed PDB sequence": pdb_sequence,
        "Parsed PDB chains": pdb_chains,
        "Parsed PDB residue numbers": pdb_residue_numbers,
        "Exported files": exported_files,
    }


def to_pyxlinkviewer(
    crosslinks: List[Dict[str, Any]],
    pdb_file: str | BinaryIO,
    gap_open: int | float = -10.0,
    gap_extension: int | float = -1.0,
    min_sequence_identity: float = 0.8,
    allow_site_mismatch: bool = False,
    ignore_chains: List[str] = [],
    filename_prefix: Optional[str] = None,
) -> Dict[str, Any]:
    r"""Exports a list of crosslinks to PyXlinkViewer format.

    Exports a list of crosslinks to PyXlinkViewer format for visualization in pyMOL. The tool
    PyXlinkViewer is available from
    `github.com/BobSchiffrin/PyXlinkViewer <https://github.com/BobSchiffrin/PyXlinkViewer>`_.
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
        Returns a dictionary with key ``PyXlinkViewer`` containing the formatted text for PyXlinkViewer,
        with key ``PyXlinkViewer DataFrame`` containing the information from ``PyXlinkViewer`` but as a
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

    Examples
    --------
    >>> from pyXLMS.exporter import to_pyxlinkviewer
    >>> from pyXLMS.parser import read_custom
    >>> pr = read_custom("data/_test/exporter/pyxlinkviewer/unique_links_all_pyxlms.csv")
    >>> crosslinks = pr["crosslinks"]
    >>> pyxlinkviewer_result = to_pyxlinkviewer(crosslinks, pdb_file="6YHU", filename_prefix="6YHU")
    >>> pyxlinkviewer_output_file_str = pyxlinkviewer_result["PyXlinkViewer"]
    >>> pyxlinkviewer_dataframe = pyxlinkviewer_result["PyXlinkViewer DataFrame"]
    >>> nr_mapped_crosslinks = pyxlinkviewer_result["Number of mapped crosslinks"]
    >>> crosslink_mapping = pyxlinkviewer_result["Mapping"]
    >>> parsed_pdb_sequenece = pyxlinkviewer_result["Parsed PDB sequence"]
    >>> parsed_pdb_chains = pyxlinkviewer_result["Parsed PDB chains"]
    >>> parsed_pdb_residue_numbers = pyxlinkviewer_result["Parsed PDB residue numbers"]
    >>> exported_files = pyxlinkviewer_result["Exported files"]
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
    return __to_pyxlinkviewer(
        crosslinks,
        pdb_file,
        float(gap_open),
        float(gap_extension),
        min_sequence_identity,
        allow_site_mismatch,
        ignore_chains,
        filename_prefix,
    )
