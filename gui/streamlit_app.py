#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "streamlit>=1.50.0",
#   "pyxlms>=1.7.0",
#   "xlsxwriter",
# ]
# ///

# pyXLMS GUI
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# ruff: noqa: F841

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

from __future__ import annotations

import io
import os
import gzip
import json
import pickle
import pandas as pd
from tempfile import NamedTemporaryFile
from tempfile import TemporaryDirectory

from pyXLMS import parser
from pyXLMS import transform
from pyXLMS import constants
from pyXLMS import plotting
from pyXLMS import exporter
from pyXLMS import __version__ as __pyxlms_version__

import streamlit as st
from streamlit.runtime.uploaded_file_manager import UploadedFile

from typing import Optional
from typing import Dict
from typing import List
from typing import Set
from typing import Any

# legacy
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


__version__ = "1.3.4"

HELP_URL = "https://pyxlms.dev/docs/webapp"


@st.cache_data
def to_text(data: str) -> bytes:
    return data.encode("utf-8")


@st.cache_data
def to_json(data: Dict[str, Any]) -> bytes:
    return json.dumps(data).encode("utf-8")


@st.cache_data
def dataframe_to_csv_stream(
    dataframe: pd.DataFrame, sep: str, index: bool, header: bool = True
) -> bytes:
    return dataframe.to_csv(sep=sep, index=index, header=header).encode("utf-8")


@st.cache_data
def dataframe_to_xlsx_stream(
    dataframe: pd.DataFrame, sheet_name: str, index: bool
) -> io.BytesIO:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:  # pyright: ignore[reportArgumentType]
        dataframe.to_excel(writer, index=index, sheet_name=sheet_name)
        writer.close()
    return buffer


@st.cache_data
def read_files(
    uploaded_files: List[UploadedFile],
    engine: str,
    crosslinker: str,
    parse_modifications: bool,
    crosslinker_mass: Optional[float],
) -> Dict[str, Any]:
    #
    with TemporaryDirectory() as d:  # pyright: ignore[reportCallIssue]
        filenames = list()
        for uploaded_file in uploaded_files:
            filename = os.path.join(d, uploaded_file.name)
            with open(filename, "wb") as f:
                f.write(uploaded_file.getbuffer())
                f.close()
            filenames.append(filename)
        if crosslinker_mass is not None:
            try:
                return parser.read(
                    filenames,
                    engine=engine,  # pyright: ignore[reportArgumentType] # ty: ignore[invalid-argument-type]
                    crosslinker=crosslinker,
                    parse_modifications=parse_modifications,
                    crosslinker_mass=crosslinker_mass,
                )
            except Exception as _e:
                parser.read(
                    filenames,
                    engine=engine,  # pyright: ignore[reportArgumentType] # ty: ignore[invalid-argument-type]
                    crosslinker=crosslinker,
                    parse_modifications=parse_modifications,
                )
        return parser.read(
            filenames,
            engine=engine,  # pyright: ignore[reportArgumentType] # ty: ignore[invalid-argument-type]
            crosslinker=crosslinker,
            parse_modifications=parse_modifications,
        )


@st.cache_data
def reannotating_positions(
    pr: List[Dict[str, Any]] | Dict[str, Any], uploaded_fasta: io.BytesIO
) -> List[Dict[str, Any]] | Dict[str, Any]:
    #
    with NamedTemporaryFile(
        suffix=os.path.splitext(uploaded_fasta.name)[1], delete_on_close=False
    ) as f:  # pyright: ignore[reportCallIssue]
        f.write(uploaded_fasta.getbuffer())
        f.close()
        return transform.reannotate_positions(pr, f.name)


@st.cache_data
def export_pyxlinkviewer_using_pdbfile(
    crosslinks: List[Dict[str, Any]], uploaded_pdb_file: io.BytesIO
) -> Dict[str, Any]:
    #
    with NamedTemporaryFile(
        suffix=os.path.splitext(uploaded_pdb_file.name)[1], delete_on_close=False
    ) as f:  # pyright: ignore[reportCallIssue]
        f.write(uploaded_pdb_file.getbuffer())
        f.close()
        return exporter.to_pyxlinkviewer(crosslinks, f.name, filename_prefix=None)


def pyxlinkviewer_get_fasta(sequence: str) -> str:
    return f">db|PARSEDPDB|sequence parsed from PDB file\n{sequence}"


def pyxlinkviewer_get_annotation(
    sequence: str, chains: str, residue_numbers: List[Any]
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Amino Acid": [c for c in sequence],
            "Chain": [c for c in chains],
            "Residue Number": residue_numbers,
        }
    )


@st.cache_data
def export_xlmstools_using_pdbfile(
    crosslinks: List[Dict[str, Any]], uploaded_pdb_file: io.BytesIO
) -> Dict[str, Any]:
    #
    with NamedTemporaryFile(
        suffix=os.path.splitext(uploaded_pdb_file.name)[1], delete_on_close=False
    ) as f:  # pyright: ignore[reportCallIssue]
        f.write(uploaded_pdb_file.getbuffer())
        f.close()
        return exporter.to_xlmstools(crosslinks, f.name, filename_prefix=None)


@st.cache_data
def filter_proteins(
    data: List[Dict[str, Any]], proteins: Set[str] | List[str]
) -> List[Dict[str, Any]]:
    filtered = transform.filter_proteins(data, proteins)
    return filtered["Both"] + filtered["One"]


@st.cache_data
def export_alphalink2(
    crosslinks: List[Dict[str, Any]], fasta: io.BytesIO, annotated_fdr: float
) -> Dict[str, Any]:
    #
    with NamedTemporaryFile(
        suffix=os.path.splitext(fasta.name)[1], delete_on_close=False
    ) as f:  # pyright: ignore[reportCallIssue]
        f.write(fasta.getbuffer())
        f.close()
        return exporter.to_alphalink2(
            crosslinks, f.name, annotated_fdr, filename_prefix=None
        )


@st.cache_data
def export_proxl(
    csms: List[Dict[str, Any]],
    fasta_file: io.BytesIO,
    search_engine: str,
    search_engine_version: str,
    score: Literal["higher_better", "lower_better"],
    crosslinker: str,
    crosslinker_mass: Optional[float],
) -> str:
    #
    with NamedTemporaryFile(
        suffix=os.path.splitext(fasta_file.name)[1], delete_on_close=False
    ) as f:  # pyright: ignore[reportCallIssue]
        f.write(fasta_file.getbuffer())
        f.close()
        return exporter.to_proxl(
            csms,
            f.name,
            search_engine,
            search_engine_version,
            score,
            crosslinker,
            crosslinker_mass,
            fasta_filename_override=fasta_file.name,
        )


@st.cache_data
def pickle_and_gzip(data: Any) -> io.BytesIO:
    buffer = io.BytesIO()
    with gzip.open(buffer, "wb") as f:
        pickle.dump(data, f)
    return buffer


def layout_plots(plots: List[Any]) -> None:
    for i in range(0, len(plots), 2):
        l_col, r_col = st.columns(2)
        if i < len(plots):
            with l_col:
                plot1 = st.pyplot(
                    plots[i],
                    transparent=True,
                    bbox_inches="tight",
                    width="stretch",
                )
        if i + 1 < len(plots):
            with r_col:
                plot1 = st.pyplot(
                    plots[i + 1],
                    transparent=True,
                    bbox_inches="tight",
                    width="stretch",
                )
    return


def reset_exports() -> None:
    # CSMs
    st.session_state["export_csms_impxfdr"] = None
    st.session_state["export_csms_msannika"] = None
    st.session_state["export_csms_proxl"] = None
    st.session_state["export_csms_xifdr"] = None
    # crosslinks
    st.session_state["export_crosslinks_alphalink2"] = None
    st.session_state["export_crosslinks_impxfdr"] = None
    st.session_state["export_crosslinks_msannika"] = None
    st.session_state["export_crosslinks_pyxlinkviewer"] = None
    st.session_state["export_crosslinks_xinet"] = None
    st.session_state["export_crosslinks_xiview"] = None
    st.session_state["export_crosslinks_xlinkdb"] = None
    st.session_state["export_crosslinks_xlmstools"] = None
    st.session_state["export_crosslinks_xmas"] = None
    # aggregated crosslinks
    st.session_state["export_aggregated_crosslinks_alphalink2"] = None
    st.session_state["export_aggregated_crosslinks_impxfdr"] = None
    st.session_state["export_aggregated_crosslinks_msannika"] = None
    st.session_state["export_aggregated_crosslinks_pyxlinkviewer"] = None
    st.session_state["export_aggregated_crosslinks_xinet"] = None
    st.session_state["export_aggregated_crosslinks_xiview"] = None
    st.session_state["export_aggregated_crosslinks_xlinkdb"] = None
    st.session_state["export_aggregated_crosslinks_xlmstools"] = None
    st.session_state["export_aggregated_crosslinks_xmas"] = None


# input tab
def input_tab():
    general_description = """
    _a python package to process protein cross-linking data_

    **pyXLMS** is a python package and web application with graphical user interface that aims to simplify and streamline the intermediate step of
    connecting crosslink search engine results with down-stream analysis tools, enabling researchers even without bioinformatics knowledge to
    conduct in-depth crosslink analyses and shifting the focus from data transformation to data interpretation and therefore gaining biological
    insight.

    Currently pyXLMS supports input from seven different crosslink search engines:
    [MaxLynx (part of MaxQuant)](https://www.maxquant.org/),
    [MeroX](https://www.stavrox.com/),
    [MS Annika](https://github.com/hgb-bin-proteomics/MSAnnika),
    [pLink 2 and pLink 3](http://pfind.ict.ac.cn/se/plink/),
    [Scout](https://github.com/diogobor/Scout),
    [xiSearch](https://www.rappsilberlab.org/software/xisearch/) and [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
    [XlinkX](https://docs.thermofisher.com/r/XlinkX-3.2-Quick-Start-Guide/),
    as well as the [mzIdentML format](https://www.psidev.info/mzidentml) of the HUPO Proteomics Standards Initiative,
    and a well-documented and [human-readable custom tabular format](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md).

    Down-stream analysis is facilitated by functionality that is directly available within pyXLMS such as validation, annotation, aggregation,
    filtering, and visualization - and [much more](https://hgb-bin-proteomics.github.io/pyXLMS/modules.html) - of crosslink-spectrum-matches and crosslinks.

    In addition, the data can easily be exported to the required data format of the various available down-stream analysis tools such as
    [AlphaLink2](https://github.com/Rappsilber-Laboratory/AlphaLink2),
    [ProXL](https://www.yeastrc.org/proxl_public/),
    [xiNET](https://crosslinkviewer.org/index.php),
    [xiVIEW](https://www.xiview.org/index.php),
    [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
    [XlinkDB](https://xlinkdb.gs.washington.edu/xlinkdb/),
    [xlms-tools](https://gitlab.com/topf-lab/xlms-tools),
    PyMOL (via [PyXlinkViewer](https://github.com/BobSchiffrin/PyXlinkViewer)),
    ChimeraX (via [XMAS](https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle)),
    or [IMP-X-FDR](https://github.com/vbc-proteomics-org/imp-x-fdr).

    **Try it yourself below!** 😉
    """
    description = st.markdown(general_description)

    header_1 = st.subheader("File Upload", divider="grey")

    uploaded_files = st.file_uploader(
        "Upload one or more cross-linking result files from any of the supported search engines or formats:",
        type=None,
        accept_multiple_files=True,
        key="uploaded_files",
        help="Upload one or more cross-linking result files from any of the supported search engines or formats.",
    )

    with st.popover(
        "Unsure which files to upload? Click me!",
        help="Display help for file selection.",
        icon="💡",
        width="stretch",
    ):
        uploaded_files_description = (
            """
            - ➡️ Check out this 🔗[**overview**](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/supported_io.md)
            on supported input and output formats.
            - ➡️ In our GitHub repository you can also find 🔗[**example files**](https://github.com/hgb-bin-proteomics/pyXLMS/tree/master/data).
            - ➡️ If you need any further help, visit the web application guide which you can find 🔗[**here**]("""
            + HELP_URL
            + """).
            - ➡️ Still stuck? Please 📩[**send us a message**](https://github.com/hgb-bin-proteomics/pyXLMS?tab=readme-ov-file#contact).
            """
        )
        uploaded_files_helper = st.markdown(uploaded_files_description)

    l1, r1 = st.columns(2)

    with l1:
        search_engine = st.selectbox(
            "Select a crosslink search engine or file format:",
            options=[
                "Custom",
                "MaxQuant",
                "MaxLynx",
                "MeroX",
                "MS Annika",
                "mzIdentML",
                "pLink",
                "Scout",
                "xiSearch/xiFDR",
                "XlinkX",
            ],
            index=None,
            key="search_engine",
            help="The crosslink search engine or file format of the uploaded file.",
        )

    with r1:
        crosslinker = st.selectbox(
            "Select the used cross-linking reagent:",
            options=["Custom"] + list(constants.CROSSLINKERS.keys()),
            index=None,
            key="crosslinker",
            help="The crosslinker used in the experiment of the uploaded result file.",
        )

    l2, center_2, r2 = st.columns(3)

    with l2:
        parse_modifications = st.toggle(
            "Parse modifications",
            value=False,
            key="parse_modifications",
            help="If post-translational-modifications should be parsed or not.",
        )

    with center_2:
        reannotate_positions = st.toggle(
            "(Re-)Annotate crosslink positions",
            value=False,
            key="reannotate_positions",
            help="If crosslink positions in proteins should be (re-)annotated.",
        )

    with r2:
        unique = st.toggle(
            "Filter for unique crosslink spectrum matches and crosslinks",
            value=False,
            key="unique",
            help="If crosslink spectrum matches and crosslinks should be filtered for unique matches only.",
        )

    l3, center_3, r3 = st.columns(3)

    with l3:
        validate = st.toggle(
            "Validate via FDR estimation",
            value=False,
            key="validate",
            help="If crosslink spectrum matches and crosslinks should be validated via false discovery rate estimation.",
        )

    with center_3:
        targets_only = st.toggle(
            "Filter for target matches only",
            value=False,
            key="targets_only",
            help="If crosslink spectrum matches and crosslinks should be filtered to contain only target-target matches.",
        )

    with r3:
        aggregate = st.toggle(
            "Aggregate crosslink spectrum matches to crosslinks",
            value=False,
            key="aggregate",
            help="If crosslink spectrum matches should be aggregated and grouped to crosslinks/residue pairs.",
        )

    l4, center_4, r4 = st.columns(3)

    with l4:
        validate_aggregated = st.toggle(
            "Validate aggregated crosslinks via FDR estimation",
            value=False,
            key="validate_aggregated",
            help="If aggregated crosslinks should be validated via false discovery rate estimation.",
        )

    if parse_modifications:
        parse_modifications_warning = st.warning(
            "pyXLMS is only able to parse a limited number of the most common post-translational modifications! "
            + "If parsing fails it is safer to leave option 'parse_modifications' turned off! More nuanced control "
            + "for additional and custom modifications is available in the pyXLMS python package!"
        )
        parse_modifications_info = st.info(
            "Lists of currently supported modifications for "
            + "[MeroX](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.html#pyXLMS.constants.MEROX_MODIFICATION_MAPPING), "
            + "[Scout](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.html#pyXLMS.constants.SCOUT_MODIFICATION_MAPPING), "
            + "[xiSearch/xiFDR](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.html#pyXLMS.constants.XI_MODIFICATION_MAPPING), "
            + "and [all other search engines](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.html#pyXLMS.constants.MODIFICATIONS)."
        )

    crosslinker_name = None
    crosslinker_mass = None
    if crosslinker == "Custom":
        l5, r5 = st.columns(2)

        with l5:
            crosslinker_name = st.text_input(
                "Name of the crosslinker:",
                value=None,
                max_chars=50,
                placeholder="DSSO",
                key="crosslinker_name",
                help="Name of the crosslinker used in the experiment of the uploaded result file.",
            )

        with r5:
            crosslinker_mass = st.number_input(
                "Mass of the crosslinker:",
                value=None,
                step=0.00001,
                format="%0.5f",
                placeholder="158.00376",
                key="crosslinker_mass",
                help="Monoisotopic delta mass of the crosslinker used in the experiment of the uploaded result file.",
            )

    if reannotate_positions:
        uploaded_fasta = st.file_uploader(
            "Upload a FASTA file for (re-)annotation:",
            type="fasta",
            accept_multiple_files=False,
            key="uploaded_fasta",
            help="Upload a FASTA file containing protein sequences for the provided crosslink spectrum matches and crosslinks.",
        )

        reannotate_positions_info = st.info(
            "In order to (re-)annotate your crosslink-spectrum-matches or crosslinks, your results should either already be pre-filtered to only contain "
            + "target-target matches, or you should also enable the **'Filter for target matches only'** option. The (re-)annotation of decoy "
            + "matches is not supported by pyXLMS."
        )

    if unique or aggregate:
        group_by = st.selectbox(
            "Group crosslinks by:",
            options=[
                "Peptide sequence and peptide crosslink position",
                "Protein Crosslink Position",
            ],
            index=0,
            key="group_by",
            help="If crosslinks should be grouped by peptide sequence and peptide crosslink position, or by protein crosslink position.",
        )

    if validate or validate_aggregated:
        l6, center_6, r6 = st.columns(3)

        with l6:
            fdr = st.number_input(
                "Target FDR:",
                value=0.01,
                min_value=0.0,
                max_value=1.0,
                step=0.001,
                format="%0.3f",
                key="fdr",
                help="The target FDR, must be given as a real number between 0 and 1. The default of 0.01 corresponds to 1% FDR.",
            )

        with center_6:
            formula = st.selectbox(
                "FDR formula:",
                options=["(TD+DD)/TT", "(TD-DD)/TT"],
                index=0,
                key="formula",
                help="Which formula to use to estimate FDR. D and DD denote decoy matches, T and TT denote target matches, and TD denotes target-decoy and decoy-target matches.",
            )

        with r6:
            separate = st.selectbox(
                "Separate intra and inter FDR?",
                options=[
                    "Concatenated FDR for intra and inter matches",
                    "Separate FDR for intra and inter matches",
                ],
                index=0,
                key="separate",
                help="If FDR should be estimated separately for intra and inter matches.",
            )

    if unique or aggregate or validate or validate_aggregated:
        score = st.selectbox(
            "Is a higher identification score considered better?",
            options=["Higher better", "Lower better"],
            index=0,
            key="score",
            help="If a higher score is considered better, or a lower score is considered better.",
        )

    l7, center_7, r7 = st.columns(3)

    with center_7:
        read_files_button = st.button("Read file(s)!", type="primary", width="stretch")

    # read in all inputs
    if read_files_button:
        # reset meta information
        if "meta_info" in st.session_state:
            del st.session_state["meta_info"]
        # reset pr and aggregated on file read
        if "pr" in st.session_state:
            del st.session_state["pr"]
        if "aggregated" in st.session_state:
            del st.session_state["aggregated"]
        # reset any exported files
        reset_exports()
        # reset proteins
        st.session_state["possible_proteins"] = None
        # check what is uploaded and set
        if uploaded_files is None or len(uploaded_files) == 0:
            _ = st.error("You need to upload at least one result file first!")
        if search_engine is None:
            _ = st.error("You need to select a search engine or format first!")
        if crosslinker is None:
            _ = st.error("You need to select a crosslinker first!")
        if crosslinker == "Custom":
            if crosslinker_name is None:
                _ = st.error("You need to specify a name for your custom crosslinker!")
            if crosslinker_mass is None:
                _ = st.error(
                    "You need to specify the crosslinker mass for your custom crosslinker!"
                )
        if (
            uploaded_files is not None
            and len(uploaded_files) > 0
            and search_engine is not None
            and crosslinker is not None
            and crosslinker != "Custom"
        ):
            with st.spinner("Parsing file...", show_time=True):
                try:
                    st.session_state["meta_info"] = {
                        "search_engine": search_engine,
                        "crosslinker_name": crosslinker,
                        "crosslinker_mass": constants.CROSSLINKERS[crosslinker],
                    }
                    st.session_state["pr"] = read_files(
                        uploaded_files,
                        search_engine,
                        crosslinker,
                        parse_modifications,
                        crosslinker_mass,
                    )
                    if aggregate:
                        if st.session_state["pr"]["crosslink-spectrum-matches"] is None:
                            st.session_state["aggregated"] = None
                        else:
                            st.session_state["aggregated"] = transform.aggregate(
                                st.session_state["pr"]["crosslink-spectrum-matches"],
                                by="peptide"
                                if group_by  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                == "Peptide sequence and peptide crosslink position"
                                else "protein",
                                score="higher_better"
                                if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                else "lower_better",
                            )
                    if validate_aggregated:
                        if st.session_state["aggregated"] is not None:
                            st.session_state["aggregated"] = transform.validate(
                                st.session_state["aggregated"],
                                fdr=fdr,  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                formula=formula,  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType] # ty: ignore[possibly-unresolved-reference, invalid-argument-type]
                                score="higher_better"
                                if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                else "lower_better",
                                separate_intra_inter=separate  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                == "Separate FDR for intra and inter matches",
                            )
                    if unique:
                        st.session_state["pr"] = transform.unique(
                            st.session_state["pr"],
                            by="peptide"
                            if group_by  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            == "Peptide sequence and peptide crosslink position"
                            else "protein",
                            score="higher_better"
                            if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            else "lower_better",
                        )
                    if validate:
                        st.session_state["pr"] = transform.validate(
                            st.session_state["pr"],
                            fdr=fdr,  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            formula=formula,  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType] # ty: ignore[possibly-unresolved-reference, invalid-argument-type]
                            score="higher_better"
                            if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            else "lower_better",
                            separate_intra_inter=separate  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            == "Separate FDR for intra and inter matches",
                        )
                    if targets_only:
                        st.session_state["pr"] = transform.targets_only(
                            st.session_state["pr"]
                        )
                        if (
                            "aggregated" in st.session_state
                            and st.session_state["aggregated"] is not None
                        ):
                            st.session_state["aggregated"] = transform.targets_only(
                                st.session_state["aggregated"]
                            )
                    if reannotate_positions:
                        if uploaded_fasta is None:  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            raise ValueError(
                                "Can't annotate crosslink position when no FASTA file is uploaded!"
                            )
                        if not targets_only:
                            _ = st.warning(
                                "Might not be able to (re-)annotate positions if results contain decoy matches!",
                                icon="⚠️",
                            )
                        st.session_state["pr"] = reannotating_positions(
                            st.session_state["pr"],
                            uploaded_fasta,  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                        )
                        if (
                            "aggregated" in st.session_state
                            and st.session_state["aggregated"] is not None
                        ):
                            st.session_state["aggregated"] = reannotating_positions(
                                st.session_state["aggregated"],
                                uploaded_fasta,  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            )
                    st.rerun()
                except Exception as e:
                    # reset meta information
                    if "meta_info" in st.session_state:
                        del st.session_state["meta_info"]
                    # reset pr and aggregated on file read
                    if "pr" in st.session_state:
                        del st.session_state["pr"]
                    if "aggregated" in st.session_state:
                        del st.session_state["aggregated"]
                    # reset any exported files
                    reset_exports()
                    # reset proteins
                    st.session_state["possible_proteins"] = None
                    _ = st.error(
                        "Something went wrong! This is most likely due to missing information in the results!",
                        icon="⚠️",
                    )
                    with st.expander("Show exception"):
                        _ = st.exception(e)
        elif (
            uploaded_files is not None
            and len(uploaded_files) > 0
            and search_engine is not None
            and crosslinker is not None
            and crosslinker == "Custom"
            and crosslinker_name is not None
            and crosslinker_mass is not None
        ):
            with st.spinner("Parsing file...", show_time=True):
                try:
                    st.session_state["meta_info"] = {
                        "search_engine": search_engine,
                        "crosslinker_name": crosslinker_name,
                        "crosslinker_mass": crosslinker_mass,
                    }
                    st.session_state["pr"] = read_files(
                        uploaded_files,
                        search_engine,
                        crosslinker_name,
                        parse_modifications,
                        crosslinker_mass,
                    )
                    if aggregate:
                        if st.session_state["pr"]["crosslink-spectrum-matches"] is None:
                            st.session_state["aggregated"] = None
                        else:
                            st.session_state["aggregated"] = transform.aggregate(
                                st.session_state["pr"]["crosslink-spectrum-matches"],
                                by="peptide"
                                if group_by  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                == "Peptide sequence and peptide crosslink position"
                                else "protein",
                                score="higher_better"
                                if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                else "lower_better",
                            )
                    if validate_aggregated:
                        if st.session_state["aggregated"] is not None:
                            st.session_state["aggregated"] = transform.validate(
                                st.session_state["aggregated"],
                                fdr=fdr,  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                formula=formula,  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType] # ty: ignore[possibly-unresolved-reference, invalid-argument-type]
                                score="higher_better"
                                if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                else "lower_better",
                                separate_intra_inter=separate  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                                == "Separate FDR for intra and inter matches",
                            )
                    if unique:
                        st.session_state["pr"] = transform.unique(
                            st.session_state["pr"],
                            by="peptide"
                            if group_by  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            == "Peptide sequence and peptide crosslink position"
                            else "protein",
                            score="higher_better"
                            if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            else "lower_better",
                        )
                    if validate:
                        st.session_state["pr"] = transform.validate(
                            st.session_state["pr"],
                            fdr=fdr,  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            formula=formula,  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType] # ty: ignore[possibly-unresolved-reference, invalid-argument-type]
                            score="higher_better"
                            if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            else "lower_better",
                            separate_intra_inter=separate  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            == "Separate FDR for intra and inter matches",
                        )
                    if targets_only:
                        st.session_state["pr"] = transform.targets_only(
                            st.session_state["pr"]
                        )
                        if (
                            "aggregated" in st.session_state
                            and st.session_state["aggregated"] is not None
                        ):
                            st.session_state["aggregated"] = transform.targets_only(
                                st.session_state["aggregated"]
                            )
                    if reannotate_positions:
                        if uploaded_fasta is None:  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            raise ValueError(
                                "Can't annotate crosslink position when no FASTA file is uploaded!"
                            )
                        if not targets_only:
                            _ = st.warning(
                                "Might not be able to (re-)annotate positions if results contain decoy matches!",
                                icon="⚠️",
                            )
                        st.session_state["pr"] = reannotating_positions(
                            st.session_state["pr"],
                            uploaded_fasta,  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                        )
                        if (
                            "aggregated" in st.session_state
                            and st.session_state["aggregated"] is not None
                        ):
                            st.session_state["aggregated"] = reannotating_positions(
                                st.session_state["aggregated"],
                                uploaded_fasta,  # pyright: ignore[reportPossiblyUnboundVariable] # ty: ignore[possibly-unresolved-reference]
                            )
                    st.rerun()
                except Exception as e:
                    # reset meta information
                    if "meta_info" in st.session_state:
                        del st.session_state["meta_info"]
                    # reset pr and aggregated on file read
                    if "pr" in st.session_state:
                        del st.session_state["pr"]
                    if "aggregated" in st.session_state:
                        del st.session_state["aggregated"]
                    # reset any exported files
                    reset_exports()
                    # reset proteins
                    st.session_state["possible_proteins"] = None
                    _ = st.error(
                        "Something went wrong! This is most likely due to missing information in the results!",
                        icon="⚠️",
                    )
                    with st.expander("Show exception"):
                        _ = st.exception(e)

    # display read data and summary [CSMs]
    if "pr" in st.session_state:
        if (
            st.session_state["pr"]["crosslink-spectrum-matches"] is not None
            and len(st.session_state["pr"]["crosslink-spectrum-matches"]) == 0
        ):
            _ = st.error(
                "Filtering criteria too strict! No crosslink-spectrum-matches passed the filter! Please reload your data from the 'Load Data' tab!"
            )
        if (
            st.session_state["pr"]["crosslink-spectrum-matches"] is not None
            and len(st.session_state["pr"]["crosslink-spectrum-matches"]) > 0
        ):
            csms_header = st.subheader(
                "Read Crosslink-Spectrum-Matches", divider="grey"
            )
            csms = st.session_state["pr"]["crosslink-spectrum-matches"]
            csms_info = st.markdown(f"**Read {len(csms)} crosslink-spectrum-matches:**")
            csms_df = st.dataframe(transform.to_dataframe(csms), width="stretch")
            summary_stats = transform.summary(csms)
            summary_stats_md = st.markdown("**Summary Statistics:**")
            summary_stats_df = st.dataframe(
                pd.DataFrame(pd.Series(summary_stats)).T, hide_index=True
            )

            l8, center_8, r8 = st.columns(3)

            with l8:
                csms_dl_csv = st.download_button(
                    label="Download crosslink-spectrum-matches as .csv!",
                    data=dataframe_to_csv_stream(
                        transform.to_dataframe(csms),
                        sep=",",
                        index=False,
                    ),
                    file_name="crosslink-spectrum-matches.csv",
                    on_click="ignore",
                    type="primary",
                    mime="text/csv",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslink-spectrum-matches in comma-separated format.",
                    key="csms_dl_csv",
                )

            with center_8:
                csms_dl_excel = st.download_button(
                    label="Download crosslink-spectrum-matches as .xlsx!",
                    data=dataframe_to_xlsx_stream(
                        transform.to_dataframe(csms),
                        sheet_name="crosslink-spectrum-matches",
                        index=False,
                    ),
                    file_name="crosslink-spectrum-matches.xlsx",
                    on_click="ignore",
                    type="primary",
                    mime="application/vnd.ms-excel",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslink-spectrum-matches in Microsoft Excel format.",
                    key="csms_dl_excel",
                )

            with r8:
                csms_dl_json = st.download_button(
                    label="Download crosslink-spectrum-matches as .json!",
                    data=to_json(csms),
                    file_name="crosslink-spectrum-matches.json",
                    on_click="ignore",
                    type="primary",
                    mime="application/json",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslink-spectrum-matches in JavaScript Object Notation (JSON) format.",
                    key="csms_dl_json",
                )

        # display read data and summary [crosslinks]
        if (
            st.session_state["pr"]["crosslinks"] is not None
            and len(st.session_state["pr"]["crosslinks"]) == 0
        ):
            _ = st.error(
                "Filtering criteria too strict! No crosslinks passed the filter! Please reload your data from the 'Load Data' tab!"
            )
        if (
            st.session_state["pr"]["crosslinks"] is not None
            and len(st.session_state["pr"]["crosslinks"]) > 0
        ):
            crosslinks_header = st.subheader("Read Crosslinks", divider="grey")
            crosslinks = st.session_state["pr"]["crosslinks"]
            crosslinks_info = st.markdown(f"**Read {len(crosslinks)} crosslinks:**")
            crosslinks_df = st.dataframe(
                transform.to_dataframe(crosslinks), width="stretch"
            )
            summary_stats = transform.summary(crosslinks)
            summary_stats_md = st.markdown("**Summary Statistics:**")
            summary_stats_df = st.dataframe(
                pd.DataFrame(pd.Series(summary_stats)).T, hide_index=True
            )

            l9, center_9, r9 = st.columns(3)

            with l9:
                crosslinks_dl_csv = st.download_button(
                    label="Download crosslinks as .csv!",
                    data=dataframe_to_csv_stream(
                        transform.to_dataframe(crosslinks),
                        sep=",",
                        index=False,
                    ),
                    file_name="crosslinks.csv",
                    on_click="ignore",
                    type="primary",
                    mime="text/csv",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslinks in comma-separated format.",
                    key="crosslinks_dl_csv",
                )

            with center_9:
                crosslinks_dl_excel = st.download_button(
                    label="Download crosslinks as .xlsx!",
                    data=dataframe_to_xlsx_stream(
                        transform.to_dataframe(crosslinks),
                        sheet_name="crosslinks",
                        index=False,
                    ),
                    file_name="crosslinks.xlsx",
                    on_click="ignore",
                    type="primary",
                    mime="application/vnd.ms-excel",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslinks in Microsoft Excel format.",
                    key="crosslinks_dl_excel",
                )

            with r9:
                crosslinks_dl_json = st.download_button(
                    label="Download crosslinks as .json!",
                    data=to_json(crosslinks),
                    file_name="crosslinks.json",
                    on_click="ignore",
                    type="primary",
                    mime="application/json",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslinks in JavaScript Object Notation (JSON) format.",
                    key="crosslinks_dl_json",
                )

    # display read data and summary [aggregated crosslinks]
    if (
        "aggregated" in st.session_state
        and st.session_state["aggregated"] is not None
        and len(st.session_state["aggregated"]) == 0
    ):
        _ = st.error(
            "Filtering criteria too strict! None of the aggregated crosslinks passed the filter! Please reload your data from the 'Load Data' tab!"
        )
    if (
        "aggregated" in st.session_state
        and st.session_state["aggregated"] is not None
        and len(st.session_state["aggregated"]) > 0
    ):
        aggregated_crosslinks_header = st.subheader(
            "Aggregated Crosslinks", divider="grey"
        )
        aggregated_crosslinks = st.session_state["aggregated"]
        aggregated_crosslinks_info = st.markdown(
            f"**Aggregated {len(aggregated_crosslinks)} crosslinks:**"
        )
        aggregated_crosslinks_df = st.dataframe(
            transform.to_dataframe(aggregated_crosslinks), width="stretch"
        )
        summary_stats = transform.summary(aggregated_crosslinks)
        summary_stats_md = st.markdown("**Summary Statistics:**")
        summary_stats_df = st.dataframe(
            pd.DataFrame(pd.Series(summary_stats)).T, hide_index=True
        )

        l10, center_10, r10 = st.columns(3)

        with l10:
            aggregated_crosslinks_dl_csv = st.download_button(
                label="Download aggregated crosslinks as .csv!",
                data=dataframe_to_csv_stream(
                    transform.to_dataframe(aggregated_crosslinks),
                    sep=",",
                    index=False,
                ),
                file_name="aggregated_crosslinks.csv",
                on_click="ignore",
                type="primary",
                mime="text/csv",
                icon=":material/download:",
                width="stretch",
                help="Download aggregated crosslinks in comma-separated format.",
                key="aggregated_crosslinks_dl_csv",
            )

        with center_10:
            aggregated_crosslinks_dl_excel = st.download_button(
                label="Download aggregated crosslinks as .xlsx!",
                data=dataframe_to_xlsx_stream(
                    transform.to_dataframe(aggregated_crosslinks),
                    sheet_name="aggregated crosslinks",
                    index=False,
                ),
                file_name="aggregated_crosslinks.xlsx",
                on_click="ignore",
                type="primary",
                mime="application/vnd.ms-excel",
                icon=":material/download:",
                width="stretch",
                help="Download aggregated crosslinks in Microsoft Excel format.",
                key="aggregated_crosslinks_dl_excel",
            )

        with r10:
            aggregated_crosslinks_dl_json = st.download_button(
                label="Download aggregated crosslinks as .json!",
                data=to_json(aggregated_crosslinks),
                file_name="aggregated_crosslinks.json",
                on_click="ignore",
                type="primary",
                mime="application/json",
                icon=":material/download:",
                width="stretch",
                help="Download aggregated crosslinks in JavaScript Object Notation (JSON) format.",
                key="aggregated_crosslinks_dl_json",
            )


# filter tab
def filter_tab():
    if "pr" not in st.session_state and "aggregated" not in st.session_state:
        no_data = st.info("You need to upload at least one result file first!")
    else:
        # filters
        # protein filter
        if (
            "possible_proteins" not in st.session_state
            or st.session_state["possible_proteins"] is None
        ):
            possible_proteins = set()
            _ = st.toast("Loading proteins...", icon="🔄")
            if "pr" in st.session_state and st.session_state["pr"] is not None:
                if st.session_state["pr"]["crosslink-spectrum-matches"] is not None:
                    for csm in st.session_state["pr"]["crosslink-spectrum-matches"]:
                        if csm["alpha_proteins"] is not None:
                            for protein in csm["alpha_proteins"]:
                                possible_proteins.add(protein)
                        if csm["beta_proteins"] is not None:
                            for protein in csm["beta_proteins"]:
                                possible_proteins.add(protein)
                if st.session_state["pr"]["crosslinks"] is not None:
                    for xl in st.session_state["pr"]["crosslinks"]:
                        if xl["alpha_proteins"] is not None:
                            for protein in xl["alpha_proteins"]:
                                possible_proteins.add(protein)
                        if xl["beta_proteins"] is not None:
                            for protein in xl["beta_proteins"]:
                                possible_proteins.add(protein)
            st.session_state["possible_proteins"] = possible_proteins
            _ = st.toast("Successfully loaded proteins!", icon="✅")
        with st.form("filter_form", enter_to_submit=False, border=False):
            protein_filter_header = st.subheader(
                "Filter by Protein Accession", divider="grey"
            )
            protein_filter = st.multiselect(
                "Select the protein accessions that you want to keep:",
                options=st.session_state["possible_proteins"],
                default=None,
                key="protein_filter",
                help="Select the protein accessions that you want to keep. "
                + "Crosslink-spectrum-matches and crosslinks containing none of the proteins will be filtered out. "
                + "Leaving this filter blank will keep everything.",
                max_selections=25,
            )
            # crosslink type filter
            crosslink_type_filter_header = st.subheader(
                "Filter by Crosslink Type", divider="grey"
            )
            crosslink_type_filter = st.multiselect(
                "Select the crosslink types that you want to keep:",
                options=["Intra", "Inter"],
                default=["Intra", "Inter"],
                key="crosslink_type_filter",
                help="Select the crosslink types that you want to keep. "
                + "Crosslink-spectrum-matches and crosslinks that are not of these types will be filtered out.",
                max_selections=2,
            )
            # target decoy filter
            target_decoy_filter_header = st.subheader(
                "Filter by Target-Decoy Type", divider="grey"
            )
            target_decoy_filter = st.multiselect(
                "Select the crosslink types that you want to keep:",
                options=["Target-Target", "Target-Decoy", "Decoy-Decoy"],
                default=[],
                key="target_decoy_filter",
                help="Select the target-decoy types that you want to keep. "
                + "Crosslink-spectrum-matches and crosslinks that are not of these types will be filtered out.",
                max_selections=3,
            )
            # filter action
            left_c, center_c, right_c = st.columns(3)

            with center_c:
                filter_button = st.form_submit_button(
                    "Filter results!", type="primary", width="stretch"
                )

            # filter in all inputs
            if filter_button:
                # reset any exported files
                reset_exports()
                # reset proteins
                st.session_state["possible_proteins"] = None

                with st.spinner("Filtering results...", show_time=True):
                    try:
                        if protein_filter is not None and len(protein_filter) > 0:
                            if "pr" in st.session_state:
                                if (
                                    st.session_state["pr"]["crosslink-spectrum-matches"]
                                    is not None
                                ):
                                    st.session_state["pr"][
                                        "crosslink-spectrum-matches"
                                    ] = filter_proteins(
                                        st.session_state["pr"][
                                            "crosslink-spectrum-matches"
                                        ],
                                        protein_filter,
                                    )
                                if st.session_state["pr"]["crosslinks"] is not None:
                                    st.session_state["pr"]["crosslinks"] = (
                                        filter_proteins(
                                            st.session_state["pr"]["crosslinks"],
                                            protein_filter,
                                        )
                                    )
                            if (
                                "aggregated" in st.session_state
                                and st.session_state["aggregated"] is not None
                            ):
                                st.session_state["aggregated"] = filter_proteins(
                                    st.session_state["aggregated"], protein_filter
                                )
                        if (
                            crosslink_type_filter is not None
                            and len(crosslink_type_filter) > 0
                        ):
                            if "pr" in st.session_state:
                                if (
                                    st.session_state["pr"]["crosslink-spectrum-matches"]
                                    is not None
                                ):
                                    intra_inter = transform.filter_crosslink_type(
                                        st.session_state["pr"][
                                            "crosslink-spectrum-matches"
                                        ]
                                    )
                                    keep = list()
                                    if "Intra" in crosslink_type_filter:
                                        keep += intra_inter["Intra"]
                                    if "Inter" in crosslink_type_filter:
                                        keep += intra_inter["Inter"]
                                    st.session_state["pr"][
                                        "crosslink-spectrum-matches"
                                    ] = keep
                                if st.session_state["pr"]["crosslinks"] is not None:
                                    intra_inter = transform.filter_crosslink_type(
                                        st.session_state["pr"]["crosslinks"]
                                    )
                                    keep = list()
                                    if "Intra" in crosslink_type_filter:
                                        keep += intra_inter["Intra"]
                                    if "Inter" in crosslink_type_filter:
                                        keep += intra_inter["Inter"]
                                    st.session_state["pr"]["crosslinks"] = keep
                            if (
                                "aggregated" in st.session_state
                                and st.session_state["aggregated"] is not None
                            ):
                                intra_inter = transform.filter_crosslink_type(
                                    st.session_state["aggregated"]
                                )
                                keep = list()
                                if "Intra" in crosslink_type_filter:
                                    keep += intra_inter["Intra"]
                                if "Inter" in crosslink_type_filter:
                                    keep += intra_inter["Inter"]
                                st.session_state["aggregated"] = keep
                        if (
                            target_decoy_filter is not None
                            and len(target_decoy_filter) > 0
                        ):
                            if "pr" in st.session_state:
                                if (
                                    st.session_state["pr"]["crosslink-spectrum-matches"]
                                    is not None
                                ):
                                    tt_td_dd = transform.filter_target_decoy(
                                        st.session_state["pr"][
                                            "crosslink-spectrum-matches"
                                        ]
                                    )
                                    keep = list()
                                    if "Target-Target" in target_decoy_filter:
                                        keep += tt_td_dd["Target-Target"]
                                    if "Target-Decoy" in target_decoy_filter:
                                        keep += tt_td_dd["Target-Decoy"]
                                    if "Decoy-Decoy" in target_decoy_filter:
                                        keep += tt_td_dd["Decoy-Decoy"]
                                    st.session_state["pr"][
                                        "crosslink-spectrum-matches"
                                    ] = keep
                                if st.session_state["pr"]["crosslinks"] is not None:
                                    tt_td_dd = transform.filter_target_decoy(
                                        st.session_state["pr"]["crosslinks"]
                                    )
                                    keep = list()
                                    if "Target-Target" in target_decoy_filter:
                                        keep += tt_td_dd["Target-Target"]
                                    if "Target-Decoy" in target_decoy_filter:
                                        keep += tt_td_dd["Target-Decoy"]
                                    if "Decoy-Decoy" in target_decoy_filter:
                                        keep += tt_td_dd["Decoy-Decoy"]
                                    st.session_state["pr"]["crosslinks"] = keep
                            if (
                                "aggregated" in st.session_state
                                and st.session_state["aggregated"] is not None
                            ):
                                tt_td_dd = transform.filter_target_decoy(
                                    st.session_state["aggregated"]
                                )
                                keep = list()
                                if "Target-Target" in target_decoy_filter:
                                    keep += tt_td_dd["Target-Target"]
                                if "Target-Decoy" in target_decoy_filter:
                                    keep += tt_td_dd["Target-Decoy"]
                                if "Decoy-Decoy" in target_decoy_filter:
                                    keep += tt_td_dd["Decoy-Decoy"]
                                st.session_state["aggregated"] = keep
                        st.rerun()
                    except Exception as e:
                        # reset meta information
                        if "meta_info" in st.session_state:
                            del st.session_state["meta_info"]
                        # reset pr and aggregated on file read
                        if "pr" in st.session_state:
                            del st.session_state["pr"]
                        if "aggregated" in st.session_state:
                            del st.session_state["aggregated"]
                        # reset any exported files
                        reset_exports()
                        # reset proteins
                        st.session_state["possible_proteins"] = None
                        _ = st.error(
                            "Something went wrong! This is most likely due to missing information in the results! All results have been reset!",
                            icon="⚠️",
                        )
                        with st.expander("Show exception"):
                            _ = st.exception(e)

    # display filtered data and summary [CSMs]
    if "pr" in st.session_state:
        if (
            st.session_state["pr"]["crosslink-spectrum-matches"] is not None
            and len(st.session_state["pr"]["crosslink-spectrum-matches"]) == 0
        ):
            _ = st.error(
                "Filtering criteria too strict! No crosslink-spectrum-matches passed the filter! Please reload your data from the 'Load Data' tab!"
            )
        if (
            st.session_state["pr"]["crosslink-spectrum-matches"] is not None
            and len(st.session_state["pr"]["crosslink-spectrum-matches"]) > 0
        ):
            csms_header = st.subheader(
                "Current Crosslink-Spectrum-Matches", divider="grey"
            )
            csms = st.session_state["pr"]["crosslink-spectrum-matches"]
            csms_info = st.markdown(
                f"**Currently {len(csms)} crosslink-spectrum-matches:**"
            )
            csms_df = st.dataframe(transform.to_dataframe(csms), width="stretch")
            summary_stats = transform.summary(csms)
            summary_stats_md = st.markdown("**Summary Statistics:**")
            summary_stats_df = st.dataframe(
                pd.DataFrame(pd.Series(summary_stats)).T, hide_index=True
            )

            l1, center_1, r1 = st.columns(3)

            with l1:
                filter_csms_dl_csv = st.download_button(
                    label="Download crosslink-spectrum-matches as .csv!",
                    data=dataframe_to_csv_stream(
                        transform.to_dataframe(csms),
                        sep=",",
                        index=False,
                    ),
                    file_name="crosslink-spectrum-matches.csv",
                    on_click="ignore",
                    type="primary",
                    mime="text/csv",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslink-spectrum-matches in comma-separated format.",
                    key="filter_csms_dl_csv",
                )

            with center_1:
                filter_csms_dl_excel = st.download_button(
                    label="Download crosslink-spectrum-matches as .xlsx!",
                    data=dataframe_to_xlsx_stream(
                        transform.to_dataframe(csms),
                        sheet_name="crosslink-spectrum-matches",
                        index=False,
                    ),
                    file_name="crosslink-spectrum-matches.xlsx",
                    on_click="ignore",
                    type="primary",
                    mime="application/vnd.ms-excel",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslink-spectrum-matches in Microsoft Excel format.",
                    key="filter_csms_dl_excel",
                )

            with r1:
                filter_csms_dl_json = st.download_button(
                    label="Download crosslink-spectrum-matches as .json!",
                    data=to_json(csms),
                    file_name="crosslink-spectrum-matches.json",
                    on_click="ignore",
                    type="primary",
                    mime="application/json",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslink-spectrum-matches in JavaScript Object Notation (JSON) format.",
                    key="filter_csms_dl_json",
                )

        # display filtered data and summary [crosslinks]
        if (
            st.session_state["pr"]["crosslinks"] is not None
            and len(st.session_state["pr"]["crosslinks"]) == 0
        ):
            _ = st.error(
                "Filtering criteria too strict! No crosslinks passed the filter! Please reload your data from the 'Load Data' tab!"
            )
        if (
            st.session_state["pr"]["crosslinks"] is not None
            and len(st.session_state["pr"]["crosslinks"]) > 0
        ):
            crosslinks_header = st.subheader("Current Crosslinks", divider="grey")
            crosslinks = st.session_state["pr"]["crosslinks"]
            crosslinks_info = st.markdown(
                f"**Currently {len(crosslinks)} crosslinks:**"
            )
            crosslinks_df = st.dataframe(
                transform.to_dataframe(crosslinks), width="stretch"
            )
            summary_stats = transform.summary(crosslinks)
            summary_stats_md = st.markdown("**Summary Statistics:**")
            summary_stats_df = st.dataframe(
                pd.DataFrame(pd.Series(summary_stats)).T, hide_index=True
            )

            l2, center_2, r2 = st.columns(3)

            with l2:
                filter_crosslinks_dl_csv = st.download_button(
                    label="Download crosslinks as .csv!",
                    data=dataframe_to_csv_stream(
                        transform.to_dataframe(crosslinks),
                        sep=",",
                        index=False,
                    ),
                    file_name="crosslinks.csv",
                    on_click="ignore",
                    type="primary",
                    mime="text/csv",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslinks in comma-separated format.",
                    key="filter_crosslinks_dl_csv",
                )

            with center_2:
                filter_crosslinks_dl_excel = st.download_button(
                    label="Download crosslinks as .xlsx!",
                    data=dataframe_to_xlsx_stream(
                        transform.to_dataframe(crosslinks),
                        sheet_name="crosslinks",
                        index=False,
                    ),
                    file_name="crosslinks.xlsx",
                    on_click="ignore",
                    type="primary",
                    mime="application/vnd.ms-excel",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslinks in Microsoft Excel format.",
                    key="filter_crosslinks_dl_excel",
                )

            with r2:
                filter_crosslinks_dl_json = st.download_button(
                    label="Download crosslinks as .json!",
                    data=to_json(crosslinks),
                    file_name="crosslinks.json",
                    on_click="ignore",
                    type="primary",
                    mime="application/json",
                    icon=":material/download:",
                    width="stretch",
                    help="Download crosslinks in JavaScript Object Notation (JSON) format.",
                    key="filter_crosslinks_dl_json",
                )

    # display filtered data and summary [aggregated crosslinks]
    if (
        "aggregated" in st.session_state
        and st.session_state["aggregated"] is not None
        and len(st.session_state["aggregated"]) == 0
    ):
        _ = st.error(
            "Filtering criteria too strict! None of the aggregated crosslinks passed the filter! Please reload your data from the 'Load Data' tab!"
        )
    if (
        "aggregated" in st.session_state
        and st.session_state["aggregated"] is not None
        and len(st.session_state["aggregated"]) > 0
    ):
        aggregated_crosslinks_header = st.subheader(
            "Current Aggregated Crosslinks", divider="grey"
        )
        aggregated_crosslinks = st.session_state["aggregated"]
        aggregated_crosslinks_info = st.markdown(
            f"**Currently aggregated {len(aggregated_crosslinks)} crosslinks:**"
        )
        aggregated_crosslinks_df = st.dataframe(
            transform.to_dataframe(aggregated_crosslinks), width="stretch"
        )
        summary_stats = transform.summary(aggregated_crosslinks)
        summary_stats_md = st.markdown("**Summary Statistics:**")
        summary_stats_df = st.dataframe(
            pd.DataFrame(pd.Series(summary_stats)).T, hide_index=True
        )

        l3, center_3, r3 = st.columns(3)

        with l3:
            filter_aggregated_crosslinks_dl_csv = st.download_button(
                label="Download aggregated crosslinks as .csv!",
                data=dataframe_to_csv_stream(
                    transform.to_dataframe(aggregated_crosslinks),
                    sep=",",
                    index=False,
                ),
                file_name="aggregated_crosslinks.csv",
                on_click="ignore",
                type="primary",
                mime="text/csv",
                icon=":material/download:",
                width="stretch",
                help="Download aggregated crosslinks in comma-separated format.",
                key="filter_aggregated_crosslinks_dl_csv",
            )

        with center_3:
            filter_aggregated_crosslinks_dl_excel = st.download_button(
                label="Download aggregated crosslinks as .xlsx!",
                data=dataframe_to_xlsx_stream(
                    transform.to_dataframe(aggregated_crosslinks),
                    sheet_name="aggregated crosslinks",
                    index=False,
                ),
                file_name="aggregated_crosslinks.xlsx",
                on_click="ignore",
                type="primary",
                mime="application/vnd.ms-excel",
                icon=":material/download:",
                width="stretch",
                help="Download aggregated crosslinks in Microsoft Excel format.",
                key="filter_aggregated_crosslinks_dl_excel",
            )

        with r3:
            filter_aggregated_crosslinks_dl_json = st.download_button(
                label="Download aggregated crosslinks as .json!",
                data=to_json(aggregated_crosslinks),
                file_name="aggregated_crosslinks.json",
                on_click="ignore",
                type="primary",
                mime="application/json",
                icon=":material/download:",
                width="stretch",
                help="Download aggregated crosslinks in JavaScript Object Notation (JSON) format.",
                key="filter_aggregated_crosslinks_dl_json",
            )


# visualize tab
def visualize_tab():
    visualization_header = st.subheader("Visualization Parameters", divider="grey")
    bins = st.select_slider(
        "Number of histogram bins:",
        range(5, 105, 5),
        value=25,
        help="The number of histogram bins to use for the score distribution plot.",
    )
    top_n = st.number_input(
        "Maximum number of proteins and peptide pairs to display:",
        min_value=1,
        max_value=None,
        value=25,
        step=1,
        help="Maximum number of proteins and peptide pairs to display. Proteins and peptide pairs are sorted by the number of associated elements.",
    )
    if "pr" not in st.session_state and "aggregated" not in st.session_state:
        no_data = st.info("You need to upload at least one result file first!")
    if "pr" in st.session_state and st.session_state["pr"] is not None:
        if (
            st.session_state["pr"]["crosslink-spectrum-matches"] is not None
            and len(st.session_state["pr"]["crosslink-spectrum-matches"]) == 0
        ):
            _ = st.error(
                "Filtering criteria too strict! No crosslink-spectrum-matches passed the filter! Please reload your data from the 'Load Data' tab!"
            )
        if (
            st.session_state["pr"]["crosslink-spectrum-matches"] is not None
            and len(st.session_state["pr"]["crosslink-spectrum-matches"]) > 0
        ):
            csms = st.session_state["pr"]["crosslink-spectrum-matches"]
            csms_viz_header = st.subheader(
                "Visualizations for Crosslink-Spectrum-Matches", divider="grey"
            )
            available_keys = transform.get_available_keys(csms)
            plots = list()
            # target decoy distribution
            if available_keys["alpha_decoy"] and available_keys["beta_decoy"]:
                fig, ax = plotting.plot_target_decoy_distribution(
                    csms, figsize=(8.0, 4.5)
                )
                plots.append(fig)
            # score distribution
            if (
                available_keys["score"]
                and available_keys["alpha_decoy"]
                and available_keys["beta_decoy"]
            ):
                fig, ax = plotting.plot_score_distribution(
                    csms, bins=bins, figsize=(8.0, 4.5)
                )
                plots.append(fig)
            # crosslink type distribution plot is always possible
            fig, ax = plotting.plot_crosslink_type_distribution(
                csms, plot_type="bar", figsize=(8.0, 4.5)
            )
            plots.append(fig)
            # protein distribution
            if available_keys["alpha_proteins"] and available_keys["beta_proteins"]:
                fig, ax = plotting.plot_protein_distribution(
                    csms, top_n=top_n, figsize=(8.0, 4.5)
                )
                plots.append(fig)
            # peptide pair distribution plot is always possible
            fig, ax = plotting.plot_peptide_pair_distribution(
                csms, top_n=top_n, figsize=(8.0, 4.5)
            )
            plots.append(fig)
            if len(plots) > 0:
                layout_plots(plots)
            else:
                # technically impossible
                csms_not_enough_data = st.info(
                    "Not enough data to plot anything for crosslink-spectrum-matches!"
                )
        if (
            st.session_state["pr"]["crosslinks"] is not None
            and len(st.session_state["pr"]["crosslinks"]) == 0
        ):
            _ = st.error(
                "Filtering criteria too strict! No crosslinks passed the filter! Please reload your data from the 'Load Data' tab!"
            )
        if (
            st.session_state["pr"]["crosslinks"] is not None
            and len(st.session_state["pr"]["crosslinks"]) > 0
        ):
            crosslinks = st.session_state["pr"]["crosslinks"]
            crosslinks_viz_header = st.subheader(
                "Visualizations for Crosslinks", divider="grey"
            )
            available_keys = transform.get_available_keys(crosslinks)
            plots = list()
            # target decoy distribution
            if available_keys["alpha_decoy"] and available_keys["beta_decoy"]:
                fig, ax = plotting.plot_target_decoy_distribution(
                    crosslinks, figsize=(8.0, 4.5)
                )
                plots.append(fig)
            # score distribution
            if (
                available_keys["score"]
                and available_keys["alpha_decoy"]
                and available_keys["beta_decoy"]
            ):
                fig, ax = plotting.plot_score_distribution(
                    crosslinks, bins=bins, figsize=(8.0, 4.5)
                )
                plots.append(fig)
            # crosslink type distribution plot is always possible
            fig, ax = plotting.plot_crosslink_type_distribution(
                crosslinks, plot_type="bar", figsize=(8.0, 4.5)
            )
            plots.append(fig)
            # protein distribution
            if available_keys["alpha_proteins"] and available_keys["beta_proteins"]:
                fig, ax = plotting.plot_protein_distribution(
                    crosslinks, top_n=top_n, figsize=(8.0, 4.5)
                )
                plots.append(fig)
            if len(plots) > 0:
                layout_plots(plots)
            else:
                # technically impossible
                crosslinks_not_enough_data = st.info(
                    "Not enough data to plot anything for crosslinks!"
                )
    if (
        "aggregated" in st.session_state
        and st.session_state["aggregated"] is not None
        and len(st.session_state["aggregated"]) == 0
    ):
        _ = st.error(
            "Filtering criteria too strict! None of the aggregated crosslinks passed the filter! Please reload your data from the 'Load Data' tab!"
        )
    if (
        "aggregated" in st.session_state
        and st.session_state["aggregated"] is not None
        and len(st.session_state["aggregated"]) > 0
    ):
        aggregated_crosslinks = st.session_state["aggregated"]
        aggregated_crosslinks_viz_header = st.subheader(
            "Visualizations for Aggregated Crosslinks", divider="grey"
        )
        available_keys = transform.get_available_keys(aggregated_crosslinks)
        plots = list()
        # target decoy distribution
        if available_keys["alpha_decoy"] and available_keys["beta_decoy"]:
            fig, ax = plotting.plot_target_decoy_distribution(
                aggregated_crosslinks, figsize=(8.0, 4.5)
            )
            plots.append(fig)
        # score distribution
        if (
            available_keys["score"]
            and available_keys["alpha_decoy"]
            and available_keys["beta_decoy"]
        ):
            fig, ax = plotting.plot_score_distribution(
                aggregated_crosslinks, bins=bins, figsize=(8.0, 4.5)
            )
            plots.append(fig)
        # crosslink type distribution plot is always possible
        fig, ax = plotting.plot_crosslink_type_distribution(
            aggregated_crosslinks, plot_type="bar", figsize=(8.0, 4.5)
        )
        plots.append(fig)
        # protein distribution
        if available_keys["alpha_proteins"] and available_keys["beta_proteins"]:
            fig, ax = plotting.plot_protein_distribution(
                aggregated_crosslinks, top_n=top_n, figsize=(8.0, 4.5)
            )
            plots.append(fig)
        if len(plots) > 0:
            layout_plots(plots)
        else:
            # technically impossible
            aggregated_crosslinks_not_enough_data = st.info(
                "Not enough data to plot anything for aggregated crosslinks!"
            )


# export tab
def export_tab():
    if "pr" not in st.session_state and "aggregated" not in st.session_state:
        no_data = st.info("You need to upload at least one result file first!")
    if "pr" in st.session_state and st.session_state["pr"] is not None:
        # exporting CSMs
        if (
            st.session_state["pr"]["crosslink-spectrum-matches"] is not None
            and len(st.session_state["pr"]["crosslink-spectrum-matches"]) == 0
        ):
            _ = st.error(
                "Filtering criteria too strict! No crosslink-spectrum-matches passed the filter! Please reload your data from the 'Load Data' tab!"
            )
        if (
            st.session_state["pr"]["crosslink-spectrum-matches"] is not None
            and len(st.session_state["pr"]["crosslink-spectrum-matches"]) > 0
        ):
            csms = st.session_state["pr"]["crosslink-spectrum-matches"]
            export_csms_header = st.subheader(
                "Export Crosslink-Spectrum-Matches", divider="grey"
            )
            export_csms_options = ["IMP-X-FDR", "MS Annika", "ProXL", "xiFDR"]
            export_csms_picker = st.selectbox(
                "Export crosslink-spectrum-matches to:",
                options=export_csms_options,
                index=None,
                help="Choose a format to export the crosslink-spectrum-matches to.",
            )
            if export_csms_picker is None:
                pass
            # IMP-X-FDR
            elif export_csms_picker == "IMP-X-FDR":
                export_csms_impxfdr_info = st.info(
                    "To export to IMP-X-FDR your crosslink-spectrum-matches should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence matches "
                    + "to compare FDR estimation to the experimentally validated FDR! You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                export_csms_impxfdr_button = st.button(
                    "Export to IMP-X-FDR!",
                    type="primary",
                    width="stretch",
                    key="export_csms_impxfdr_button",
                )
                if export_csms_impxfdr_button:
                    with st.spinner(
                        "Exporting crosslink-spectrum-matches to IMP-X-FDR...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_csms_impxfdr"] = (
                                exporter.to_impxfdr(csms, filename=None)
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_csms_impxfdr" in st.session_state
                    and st.session_state["export_csms_impxfdr"] is not None
                ):
                    export_csms_impxfdr_download_info = st.markdown(
                        "Your exported crosslink-spectrum-matches in IMP-X-FDR format are ready for download:"
                    )
                    export_csms_impxfdr_download = st.download_button(
                        label="Download in IMP-X-FDR format!",
                        data=dataframe_to_xlsx_stream(
                            st.session_state["export_csms_impxfdr"],
                            sheet_name="imp-x-fdr",
                            index=False,
                        ),
                        file_name="crosslink-spectrum-matches_imp-x-fdr.xlsx",
                        on_click="ignore",
                        type="primary",
                        mime="application/vnd.ms-excel",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslink-spectrum-matches in IMP-X-FDR format.",
                        key="export_csms_impxfdr_download",
                    )
                    export_csms_impxfdr_download_goto_tool = st.link_button(
                        "Go to IMP-X-FDR!",
                        url="https://github.com/vbc-proteomics-org/imp-x-fdr",
                        help="Go to the IMP-X-FDR download page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # MS Annika
            elif export_csms_picker == "MS Annika":
                export_csms_msannika_button = st.button(
                    "Export to MS Annika format!",
                    type="primary",
                    width="stretch",
                    key="export_csms_msannika_button",
                )
                if export_csms_msannika_button:
                    with st.spinner(
                        "Exporting crosslink-spectrum-matches to MS Annika format...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_csms_msannika"] = (
                                exporter.to_msannika(csms, filename=None)
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_csms_msannika" in st.session_state
                    and st.session_state["export_csms_msannika"] is not None
                ):
                    export_csms_msannika_download_info = st.markdown(
                        "Your exported crosslink-spectrum-matches in MS Annika format are ready for download:"
                    )
                    l_csms_msannika, r_csms_msannika = st.columns(2)
                    with l_csms_msannika:
                        export_csms_msannika_download_csv = st.download_button(
                            label="Download in MS Annika .csv format!",
                            data=dataframe_to_csv_stream(
                                st.session_state["export_csms_msannika"],
                                sep=",",
                                index=False,
                            ),
                            file_name="crosslink-spectrum-matches_ms-annika.csv",
                            on_click="ignore",
                            type="primary",
                            mime="text/csv",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the exported crosslink-spectrum-matches in MS Annika comma-separated-values (.csv) format.",
                            key="export_csms_msannika_download_csv",
                        )
                    with r_csms_msannika:
                        export_csms_msannika_download_xlsx = st.download_button(
                            label="Download in MS Annika .xlsx format!",
                            data=dataframe_to_xlsx_stream(
                                st.session_state["export_csms_msannika"],
                                sheet_name="msannika",
                                index=False,
                            ),
                            file_name="crosslink-spectrum-matches_ms-annika.xlsx",
                            on_click="ignore",
                            type="primary",
                            mime="application/vnd.ms-excel",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the exported crosslink-spectrum-matches in MS Annika Microsoft Excel (.xlsx) format.",
                            key="export_csms_msannika_download_xlsx",
                        )
            # ProXL
            elif export_csms_picker == "ProXL":
                export_csms_proxl_info = st.info(
                    "To export to ProXL your crosslink-spectrum-matches should be **unique** and **should not** "
                    + "**contain decoy matches**! It is also required that all crosslink-spectrum-matches have an "
                    + "associated score and charge! It is **not necessary to check this preemptively** as the exporter "
                    + "automatically checks that this information is available and will throw an error otherwise! "
                    + "Please however **make sure** that your crosslink-spectrum-matches are **unique** and **do not contain decoys** "
                    + "as otherwise the export to ProXL or ProXL itself will not work as intended! You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                with st.form(
                    "export_csms_proxl_form", enter_to_submit=False, border=False
                ):
                    export_csms_proxl_fasta = st.file_uploader(
                        "Upload the FASTA file containing the protein sequences of your crosslink-spectrum-matches:",
                        type="fasta",
                        accept_multiple_files=False,
                        key="export_csms_proxl_fasta",
                        help="Upload the FASTA file containing protein sequences for the provided crosslink spectrum matches.",
                    )
                    export_csms_proxl_search_engine = st.text_input(
                        "Name of the used crosslink search engine [this field has been pre-filled with your selection from the 'Load Data' tab]:",
                        value=st.session_state["meta_info"]["search_engine"],
                        max_chars=150,
                        placeholder=None,
                        key="export_csms_proxl_search_engine",
                        help="Name of the crosslink search engine used in the experiment of the uploaded result file.",
                    )
                    export_csms_proxl_search_engine_version = st.text_input(
                        "Software version of the used crosslink search engine:",
                        value=None,
                        max_chars=150,
                        placeholder="v1.0.0",
                        key="export_csms_proxl_search_engine_version",
                        help="Name of the crosslink search engine used in the experiment of the uploaded result file.",
                    )
                    export_csms_proxl_score = st.selectbox(
                        "Is a higher crosslink-spectrum-match score considered better?",
                        options=["Higher better", "Lower better"],
                        index=0,
                        key="export_csms_proxl_score",
                        help="If a higher crosslink-spectrum-match score is considered better, or a lower score is considered better.",
                    )
                    export_csms_proxl_crosslinker_name = st.text_input(
                        "Name of the used crosslinker [this field has been pre-filled with your selection from the 'Load Data' tab]:",
                        value=st.session_state["meta_info"]["crosslinker_name"],
                        max_chars=50,
                        placeholder="DSSO",
                        key="export_csms_proxl_crosslinker_name",
                        help="Name of the crosslinker used in the experiment of the uploaded result file.",
                    )
                    export_csms_proxl_crosslinker_mass = st.number_input(
                        "Mass of the used crosslinker [this field has been pre-filled with your selection from the 'Load Data' tab]:",
                        value=st.session_state["meta_info"]["crosslinker_mass"],
                        step=0.00001,
                        format="%0.5f",
                        placeholder="158.00376",
                        key="export_csms_proxl_crosslinker_mass",
                        help="Monoisotopic delta mass of the crosslinker used in the experiment of the uploaded result file.",
                    )
                    export_csms_proxl_button = st.form_submit_button(
                        "Export to ProXL!",
                        type="primary",
                        width="stretch",
                        key="export_csms_proxl_button",
                    )
                    if export_csms_proxl_button:
                        if export_csms_proxl_fasta is None:
                            _ = st.error("You need to upload a FASTA file first!")
                        if (
                            export_csms_proxl_search_engine is None
                            or export_csms_proxl_search_engine.strip() == ""
                        ):
                            _ = st.error(
                                "You need to specify the name of the crosslink search engine first!"
                            )
                        if (
                            export_csms_proxl_search_engine_version is None
                            or export_csms_proxl_search_engine_version.strip() == ""
                        ):
                            _ = st.error(
                                "You need to specify the version of the crosslink search engine first!"
                            )
                        if (
                            export_csms_proxl_crosslinker_name is None
                            or export_csms_proxl_crosslinker_name.strip() == ""
                        ):
                            _ = st.error(
                                "You need the specify the name of the used crosslink reagent first!"
                            )
                        if export_csms_proxl_crosslinker_mass is None:
                            _ = st.error(
                                "You need to specify the delta mass of the used crosslink reagent first!"
                            )
                        if (
                            export_csms_proxl_fasta is not None
                            and export_csms_proxl_search_engine is not None
                            and export_csms_proxl_search_engine.strip() != ""
                            and export_csms_proxl_search_engine_version is not None
                            and export_csms_proxl_search_engine_version.strip() != ""
                            and export_csms_proxl_crosslinker_name is not None
                            and export_csms_proxl_crosslinker_name.strip() != ""
                            and export_csms_proxl_crosslinker_mass is not None
                        ):
                            with st.spinner(
                                "Exporting crosslink-spectrum-matches to ProXL...",
                                show_time=True,
                            ):
                                try:
                                    st.session_state["export_csms_proxl"] = (
                                        export_proxl(
                                            csms,
                                            export_csms_proxl_fasta,
                                            export_csms_proxl_search_engine,
                                            export_csms_proxl_search_engine_version,
                                            "higher_better"
                                            if export_csms_proxl_score
                                            == "Higher better"
                                            else "lower_better",
                                            export_csms_proxl_crosslinker_name,
                                            export_csms_proxl_crosslinker_mass,
                                        )
                                    )
                                except Exception as e:
                                    _ = st.error(
                                        "Something went wrong! This is most likely due to missing information in the results!",
                                        icon="⚠️",
                                    )
                                    with st.expander("Show exception"):
                                        _ = st.exception(e)
                if (
                    "export_csms_proxl" in st.session_state
                    and st.session_state["export_csms_proxl"] is not None
                ):
                    export_csms_proxl_download_info = st.markdown(
                        "Your exported crosslink-spectrum-matches in ProXL format are ready for download:"
                    )
                    export_csms_proxl_download = st.download_button(
                        label="Download in ProXL format!",
                        data=to_text(st.session_state["export_csms_proxl"]),
                        file_name="crosslink-spectrum-matches_proxl.xml",
                        on_click="ignore",
                        type="primary",
                        mime="application/xml",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslink-spectrum-matches in ProXL format.",
                        key="export_csms_proxl_download",
                    )
                    export_csms_proxl_download_goto_tool = st.link_button(
                        "Go to ProXL!",
                        url="https://www.yeastrc.org/proxl_public/",
                        help="Go to the ProXL page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # xiFDR
            elif export_csms_picker == "xiFDR":
                export_csms_xifdr_info = st.info(
                    "To export to xiFDR your crosslink-spectrum-matches should be **unique** and **must** "
                    + "**contain decoy matches**! It is also required that all crosslink-spectrum-matches have "
                    + "associated proteins for the alpha and beta peptide as well as the corresponding crosslink "
                    + "positions in those proteins! Additionally, all crosslink-spectrum-matches need to specify "
                    + "if the alpha and beta peptide are target or decoy matches (separately) and need to have an "
                    + "associated score and charge! It is **not necessary to check this preemptively** as the exporter "
                    + "automatically checks that this information is available and will throw an error otherwise! "
                    + "Please however **make sure** that your crosslink-spectrum-matches are **unique** and **contain decoys** "
                    + "as otherwise xiFDR will not work as intended! You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                export_csms_xifdr_button = st.button(
                    "Export to xiFDR!",
                    type="primary",
                    width="stretch",
                    key="export_csms_xifdr_button",
                )
                if export_csms_xifdr_button:
                    with st.spinner(
                        "Exporting crosslink-spectrum-matches to xiFDR...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_csms_xifdr"] = exporter.to_xifdr(
                                csms, filename=None
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_csms_xifdr" in st.session_state
                    and st.session_state["export_csms_xifdr"] is not None
                ):
                    export_csms_xifdr_download_info = st.markdown(
                        "Your exported crosslink-spectrum-matches in xiFDR format are ready for download:"
                    )
                    export_csms_xifdr_download = st.download_button(
                        label="Download in xiFDR format!",
                        data=dataframe_to_csv_stream(
                            st.session_state["export_csms_xifdr"],
                            sep=",",
                            index=False,
                        ),
                        file_name="crosslink-spectrum-matches_xifdr.csv",
                        on_click="ignore",
                        type="primary",
                        mime="text/csv",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslink-spectrum-matches in xiFDR format.",
                        key="export_csms_xifdr_download",
                    )
                    export_csms_xifdr_download_goto_tool = st.link_button(
                        "Go to xiFDR!",
                        url="https://www.rappsilberlab.org/software/xifdr/",
                        help="Go to the xiFDR download page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            else:
                pass

        # exporting crosslinks
        if (
            st.session_state["pr"]["crosslinks"] is not None
            and len(st.session_state["pr"]["crosslinks"]) == 0
        ):
            _ = st.error(
                "Filtering criteria too strict! No crosslinks passed the filter! Please reload your data from the 'Load Data' tab!"
            )
        if (
            st.session_state["pr"]["crosslinks"] is not None
            and len(st.session_state["pr"]["crosslinks"]) > 0
        ):
            crosslinks = st.session_state["pr"]["crosslinks"]
            export_crosslinks_header = st.subheader("Export Crosslinks", divider="grey")
            export_crosslinks_options = [
                "AlphaLink2",
                "IMP-X-FDR",
                "MS Annika",
                "PyXlinkViewer",
                "xiNET",
                "xiVIEW",
                "XlinkDB",
                "xlms-tools",
                "XMAS",
            ]
            export_crosslinks_picker = st.selectbox(
                "Export crosslinks to:",
                options=export_crosslinks_options,
                index=None,
                help="Choose a format to export the crosslinks to.",
            )
            if export_crosslinks_picker is None:
                pass
            # AlphaLink2
            elif export_crosslinks_picker == "AlphaLink2":
                export_crosslinks_alphalink2_info = st.info(
                    "To export to AlphaLink2 your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                    + "You should also filter your crosslinks to only contain residue pairs of your protein(s) of interest! "
                    + "You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                with st.form(
                    "export_crosslinks_alphalink2_form",
                    enter_to_submit=False,
                    border=False,
                ):
                    crosslinks_alphalink2_fasta_file = st.file_uploader(
                        "Upload a FASTA file of proteins/chains of interest:",
                        type="fasta",
                        accept_multiple_files=False,
                        key="crosslinks_alphalink2_fasta_file",
                        help="Upload a FASTA file containing protein/chain sequences. Please keep in mind that AlphaLink2 supports a maximum of 62 proteins/chains!",
                    )
                    crosslinks_alphalink2_annotated_fdr = st.number_input(
                        "Annotated FDR:",
                        value=0.01,
                        min_value=0.0,
                        max_value=1.0,
                        step=0.001,
                        format="%0.3f",
                        key="crosslinks_alphalink2_annotated_fdr",
                        help="Value to use for the 'FDR' column in the AlphaLink2 crosslink table, must be given as a real number between 0 and 1. The default of 0.01 corresponds to 1% FDR.",
                    )
                    export_crosslinks_alphalink2_button = st.form_submit_button(
                        "Export to AlphaLink2 format!",
                        type="primary",
                        width="stretch",
                        key="export_crosslinks_alphalink2_button",
                    )
                    if export_crosslinks_alphalink2_button:
                        if (
                            crosslinks_alphalink2_fasta_file is None
                            or crosslinks_alphalink2_annotated_fdr is None
                        ):
                            _ = st.error(
                                "Can't export to AlphaLink2 when either FASTA file or annotated FDR are missing!",
                                icon="⚠️",
                            )
                        else:
                            with st.spinner(
                                "Exporting crosslinks to AlphaLink2 format...",
                                show_time=True,
                            ):
                                try:
                                    st.session_state["export_crosslinks_alphalink2"] = (
                                        export_alphalink2(
                                            crosslinks,
                                            crosslinks_alphalink2_fasta_file,
                                            float(crosslinks_alphalink2_annotated_fdr),
                                        )
                                    )
                                except Exception as e:
                                    _ = st.error(
                                        "Something went wrong! This is most likely due to missing information in the results!",
                                        icon="⚠️",
                                    )
                                    with st.expander("Show exception"):
                                        _ = st.exception(e)
                if (
                    "export_crosslinks_alphalink2" in st.session_state
                    and st.session_state["export_crosslinks_alphalink2"] is not None
                ):
                    export_crosslinks_alphalink2_download_info = st.markdown(
                        "Your exported crosslinks in AlphaLink2 format are ready for download:"
                    )
                    (
                        export_crosslinks_alphalink2_download_l,
                        export_crosslinks_alphalink2_download_m,
                        export_crosslinks_alphalink2_download_r,
                    ) = st.columns(3)

                    with export_crosslinks_alphalink2_download_l:
                        export_crosslinks_alphalink2_download_txt = st.download_button(
                            label="Download crosslinks in AlphaLink2 format!",
                            data=to_text(
                                st.session_state["export_crosslinks_alphalink2"][
                                    "AlphaLink2 crosslinks"
                                ]
                            ),
                            file_name="crosslinks_AlphaLink2.txt",
                            on_click="ignore",
                            type="primary",
                            mime="text/plain",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the exported crosslinks in AlphaLink2 format.",
                            key="export_crosslinks_alphalink2_download_txt",
                        )
                    with export_crosslinks_alphalink2_download_m:
                        export_crosslinks_alphalink2_download_fasta = st.download_button(
                            label="Download FASTA in AlphaLink2 format!",
                            data=to_text(
                                st.session_state["export_crosslinks_alphalink2"][
                                    "AlphaLink2 FASTA"
                                ]
                            ),
                            file_name="crosslinks_AlphaLink2.fasta",
                            on_click="ignore",
                            type="primary",
                            mime="chemical/seq-aa-fasta",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the uploaded FASTA file in AlphaLink2 format.",
                            key="export_crosslinks_alphalink2_download_fasta",
                        )
                    with export_crosslinks_alphalink2_download_r:
                        export_crosslinks_alphalink2_download_pickle = st.download_button(
                            label="Download pickled crosslinks in AlphaLink2 format!",
                            data=pickle_and_gzip(
                                st.session_state["export_crosslinks_alphalink2"][
                                    "AlphaLink2 Pickle"
                                ]
                            ),
                            file_name="crosslinks_AlphaLink2.pickle.gz",
                            on_click="ignore",
                            type="primary",
                            mime="application/gzip",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the exported crosslinks in pickled AlphaLink2 format.",
                            key="export_crosslinks_alphalink2_download_pickle",
                        )
                    with st.expander("Show Exported Crosslinks"):
                        export_crosslinks_alphalink2_df_info = st.markdown(
                            "**Number of mapped residue pairs:** "
                            + f"{st.session_state['export_crosslinks_alphalink2']['AlphaLink2 DataFrame'].shape[0]}"
                        )
                        export_crosslinks_alphalink2_df = st.dataframe(
                            st.session_state["export_crosslinks_alphalink2"][
                                "AlphaLink2 DataFrame"
                            ],
                            width="stretch",
                        )
                    with st.expander("Show Exported FASTA"):
                        export_crosslinks_alphalink2_nr_proteins = len(
                            [
                                line
                                for line in st.session_state[
                                    "export_crosslinks_alphalink2"
                                ]["AlphaLink2 FASTA"].split("\n")
                                if line.startswith(">")
                            ]
                        )
                        export_crosslinks_alphalink2_display_fasta_info = st.markdown(
                            "**Number of mapped proteins/chains:** "
                            + f"{export_crosslinks_alphalink2_nr_proteins}"
                        )
                        export_crosslinks_alphalink2_display_fasta = st.text(
                            st.session_state["export_crosslinks_alphalink2"][
                                "AlphaLink2 FASTA"
                            ],
                            width="stretch",
                        )

                    export_crosslinks_alphalink2_download_goto_tool = st.link_button(
                        "Go to AlphaLink2!",
                        url="https://github.com/Rappsilber-Laboratory/AlphaLink2",
                        help="Go to the AlphaLink2 page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # IMP-X-FDR
            elif export_crosslinks_picker == "IMP-X-FDR":
                export_crosslinks_impxfdr_info = st.info(
                    "To export to IMP-X-FDR your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence matches "
                    + "to compare FDR estimation to the experimentally validated FDR! You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                export_crosslinks_impxfdr_button = st.button(
                    "Export to IMP-X-FDR!",
                    type="primary",
                    width="stretch",
                    key="export_crosslinks_impxfdr_button",
                )
                if export_crosslinks_impxfdr_button:
                    with st.spinner(
                        "Exporting crosslinks to IMP-X-FDR...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_crosslinks_impxfdr"] = (
                                exporter.to_impxfdr(crosslinks, filename=None)
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_crosslinks_impxfdr" in st.session_state
                    and st.session_state["export_crosslinks_impxfdr"] is not None
                ):
                    export_crosslinks_impxfdr_download_info = st.markdown(
                        "Your exported crosslinks in IMP-X-FDR format are ready for download:"
                    )
                    export_crosslinks_impxfdr_download = st.download_button(
                        label="Download in IMP-X-FDR format!",
                        data=dataframe_to_xlsx_stream(
                            st.session_state["export_crosslinks_impxfdr"],
                            sheet_name="imp-x-fdr",
                            index=False,
                        ),
                        file_name="crosslinks_imp-x-fdr.xlsx",
                        on_click="ignore",
                        type="primary",
                        mime="application/vnd.ms-excel",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in IMP-X-FDR format.",
                        key="export_crosslinks_impxfdr_download",
                    )
                    export_crosslinks_impxfdr_download_goto_tool = st.link_button(
                        "Go to IMP-X-FDR!",
                        url="https://github.com/vbc-proteomics-org/imp-x-fdr",
                        help="Go to the IMP-X-FDR download page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # MS Annika
            elif export_crosslinks_picker == "MS Annika":
                export_crosslinks_msannika_button = st.button(
                    "Export to MS Annika format!",
                    type="primary",
                    width="stretch",
                    key="export_crosslinks_msannika_button",
                )
                if export_crosslinks_msannika_button:
                    with st.spinner(
                        "Exporting crosslinks to MS Annika format...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_crosslinks_msannika"] = (
                                exporter.to_msannika(crosslinks, filename=None)
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_crosslinks_msannika" in st.session_state
                    and st.session_state["export_crosslinks_msannika"] is not None
                ):
                    export_crosslinks_msannika_download_info = st.markdown(
                        "Your exported crosslinks in MS Annika format are ready for download:"
                    )
                    l_crosslinks_msannika, r_crosslinks_msannika = st.columns(2)
                    with l_crosslinks_msannika:
                        export_crosslinks_msannika_download_csv = st.download_button(
                            label="Download in MS Annika .csv format!",
                            data=dataframe_to_csv_stream(
                                st.session_state["export_crosslinks_msannika"],
                                sep=",",
                                index=False,
                            ),
                            file_name="crosslinks_ms-annika.csv",
                            on_click="ignore",
                            type="primary",
                            mime="text/csv",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the exported crosslinks in MS Annika comma-separated-values (.csv) format.",
                            key="export_crosslinks_msannika_download_csv",
                        )
                    with r_crosslinks_msannika:
                        export_crosslinks_msannika_download_xlsx = st.download_button(
                            label="Download in MS Annika .xlsx format!",
                            data=dataframe_to_xlsx_stream(
                                st.session_state["export_crosslinks_msannika"],
                                sheet_name="msannika",
                                index=False,
                            ),
                            file_name="crosslinks_ms-annika.xlsx",
                            on_click="ignore",
                            type="primary",
                            mime="application/vnd.ms-excel",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the exported crosslinks in MS Annika Microsoft Excel (.xlsx) format.",
                            key="export_crosslinks_msannika_download_xlsx",
                        )
            # PyXlinkViewer
            elif export_crosslinks_picker == "PyXlinkViewer":
                export_crosslinks_pyxlinkviewer_info = st.info(
                    "To export to PyXlinkViewer your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                    + "You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                with st.form(
                    "export_crosslinks_pyxlinkviewer_form",
                    enter_to_submit=False,
                    border=False,
                ):
                    crosslinks_pdb_code = st.text_input(
                        "Specify the PDB identification code of your protein(-complex) of interest:",
                        value=None,
                        max_chars=4,
                        key="crosslinks_pdb_code",
                        help="Specify a 4-letter PDB identification code of your cross-linked protein(-complex) of interest.",
                    )
                    crosslinks_pdb_file = st.file_uploader(
                        "Alternatively, upload a PDB file of your protein(-complex) of interest:",
                        type="pdb",
                        accept_multiple_files=False,
                        key="crosslinks_pdb_file",
                        help="Upload a PDB file of your cross-linked protein(-complex) of interest.",
                    )
                    export_crosslinks_pyxlinkviewer_button = st.form_submit_button(
                        "Export to PyXlinkViewer format!",
                        type="primary",
                        width="stretch",
                        key="export_crosslinks_pyxlinkviewer_button",
                    )
                    if export_crosslinks_pyxlinkviewer_button:
                        if crosslinks_pdb_code is None and crosslinks_pdb_file is None:
                            _ = st.error(
                                "Can't export to PyXlinkViewer when neither PDB code nor file are given!",
                                icon="⚠️",
                            )
                        else:
                            with st.spinner(
                                "Exporting crosslinks to PyXlinkViewer format...",
                                show_time=True,
                            ):
                                try:
                                    if crosslinks_pdb_file is not None:
                                        st.session_state[
                                            "export_crosslinks_pyxlinkviewer"
                                        ] = export_pyxlinkviewer_using_pdbfile(
                                            crosslinks, crosslinks_pdb_file
                                        )
                                    else:
                                        if crosslinks_pdb_code is not None:
                                            if len(crosslinks_pdb_code.strip()) != 4:
                                                raise ValueError(
                                                    "Specified PDB code is not a valid 4-letter PDB identification code!"
                                                )
                                            st.session_state[
                                                "export_crosslinks_pyxlinkviewer"
                                            ] = exporter.to_pyxlinkviewer(
                                                crosslinks,
                                                crosslinks_pdb_code.strip(),
                                                filename_prefix=None,
                                            )
                                        else:
                                            raise RuntimeError(
                                                "Can't export to PyXlinkViewer when neither PDB code nor file are given!"
                                            )
                                except Exception as e:
                                    _ = st.error(
                                        "Something went wrong! This is most likely due to missing information in the results!",
                                        icon="⚠️",
                                    )
                                    with st.expander("Show exception"):
                                        _ = st.exception(e)
                if (
                    "export_crosslinks_pyxlinkviewer" in st.session_state
                    and st.session_state["export_crosslinks_pyxlinkviewer"] is not None
                ):
                    export_crosslinks_pyxlinkviewer_download_info = st.markdown(
                        "Your exported crosslinks in PyXlinkViewer format are ready for download:"
                    )
                    export_crosslinks_pyxlinkviewer_download = st.download_button(
                        label="Download in PyXlinkViewer format!",
                        data=to_text(
                            st.session_state["export_crosslinks_pyxlinkviewer"][
                                "PyXlinkViewer"
                            ]
                        ),
                        file_name="crosslinks_pyxlinkviewer.txt",
                        on_click="ignore",
                        type="primary",
                        mime="text/plain",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in PyXlinkViewer format.",
                        key="export_crosslinks_pyxlinkviewer_download",
                    )
                    with st.expander("Download Meta-data"):
                        export_crosslinks_pyxlinkviewer_download_meta_nr_xl = st.markdown(
                            "**Number of mapped crosslinks:** "
                            + f"{st.session_state['export_crosslinks_pyxlinkviewer']['Number of mapped crosslinks']}"
                        )
                        export_crosslinks_pyxlinkviewer_download_meta_mapping = st.download_button(
                            label="Download crosslink mapping!",
                            data=to_text(
                                st.session_state["export_crosslinks_pyxlinkviewer"][
                                    "Mapping"
                                ]
                            ),
                            file_name="crosslinks_pyxlinkviewer_mapping.txt",
                            on_click="ignore",
                            type="primary",
                            mime="text/plain",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the mapping of crosslinks to the PDB structure.",
                            key="export_crosslinks_pyxlinkviewer_download_meta_mapping",
                        )
                        export_crosslinks_pyxlinkviewer_download_meta_pdb_sequence = st.download_button(
                            label="Download parsed PDB sequence!",
                            data=to_text(
                                pyxlinkviewer_get_fasta(
                                    st.session_state["export_crosslinks_pyxlinkviewer"][
                                        "Parsed PDB sequence"
                                    ]
                                )
                            ),
                            file_name="crosslinks_pyxlinkviewer_pdb_sequence.fasta",
                            on_click="ignore",
                            type="primary",
                            mime="chemical/seq-aa-fasta",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the parsed PDB sequence.",
                            key="export_crosslinks_pyxlinkviewer_download_meta_pdb_sequence",
                        )
                        export_crosslinks_pyxlinkviewer_download_meta_pdb_annotation = st.download_button(
                            label="Download parsed PDB annotation!",
                            data=dataframe_to_csv_stream(
                                pyxlinkviewer_get_annotation(
                                    st.session_state["export_crosslinks_pyxlinkviewer"][
                                        "Parsed PDB sequence"
                                    ],
                                    st.session_state["export_crosslinks_pyxlinkviewer"][
                                        "Parsed PDB chains"
                                    ],
                                    st.session_state["export_crosslinks_pyxlinkviewer"][
                                        "Parsed PDB residue numbers"
                                    ],
                                ),
                                sep=",",
                                index=False,
                            ),
                            file_name="crosslinks_pyxlinkviewer_pdb_annotation.csv",
                            on_click="ignore",
                            type="primary",
                            mime="text/csv",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the parsed PDB annotation.",
                            key="export_crosslinks_pyxlinkviewer_download_meta_pdb_annotation",
                        )
                    export_crosslinks_pyxlinkviewer_download_goto_tool = st.link_button(
                        "Go to PyXlinkViewer!",
                        url="https://github.com/BobSchiffrin/PyXlinkViewer",
                        help="Go to the PyXlinkViewer download page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # xiNET
            elif export_crosslinks_picker == "xiNET":
                export_crosslinks_xinet_info = st.info(
                    "To export to xiNET your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks!"
                    + "It is also required that all crosslinks have "
                    + "associated proteins for the alpha and beta peptide as well as the corresponding crosslink "
                    + "positions in those proteins! It is **not necessary to check this preemptively** as the exporter "
                    + "automatically checks that this information is available and will throw an error otherwise! "
                    + "You can additionally check this yourself in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                export_crosslinks_xinet_button = st.button(
                    "Export to xiNET!",
                    type="primary",
                    width="stretch",
                    key="export_crosslinks_xinet_button",
                )
                if export_crosslinks_xinet_button:
                    with st.spinner(
                        "Exporting crosslinks to xiNET...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_crosslinks_xinet"] = (
                                exporter.to_xinet(crosslinks, filename=None)
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_crosslinks_xinet" in st.session_state
                    and st.session_state["export_crosslinks_xinet"] is not None
                ):
                    export_crosslinks_xinet_download_info = st.markdown(
                        "Your exported crosslinks in xiNET format are ready for download:"
                    )
                    export_crosslinks_xinet_download = st.download_button(
                        label="Download in xiNET format!",
                        data=dataframe_to_csv_stream(
                            st.session_state["export_crosslinks_xinet"],
                            sep=",",
                            index=False,
                        ),
                        file_name="crosslinks_xinet.csv",
                        on_click="ignore",
                        type="primary",
                        mime="text/csv",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in xiNET format.",
                        key="export_crosslinks_xinet_download",
                    )
                    export_crosslinks_xinet_download_goto_tool = st.link_button(
                        "Go to xiNET!",
                        url="https://crosslinkviewer.org/",
                        help="Go to the xiNET website.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # xiVIEW
            elif export_crosslinks_picker == "xiVIEW":
                export_crosslinks_xiview_info = st.info(
                    "To export to xiVIEW your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks!"
                    + "It is also required that all crosslinks have "
                    + "associated proteins for the alpha and beta peptide as well as the corresponding crosslink "
                    + "positions in those proteins! It is **not necessary to check this preemptively** as the exporter "
                    + "automatically checks that this information is available and will throw an error otherwise! "
                    + "You can additionally check this yourself in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                export_crosslinks_xiview_button = st.button(
                    "Export to xiVIEW!",
                    type="primary",
                    width="stretch",
                    key="export_crosslinks_xiview_button",
                )
                if export_crosslinks_xiview_button:
                    with st.spinner(
                        "Exporting crosslinks to xiVIEW...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_crosslinks_xiview"] = (
                                exporter.to_xiview(crosslinks, filename=None)
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_crosslinks_xiview" in st.session_state
                    and st.session_state["export_crosslinks_xiview"] is not None
                ):
                    export_crosslinks_xiview_download_info = st.markdown(
                        "Your exported crosslinks in xiVIEW format are ready for download:"
                    )
                    export_crosslinks_xiview_download = st.download_button(
                        label="Download in xiVIEW format!",
                        data=dataframe_to_csv_stream(
                            st.session_state["export_crosslinks_xiview"],
                            sep=",",
                            index=False,
                        ),
                        file_name="crosslinks_xiview.csv",
                        on_click="ignore",
                        type="primary",
                        mime="text/csv",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in xiVIEW format.",
                        key="export_crosslinks_xiview_download",
                    )
                    export_crosslinks_xiview_download_goto_tool = st.link_button(
                        "Go to xiVIEW!",
                        url="https://xiview.org/",
                        help="Go to the xiVIEW website.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # XlinkDB
            elif export_crosslinks_picker == "XlinkDB":
                export_crosslinks_xlinkdb_info = st.info(
                    "To export to XlinkDB your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks!"
                    + "It is also required that all crosslinks have "
                    + "associated proteins for the alpha and beta peptide! "
                    + "It is **not necessary to check this preemptively** as the exporter "
                    + "automatically checks that this information is available and will throw an error otherwise! "
                    + "You can additionally check this yourself in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                export_crosslinks_xlinkdb_button = st.button(
                    "Export to XlinkDB!",
                    type="primary",
                    width="stretch",
                    key="export_crosslinks_xlinkdb_button",
                )
                if export_crosslinks_xlinkdb_button:
                    with st.spinner(
                        "Exporting crosslinks to XlinkDB...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_crosslinks_xlinkdb"] = (
                                exporter.to_xlinkdb(crosslinks, filename=None)
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_crosslinks_xlinkdb" in st.session_state
                    and st.session_state["export_crosslinks_xlinkdb"] is not None
                ):
                    export_crosslinks_xlinkdb_download_info = st.markdown(
                        "Your exported crosslinks in XlinkDB format are ready for download:"
                    )
                    export_crosslinks_xlinkdb_download = st.download_button(
                        label="Download in XlinkDB format!",
                        data=dataframe_to_csv_stream(
                            st.session_state["export_crosslinks_xlinkdb"],
                            sep="\t",
                            index=False,
                            header=False,
                        ),
                        file_name="crosslinksForXlinkDB.tsv",
                        on_click="ignore",
                        type="primary",
                        mime="text/csv",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in XlinkDB format.",
                        key="export_crosslinks_xlinkdb_download",
                    )
                    export_crosslinks_xlinkdb_download_goto_tool = st.link_button(
                        "Go to XlinkDB!",
                        url="https://xlinkdb.gs.washington.edu/xlinkdb/index.php",
                        help="Go to the XlinkDB website.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # xlms-tools
            elif export_crosslinks_picker == "xlms-tools":
                export_crosslinks_xlmstools_info = st.info(
                    "To export to xlms-tools your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                    + "You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                with st.form(
                    "export_crosslinks_xlmstools_form",
                    enter_to_submit=False,
                    border=False,
                ):
                    xlmstools_crosslinks_pdb_code = st.text_input(
                        "Specify the PDB identification code of your protein(-complex) of interest:",
                        value=None,
                        max_chars=4,
                        key="xlmstools_crosslinks_pdb_code",
                        help="Specify a 4-letter PDB identification code of your cross-linked protein(-complex) of interest.",
                    )
                    xlmstools_crosslinks_pdb_file = st.file_uploader(
                        "Alternatively, upload a PDB file of your protein(-complex) of interest:",
                        type="pdb",
                        accept_multiple_files=False,
                        key="xlmstools_crosslinks_pdb_file",
                        help="Upload a PDB file of your cross-linked protein(-complex) of interest.",
                    )
                    export_crosslinks_xlmstools_button = st.form_submit_button(
                        "Export to xlms-tools format!",
                        type="primary",
                        width="stretch",
                        key="export_crosslinks_xlmstools_button",
                    )
                    if export_crosslinks_xlmstools_button:
                        if (
                            xlmstools_crosslinks_pdb_code is None
                            and xlmstools_crosslinks_pdb_file is None
                        ):
                            _ = st.error(
                                "Can't export to xlms-tools when neither PDB code nor file are given!",
                                icon="⚠️",
                            )
                        else:
                            with st.spinner(
                                "Exporting crosslinks to xlms-tools format...",
                                show_time=True,
                            ):
                                try:
                                    if xlmstools_crosslinks_pdb_file is not None:
                                        st.session_state[
                                            "export_crosslinks_xlmstools"
                                        ] = export_xlmstools_using_pdbfile(
                                            crosslinks, xlmstools_crosslinks_pdb_file
                                        )
                                    else:
                                        if xlmstools_crosslinks_pdb_code is not None:
                                            if (
                                                len(
                                                    xlmstools_crosslinks_pdb_code.strip()
                                                )
                                                != 4
                                            ):
                                                raise ValueError(
                                                    "Specified PDB code is not a valid 4-letter PDB identification code!"
                                                )
                                            st.session_state[
                                                "export_crosslinks_xlmstools"
                                            ] = exporter.to_xlmstools(
                                                crosslinks,
                                                xlmstools_crosslinks_pdb_code.strip(),
                                                filename_prefix=None,
                                            )
                                        else:
                                            raise RuntimeError(
                                                "Can't export to xlms-tools when neither PDB code nor file are given!"
                                            )
                                except Exception as e:
                                    _ = st.error(
                                        "Something went wrong! This is most likely due to missing information in the results!",
                                        icon="⚠️",
                                    )
                                    with st.expander("Show exception"):
                                        _ = st.exception(e)
                if (
                    "export_crosslinks_xlmstools" in st.session_state
                    and st.session_state["export_crosslinks_xlmstools"] is not None
                ):
                    export_crosslinks_xlmstools_download_info = st.markdown(
                        "Your exported crosslinks in xlms-tools format are ready for download:"
                    )
                    export_crosslinks_xlmstools_download = st.download_button(
                        label="Download in xlms-tools format!",
                        data=to_text(
                            st.session_state["export_crosslinks_xlmstools"][
                                "xlms-tools"
                            ]
                        ),
                        file_name="crosslinks_xlmstools.txt",
                        on_click="ignore",
                        type="primary",
                        mime="text/plain",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in xlms-tools format.",
                        key="export_crosslinks_xlmstools_download",
                    )
                    with st.expander("Download Meta-data"):
                        export_crosslinks_xlmstools_download_meta_nr_xl = st.markdown(
                            "**Number of mapped crosslinks:** "
                            + f"{st.session_state['export_crosslinks_xlmstools']['Number of mapped crosslinks']}"
                        )
                        export_crosslinks_xlmstools_download_meta_mapping = st.download_button(
                            label="Download crosslink mapping!",
                            data=to_text(
                                st.session_state["export_crosslinks_xlmstools"][
                                    "Mapping"
                                ]
                            ),
                            file_name="crosslinks_xlmstools_mapping.txt",
                            on_click="ignore",
                            type="primary",
                            mime="text/plain",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the mapping of crosslinks to the PDB structure.",
                            key="export_crosslinks_xlmstools_download_meta_mapping",
                        )
                        export_crosslinks_xlmstools_download_meta_pdb_sequence = st.download_button(
                            label="Download parsed PDB sequence!",
                            data=to_text(
                                pyxlinkviewer_get_fasta(
                                    st.session_state["export_crosslinks_xlmstools"][
                                        "Parsed PDB sequence"
                                    ]
                                )
                            ),
                            file_name="crosslinks_xlmstools_pdb_sequence.fasta",
                            on_click="ignore",
                            type="primary",
                            mime="chemical/seq-aa-fasta",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the parsed PDB sequence.",
                            key="export_crosslinks_xlmstools_download_meta_pdb_sequence",
                        )
                        export_crosslinks_xlmstools_download_meta_pdb_annotation = st.download_button(
                            label="Download parsed PDB annotation!",
                            data=dataframe_to_csv_stream(
                                pyxlinkviewer_get_annotation(
                                    st.session_state["export_crosslinks_xlmstools"][
                                        "Parsed PDB sequence"
                                    ],
                                    st.session_state["export_crosslinks_xlmstools"][
                                        "Parsed PDB chains"
                                    ],
                                    st.session_state["export_crosslinks_xlmstools"][
                                        "Parsed PDB residue numbers"
                                    ],
                                ),
                                sep=",",
                                index=False,
                            ),
                            file_name="crosslinks_xlmstools_pdb_annotation.csv",
                            on_click="ignore",
                            type="primary",
                            mime="text/csv",
                            icon=":material/download:",
                            width="stretch",
                            help="Downloads the parsed PDB annotation.",
                            key="export_crosslinks_xlmstools_download_meta_pdb_annotation",
                        )
                    export_crosslinks_xlmstools_download_goto_tool = st.link_button(
                        "Go to xlms-tools!",
                        url="https://gitlab.com/topf-lab/xlms-tools",
                        help="Go to the xlms-tools project page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            # XMAS
            elif export_crosslinks_picker == "XMAS":
                export_crosslinks_xmas_info = st.info(
                    "To export to XMAS your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                    + "You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                export_crosslinks_xmas_button = st.button(
                    "Export to XMAS format!",
                    type="primary",
                    width="stretch",
                    key="export_crosslinks_xmas_button",
                )
                if export_crosslinks_xmas_button:
                    with st.spinner(
                        "Exporting crosslinks to XMAS format...",
                        show_time=True,
                    ):
                        try:
                            st.session_state["export_crosslinks_xmas"] = (
                                exporter.to_xmas(crosslinks, filename=None)
                            )
                        except Exception as e:
                            _ = st.error(
                                "Something went wrong! This is most likely due to missing information in the results!",
                                icon="⚠️",
                            )
                            with st.expander("Show exception"):
                                _ = st.exception(e)
                if (
                    "export_crosslinks_xmas" in st.session_state
                    and st.session_state["export_crosslinks_xmas"] is not None
                ):
                    export_crosslinks_xmas_download_info = st.markdown(
                        "Your exported crosslinks in XMAS format are ready for download:"
                    )
                    export_crosslinks_xmas_download = st.download_button(
                        label="Download in XMAS format!",
                        data=dataframe_to_xlsx_stream(
                            st.session_state["export_crosslinks_xmas"],
                            sheet_name="xmas",
                            index=False,
                        ),
                        file_name="crosslinks_xmas.xlsx",
                        on_click="ignore",
                        type="primary",
                        mime="application/vnd.ms-excel",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in XMAS format.",
                        key="export_crosslinks_xmas_download",
                    )
                    export_crosslinks_xmas_download_goto_tool = st.link_button(
                        "Go to XMAS!",
                        url="https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle",
                        help="Go to the XMAS project page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
            else:
                pass

    # exporting aggregated crosslinks
    if (
        "aggregated" in st.session_state
        and st.session_state["aggregated"] is not None
        and len(st.session_state["aggregated"]) == 0
    ):
        _ = st.error(
            "Filtering criteria too strict! None of the aggregated crosslinks passed the filter! Please reload your data from the 'Load Data' tab!"
        )
    if (
        "aggregated" in st.session_state
        and st.session_state["aggregated"] is not None
        and len(st.session_state["aggregated"]) > 0
    ):
        aggregated_crosslinks = st.session_state["aggregated"]
        export_aggregated_crosslinks_header = st.subheader(
            "Export Aggregated Crosslinks", divider="grey"
        )
        export_aggregated_crosslinks_options = [
            "AlphaLink2",
            "IMP-X-FDR",
            "MS Annika",
            "PyXlinkViewer",
            "xiNET",
            "xiVIEW",
            "XlinkDB",
            "xlms-tools",
            "XMAS",
        ]
        export_aggregated_crosslinks_picker = st.selectbox(
            "Export aggregated crosslinks to:",
            options=export_aggregated_crosslinks_options,
            index=None,
            help="Choose a format to export the aggregated crosslinks to.",
        )
        if export_aggregated_crosslinks_picker is None:
            pass
        # AlphaLink2
        elif export_aggregated_crosslinks_picker == "AlphaLink2":
            export_aggregated_crosslinks_alphalink2_info = st.info(
                "To export to AlphaLink2 your crosslinks should be **unique** and **not** "
                + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                + "You should also filter your crosslinks to only contain residue pairs of your protein(s) of interest! "
                + "You can check this in the "
                + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
            )
            with st.form(
                "export_aggregated_crosslinks_alphalink2_form",
                enter_to_submit=False,
                border=False,
            ):
                aggregated_crosslinks_alphalink2_fasta_file = st.file_uploader(
                    "Upload a FASTA file of proteins/chains of interest:",
                    type="fasta",
                    accept_multiple_files=False,
                    key="aggregated_crosslinks_alphalink2_fasta_file",
                    help="Upload a FASTA file containing protein/chain sequences. Please keep in mind that AlphaLink2 supports a maximum of 62 proteins/chains!",
                )
                aggregated_crosslinks_alphalink2_annotated_fdr = st.number_input(
                    "Annotated FDR:",
                    value=0.01,
                    min_value=0.0,
                    max_value=1.0,
                    step=0.001,
                    format="%0.3f",
                    key="aggregated_crosslinks_alphalink2_annotated_fdr",
                    help="Value to use for the 'FDR' column in the AlphaLink2 crosslink table, must be given as a real number between 0 and 1. The default of 0.01 corresponds to 1% FDR.",
                )
                export_aggregated_crosslinks_alphalink2_button = st.form_submit_button(
                    "Export to AlphaLink2 format!",
                    type="primary",
                    width="stretch",
                    key="export_aggregated_crosslinks_alphalink2_button",
                )
                if export_aggregated_crosslinks_alphalink2_button:
                    if (
                        aggregated_crosslinks_alphalink2_fasta_file is None
                        or aggregated_crosslinks_alphalink2_annotated_fdr is None
                    ):
                        _ = st.error(
                            "Can't export to AlphaLink2 when either FASTA file or annotated FDR are missing!",
                            icon="⚠️",
                        )
                    else:
                        with st.spinner(
                            "Exporting crosslinks to AlphaLink2 format...",
                            show_time=True,
                        ):
                            try:
                                st.session_state[
                                    "export_aggregated_crosslinks_alphalink2"
                                ] = export_alphalink2(
                                    aggregated_crosslinks,
                                    aggregated_crosslinks_alphalink2_fasta_file,
                                    float(
                                        aggregated_crosslinks_alphalink2_annotated_fdr
                                    ),
                                )
                            except Exception as e:
                                _ = st.error(
                                    "Something went wrong! This is most likely due to missing information in the results!",
                                    icon="⚠️",
                                )
                                with st.expander("Show exception"):
                                    _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_alphalink2" in st.session_state
                and st.session_state["export_aggregated_crosslinks_alphalink2"]
                is not None
            ):
                export_aggregated_crosslinks_alphalink2_download_info = st.markdown(
                    "Your exported crosslinks in AlphaLink2 format are ready for download:"
                )
                (
                    export_aggregated_crosslinks_alphalink2_download_l,
                    export_aggregated_crosslinks_alphalink2_download_m,
                    export_aggregated_crosslinks_alphalink2_download_r,
                ) = st.columns(3)

                with export_aggregated_crosslinks_alphalink2_download_l:
                    export_aggregated_crosslinks_alphalink2_download_txt = st.download_button(
                        label="Download crosslinks in AlphaLink2 format!",
                        data=to_text(
                            st.session_state["export_aggregated_crosslinks_alphalink2"][
                                "AlphaLink2 crosslinks"
                            ]
                        ),
                        file_name="aggregated_crosslinks_AlphaLink2.txt",
                        on_click="ignore",
                        type="primary",
                        mime="text/plain",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in AlphaLink2 format.",
                        key="export_aggregated_crosslinks_alphalink2_download_txt",
                    )
                with export_aggregated_crosslinks_alphalink2_download_m:
                    export_aggregated_crosslinks_alphalink2_download_fasta = st.download_button(
                        label="Download FASTA in AlphaLink2 format!",
                        data=to_text(
                            st.session_state["export_aggregated_crosslinks_alphalink2"][
                                "AlphaLink2 FASTA"
                            ]
                        ),
                        file_name="aggregated_crosslinks_AlphaLink2.fasta",
                        on_click="ignore",
                        type="primary",
                        mime="chemical/seq-aa-fasta",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the uploaded FASTA file in AlphaLink2 format.",
                        key="export_aggregated_crosslinks_alphalink2_download_fasta",
                    )
                with export_aggregated_crosslinks_alphalink2_download_r:
                    export_aggregated_crosslinks_alphalink2_download_pickle = st.download_button(
                        label="Download pickled crosslinks in AlphaLink2 format!",
                        data=pickle_and_gzip(
                            st.session_state["export_aggregated_crosslinks_alphalink2"][
                                "AlphaLink2 Pickle"
                            ]
                        ),
                        file_name="aggregated_crosslinks_AlphaLink2.pickle.gz",
                        on_click="ignore",
                        type="primary",
                        mime="application/gzip",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported crosslinks in pickled AlphaLink2 format.",
                        key="export_aggregated_crosslinks_alphalink2_download_pickle",
                    )
                with st.expander("Show Exported Crosslinks"):
                    export_aggregated_crosslinks_alphalink2_df_info = st.markdown(
                        "**Number of mapped residue pairs:** "
                        + f"{st.session_state['export_aggregated_crosslinks_alphalink2']['AlphaLink2 DataFrame'].shape[0]}"
                    )
                    export_aggregated_crosslinks_alphalink2_df = st.dataframe(
                        st.session_state["export_aggregated_crosslinks_alphalink2"][
                            "AlphaLink2 DataFrame"
                        ],
                        width="stretch",
                    )
                with st.expander("Show Exported FASTA"):
                    export_aggregated_crosslinks_alphalink2_nr_proteins = len(
                        [
                            line
                            for line in st.session_state[
                                "export_aggregated_crosslinks_alphalink2"
                            ]["AlphaLink2 FASTA"].split("\n")
                            if line.startswith(">")
                        ]
                    )
                    export_aggregated_crosslinks_alphalink2_display_fasta_info = (
                        st.markdown(
                            "**Number of mapped proteins/chains:** "
                            + f"{export_aggregated_crosslinks_alphalink2_nr_proteins}"
                        )
                    )
                    export_aggregated_crosslinks_alphalink2_display_fasta = st.text(
                        st.session_state["export_aggregated_crosslinks_alphalink2"][
                            "AlphaLink2 FASTA"
                        ],
                        width="stretch",
                    )

                export_aggregated_crosslinks_alphalink2_download_goto_tool = (
                    st.link_button(
                        "Go to AlphaLink2!",
                        url="https://github.com/Rappsilber-Laboratory/AlphaLink2",
                        help="Go to the AlphaLink2 page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
                )
        # IMP-X-FDR
        elif export_aggregated_crosslinks_picker == "IMP-X-FDR":
            export_aggregated_crosslinks_impxfdr_info = st.info(
                "To export to IMP-X-FDR your aggregated crosslinks should be **unique** and **not** "
                + "**contain any decoy matches**! Usually you would also want to filter for high-confidence matches "
                + "to compare FDR estimation to the experimentally validated FDR! You can check this in the "
                + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
            )
            export_aggregated_crosslinks_impxfdr_button = st.button(
                "Export to IMP-X-FDR!",
                type="primary",
                width="stretch",
                key="export_aggregated_crosslinks_impxfdr_button",
            )
            if export_aggregated_crosslinks_impxfdr_button:
                with st.spinner(
                    "Exporting aggregated crosslinks to IMP-X-FDR...",
                    show_time=True,
                ):
                    try:
                        st.session_state["export_aggregated_crosslinks_impxfdr"] = (
                            exporter.to_impxfdr(aggregated_crosslinks, filename=None)
                        )
                    except Exception as e:
                        _ = st.error(
                            "Something went wrong! This is most likely due to missing information in the results!",
                            icon="⚠️",
                        )
                        with st.expander("Show exception"):
                            _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_impxfdr" in st.session_state
                and st.session_state["export_aggregated_crosslinks_impxfdr"] is not None
            ):
                export_aggregated_crosslinks_impxfdr_download_info = st.markdown(
                    "Your exported aggregated crosslinks in IMP-X-FDR format are ready for download:"
                )
                export_aggregated_crosslinks_impxfdr_download = st.download_button(
                    label="Download in IMP-X-FDR format!",
                    data=dataframe_to_xlsx_stream(
                        st.session_state["export_aggregated_crosslinks_impxfdr"],
                        sheet_name="imp-x-fdr",
                        index=False,
                    ),
                    file_name="aggregated_crosslinks_imp-x-fdr.xlsx",
                    on_click="ignore",
                    type="primary",
                    mime="application/vnd.ms-excel",
                    icon=":material/download:",
                    width="stretch",
                    help="Downloads the exported aggregated crosslinks in IMP-X-FDR format.",
                    key="export_aggregated_crosslinks_impxfdr_download",
                )
                export_aggregated_crosslinks_impxfdr_download_goto_tool = (
                    st.link_button(
                        "Go to IMP-X-FDR!",
                        url="https://github.com/vbc-proteomics-org/imp-x-fdr",
                        help="Go to the IMP-X-FDR download page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
                )
        # MS Annika
        elif export_aggregated_crosslinks_picker == "MS Annika":
            export_aggregated_crosslinks_msannika_button = st.button(
                "Export to MS Annika format!",
                type="primary",
                width="stretch",
                key="export_aggregated_crosslinks_msannika_button",
            )
            if export_aggregated_crosslinks_msannika_button:
                with st.spinner(
                    "Exporting aggregated crosslinks to MS Annika format...",
                    show_time=True,
                ):
                    try:
                        st.session_state["export_aggregated_crosslinks_msannika"] = (
                            exporter.to_msannika(aggregated_crosslinks, filename=None)
                        )
                    except Exception as e:
                        _ = st.error(
                            "Something went wrong! This is most likely due to missing information in the results!",
                            icon="⚠️",
                        )
                        with st.expander("Show exception"):
                            _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_msannika" in st.session_state
                and st.session_state["export_aggregated_crosslinks_msannika"]
                is not None
            ):
                export_aggregated_crosslinks_msannika_download_info = st.markdown(
                    "Your exported aggregated crosslinks in MS Annika format are ready for download:"
                )
                l_aggregated_crosslinks_msannika, r_aggregated_crosslinks_msannika = (
                    st.columns(2)
                )
                with l_aggregated_crosslinks_msannika:
                    export_aggregated_crosslinks_msannika_download_csv = st.download_button(
                        label="Download in MS Annika .csv format!",
                        data=dataframe_to_csv_stream(
                            st.session_state["export_aggregated_crosslinks_msannika"],
                            sep=",",
                            index=False,
                        ),
                        file_name="aggregated_crosslinks_ms-annika.csv",
                        on_click="ignore",
                        type="primary",
                        mime="text/csv",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported aggregated crosslinks in MS Annika comma-separated-values (.csv) format.",
                        key="export_aggregated_crosslinks_msannika_download_csv",
                    )
                with r_aggregated_crosslinks_msannika:
                    export_aggregated_crosslinks_msannika_download_xlsx = st.download_button(
                        label="Download in MS Annika .xlsx format!",
                        data=dataframe_to_xlsx_stream(
                            st.session_state["export_aggregated_crosslinks_msannika"],
                            sheet_name="msannika",
                            index=False,
                        ),
                        file_name="aggregated_crosslinks_ms-annika.xlsx",
                        on_click="ignore",
                        type="primary",
                        mime="application/vnd.ms-excel",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the exported aggregated crosslinks in MS Annika Microsoft Excel (.xlsx) format.",
                        key="export_aggregated_crosslinks_msannika_download_xlsx",
                    )
        # PyXlinkViewer
        elif export_aggregated_crosslinks_picker == "PyXlinkViewer":
            export_aggregated_crosslinks_pyxlinkviewer_info = st.info(
                "To export to PyXlinkViewer your aggregated crosslinks should be **unique** and **not** "
                + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                + "You can check this in the "
                + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
            )
            with st.form(
                "export_aggregated_crosslinks_pyxlinkviewer_form",
                enter_to_submit=False,
                border=False,
            ):
                aggregated_crosslinks_pdb_code = st.text_input(
                    "Specify the PDB identification code of your protein(-complex) of interest:",
                    value=None,
                    max_chars=4,
                    key="aggregated_crosslinks_pdb_code",
                    help="Specify a 4-letter PDB identification code of your cross-linked protein(-complex) of interest.",
                )
                aggregated_crosslinks_pdb_file = st.file_uploader(
                    "Alternatively, upload a PDB file of your protein(-complex) of interest:",
                    type="pdb",
                    accept_multiple_files=False,
                    key="aggregated_crosslinks_pdb_file",
                    help="Upload a PDB file of your cross-linked protein(-complex) of interest.",
                )
                export_aggregated_crosslinks_pyxlinkviewer_button = (
                    st.form_submit_button(
                        "Export to PyXlinkViewer format!",
                        type="primary",
                        width="stretch",
                        key="export_aggregated_crosslinks_pyxlinkviewer_button",
                    )
                )
                if export_aggregated_crosslinks_pyxlinkviewer_button:
                    if (
                        aggregated_crosslinks_pdb_code is None
                        and aggregated_crosslinks_pdb_file is None
                    ):
                        _ = st.error(
                            "Can't export to PyXlinkViewer when neither PDB code nor file are given!",
                            icon="⚠️",
                        )
                    else:
                        with st.spinner(
                            "Exporting aggregated crosslinks to PyXlinkViewer format...",
                            show_time=True,
                        ):
                            try:
                                if aggregated_crosslinks_pdb_file is not None:
                                    st.session_state[
                                        "export_aggregated_crosslinks_pyxlinkviewer"
                                    ] = export_pyxlinkviewer_using_pdbfile(
                                        aggregated_crosslinks,
                                        aggregated_crosslinks_pdb_file,
                                    )
                                else:
                                    if aggregated_crosslinks_pdb_code is not None:
                                        if (
                                            len(aggregated_crosslinks_pdb_code.strip())
                                            != 4
                                        ):
                                            raise ValueError(
                                                "Specified PDB code is not a valid 4-letter PDB identification code!"
                                            )
                                        st.session_state[
                                            "export_aggregated_crosslinks_pyxlinkviewer"
                                        ] = exporter.to_pyxlinkviewer(
                                            aggregated_crosslinks,
                                            aggregated_crosslinks_pdb_code.strip(),
                                            filename_prefix=None,
                                        )
                                    else:
                                        raise RuntimeError(
                                            "Can't export to PyXlinkViewer when neither PDB code nor file are given!"
                                        )
                            except Exception as e:
                                _ = st.error(
                                    "Something went wrong! This is most likely due to missing information in the results!",
                                    icon="⚠️",
                                )
                                with st.expander("Show exception"):
                                    _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_pyxlinkviewer" in st.session_state
                and st.session_state["export_aggregated_crosslinks_pyxlinkviewer"]
                is not None
            ):
                export_aggregated_crosslinks_pyxlinkviewer_download_info = st.markdown(
                    "Your exported aggregated crosslinks in PyXlinkViewer format are ready for download:"
                )
                export_aggregated_crosslinks_pyxlinkviewer_download = st.download_button(
                    label="Download in PyXlinkViewer format!",
                    data=to_text(
                        st.session_state["export_aggregated_crosslinks_pyxlinkviewer"][
                            "PyXlinkViewer"
                        ]
                    ),
                    file_name="aggregated_crosslinks_pyxlinkviewer.txt",
                    on_click="ignore",
                    type="primary",
                    mime="text/plain",
                    icon=":material/download:",
                    width="stretch",
                    help="Downloads the exported aggregated crosslinks in PyXlinkViewer format.",
                    key="export_aggregated_crosslinks_pyxlinkviewer_download",
                )
                with st.expander("Download Meta-data"):
                    export_aggregated_crosslinks_pyxlinkviewer_download_meta_nr_xl = st.markdown(
                        "**Number of mapped crosslinks:** "
                        + f"{st.session_state['export_aggregated_crosslinks_pyxlinkviewer']['Number of mapped crosslinks']}"
                    )
                    export_aggregated_crosslinks_pyxlinkviewer_download_meta_mapping = st.download_button(
                        label="Download crosslink mapping!",
                        data=to_text(
                            st.session_state[
                                "export_aggregated_crosslinks_pyxlinkviewer"
                            ]["Mapping"]
                        ),
                        file_name="aggregated_crosslinks_pyxlinkviewer_mapping.txt",
                        on_click="ignore",
                        type="primary",
                        mime="text/plain",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the mapping of aggregated crosslinks to the PDB structure.",
                        key="export_aggregated_crosslinks_pyxlinkviewer_download_meta_mapping",
                    )
                    export_aggregated_crosslinks_pyxlinkviewer_download_meta_pdb_sequence = st.download_button(
                        label="Download parsed PDB sequence!",
                        data=to_text(
                            pyxlinkviewer_get_fasta(
                                st.session_state[
                                    "export_aggregated_crosslinks_pyxlinkviewer"
                                ]["Parsed PDB sequence"]
                            )
                        ),
                        file_name="aggregated_crosslinks_pyxlinkviewer_pdb_sequence.fasta",
                        on_click="ignore",
                        type="primary",
                        mime="chemical/seq-aa-fasta",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the parsed PDB sequence.",
                        key="export_aggregated_crosslinks_pyxlinkviewer_download_meta_pdb_sequence",
                    )
                    export_aggregated_crosslinks_pyxlinkviewer_download_meta_pdb_annotation = st.download_button(
                        label="Download parsed PDB annotation!",
                        data=dataframe_to_csv_stream(
                            pyxlinkviewer_get_annotation(
                                st.session_state[
                                    "export_aggregated_crosslinks_pyxlinkviewer"
                                ]["Parsed PDB sequence"],
                                st.session_state[
                                    "export_aggregated_crosslinks_pyxlinkviewer"
                                ]["Parsed PDB chains"],
                                st.session_state[
                                    "export_aggregated_crosslinks_pyxlinkviewer"
                                ]["Parsed PDB residue numbers"],
                            ),
                            sep=",",
                            index=False,
                        ),
                        file_name="aggregated_crosslinks_pyxlinkviewer_pdb_annotation.csv",
                        on_click="ignore",
                        type="primary",
                        mime="text/csv",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the parsed PDB annotation.",
                        key="export_aggregated_crosslinks_pyxlinkviewer_download_meta_pdb_annotation",
                    )
                export_aggregated_crosslinks_pyxlinkviewer_download_goto_tool = (
                    st.link_button(
                        "Go to PyXlinkViewer!",
                        url="https://github.com/BobSchiffrin/PyXlinkViewer",
                        help="Go to the PyXlinkViewer download page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
                )
        # xiNET
        elif export_aggregated_crosslinks_picker == "xiNET":
            export_aggregated_crosslinks_xinet_info = st.info(
                "To export to xiNET your aggregated crosslinks should be **unique** and **not** "
                + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks!"
                + "It is also required that all crosslinks have "
                + "associated proteins for the alpha and beta peptide as well as the corresponding crosslink "
                + "positions in those proteins! It is **not necessary to check this preemptively** as the exporter "
                + "automatically checks that this information is available and will throw an error otherwise! "
                + "You can additionally check this yourself in the "
                + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
            )
            export_aggregated_crosslinks_xinet_button = st.button(
                "Export to xiNET!",
                type="primary",
                width="stretch",
                key="export_aggregated_crosslinks_xinet_button",
            )
            if export_aggregated_crosslinks_xinet_button:
                with st.spinner(
                    "Exporting aggregated crosslinks to xiNET...",
                    show_time=True,
                ):
                    try:
                        st.session_state["export_aggregated_crosslinks_xinet"] = (
                            exporter.to_xinet(aggregated_crosslinks, filename=None)
                        )
                    except Exception as e:
                        _ = st.error(
                            "Something went wrong! This is most likely due to missing information in the results!",
                            icon="⚠️",
                        )
                        with st.expander("Show exception"):
                            _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_xinet" in st.session_state
                and st.session_state["export_aggregated_crosslinks_xinet"] is not None
            ):
                export_aggregated_crosslinks_xinet_download_info = st.markdown(
                    "Your exported aggregated crosslinks in xiNET format are ready for download:"
                )
                export_aggregated_crosslinks_xinet_download = st.download_button(
                    label="Download in xiNET format!",
                    data=dataframe_to_csv_stream(
                        st.session_state["export_aggregated_crosslinks_xinet"],
                        sep=",",
                        index=False,
                    ),
                    file_name="aggregated_crosslinks_xinet.csv",
                    on_click="ignore",
                    type="primary",
                    mime="text/csv",
                    icon=":material/download:",
                    width="stretch",
                    help="Downloads the exported aggregated crosslinks in xiNET format.",
                    key="export_aggregated_crosslinks_xinet_download",
                )
                export_aggregated_crosslinks_xinet_download_goto_tool = st.link_button(
                    "Go to xiNET!",
                    url="https://crosslinkviewer.org/",
                    help="Go to the xiNET website.",
                    type="primary",
                    icon="🔗",
                    width="stretch",
                )
        # xiVIEW
        elif export_aggregated_crosslinks_picker == "xiVIEW":
            export_aggregated_crosslinks_xiview_info = st.info(
                "To export to xiVIEW your aggregated crosslinks should be **unique** and **not** "
                + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks!"
                + "It is also required that all crosslinks have "
                + "associated proteins for the alpha and beta peptide as well as the corresponding crosslink "
                + "positions in those proteins! It is **not necessary to check this preemptively** as the exporter "
                + "automatically checks that this information is available and will throw an error otherwise! "
                + "You can additionally check this yourself in the "
                + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
            )
            export_aggregated_crosslinks_xiview_button = st.button(
                "Export to xiVIEW!",
                type="primary",
                width="stretch",
                key="export_aggregated_crosslinks_xiview_button",
            )
            if export_aggregated_crosslinks_xiview_button:
                with st.spinner(
                    "Exporting aggregated crosslinks to xiVIEW...",
                    show_time=True,
                ):
                    try:
                        st.session_state["export_aggregated_crosslinks_xiview"] = (
                            exporter.to_xiview(aggregated_crosslinks, filename=None)
                        )
                    except Exception as e:
                        _ = st.error(
                            "Something went wrong! This is most likely due to missing information in the results!",
                            icon="⚠️",
                        )
                        with st.expander("Show exception"):
                            _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_xiview" in st.session_state
                and st.session_state["export_aggregated_crosslinks_xiview"] is not None
            ):
                export_aggregated_crosslinks_xiview_download_info = st.markdown(
                    "Your exported aggregated crosslinks in xiVIEW format are ready for download:"
                )
                export_aggregated_crosslinks_xiview_download = st.download_button(
                    label="Download in xiVIEW format!",
                    data=dataframe_to_csv_stream(
                        st.session_state["export_aggregated_crosslinks_xiview"],
                        sep=",",
                        index=False,
                    ),
                    file_name="aggregated_crosslinks_xiview.csv",
                    on_click="ignore",
                    type="primary",
                    mime="text/csv",
                    icon=":material/download:",
                    width="stretch",
                    help="Downloads the exported aggregated crosslinks in xiVIEW format.",
                    key="export_aggregated_crosslinks_xiview_download",
                )
                export_aggregated_crosslinks_xiview_download_goto_tool = st.link_button(
                    "Go to xiVIEW!",
                    url="https://xiview.org/",
                    help="Go to the xiVIEW website.",
                    type="primary",
                    icon="🔗",
                    width="stretch",
                )
        # XlinkDB
        elif export_aggregated_crosslinks_picker == "XlinkDB":
            export_aggregated_crosslinks_xlinkdb_info = st.info(
                "To export to XlinkDB your aggregated crosslinks should be **unique** and **not** "
                + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks!"
                + "It is also required that all crosslinks have "
                + "associated proteins for the alpha and beta peptide! "
                + "It is **not necessary to check this preemptively** as the exporter "
                + "automatically checks that this information is available and will throw an error otherwise! "
                + "You can additionally check this yourself in the "
                + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
            )
            export_aggregated_crosslinks_xlinkdb_button = st.button(
                "Export to XlinkDB!",
                type="primary",
                width="stretch",
                key="export_aggregated_crosslinks_xlinkdb_button",
            )
            if export_aggregated_crosslinks_xlinkdb_button:
                with st.spinner(
                    "Exporting aggregated crosslinks to XlinkDB...",
                    show_time=True,
                ):
                    try:
                        st.session_state["export_aggregated_crosslinks_xlinkdb"] = (
                            exporter.to_xlinkdb(aggregated_crosslinks, filename=None)
                        )
                    except Exception as e:
                        _ = st.error(
                            "Something went wrong! This is most likely due to missing information in the results!",
                            icon="⚠️",
                        )
                        with st.expander("Show exception"):
                            _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_xlinkdb" in st.session_state
                and st.session_state["export_aggregated_crosslinks_xlinkdb"] is not None
            ):
                export_aggregated_crosslinks_xlinkdb_download_info = st.markdown(
                    "Your exported aggregated crosslinks in XlinkDB format are ready for download:"
                )
                export_aggregated_crosslinks_xlinkdb_download = st.download_button(
                    label="Download in XlinkDB format!",
                    data=dataframe_to_csv_stream(
                        st.session_state["export_aggregated_crosslinks_xlinkdb"],
                        sep="\t",
                        index=False,
                        header=False,
                    ),
                    file_name="aggregatedCrosslinksForXlinkDB.tsv",
                    on_click="ignore",
                    type="primary",
                    mime="text/csv",
                    icon=":material/download:",
                    width="stretch",
                    help="Downloads the exported aggregated crosslinks in XlinkDB format.",
                    key="export_aggregated_crosslinks_xlinkdb_download",
                )
                export_aggregated_crosslinks_xlinkdb_download_goto_tool = (
                    st.link_button(
                        "Go to XlinkDB!",
                        url="https://xlinkdb.gs.washington.edu/xlinkdb/index.php",
                        help="Go to the XlinkDB website.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
                )
        # xlms-tools
        elif export_aggregated_crosslinks_picker == "xlms-tools":
            export_aggregated_crosslinks_xlmstools_info = st.info(
                "To export to xlms-tools your aggregated crosslinks should be **unique** and **not** "
                + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                + "You can check this in the "
                + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
            )
            with st.form(
                "export_aggregated_crosslinks_xlmstools_form",
                enter_to_submit=False,
                border=False,
            ):
                xlmstools_aggregated_crosslinks_pdb_code = st.text_input(
                    "Specify the PDB identification code of your protein(-complex) of interest:",
                    value=None,
                    max_chars=4,
                    key="xlmstools_aggregated_crosslinks_pdb_code",
                    help="Specify a 4-letter PDB identification code of your cross-linked protein(-complex) of interest.",
                )
                xlmstools_aggregated_crosslinks_pdb_file = st.file_uploader(
                    "Alternatively, upload a PDB file of your protein(-complex) of interest:",
                    type="pdb",
                    accept_multiple_files=False,
                    key="xlmstools_aggregated_crosslinks_pdb_file",
                    help="Upload a PDB file of your cross-linked protein(-complex) of interest.",
                )
                export_aggregated_crosslinks_xlmstools_button = st.form_submit_button(
                    "Export to xlms-tools format!",
                    type="primary",
                    width="stretch",
                    key="export_aggregated_crosslinks_xlmstools_button",
                )
                if export_aggregated_crosslinks_xlmstools_button:
                    if (
                        xlmstools_aggregated_crosslinks_pdb_code is None
                        and xlmstools_aggregated_crosslinks_pdb_file is None
                    ):
                        _ = st.error(
                            "Can't export to xlms-tools when neither PDB code nor file are given!",
                            icon="⚠️",
                        )
                    else:
                        with st.spinner(
                            "Exporting aggregated crosslinks to xlms-tools format...",
                            show_time=True,
                        ):
                            try:
                                if xlmstools_aggregated_crosslinks_pdb_file is not None:
                                    st.session_state[
                                        "export_aggregated_crosslinks_xlmstools"
                                    ] = export_xlmstools_using_pdbfile(
                                        aggregated_crosslinks,
                                        xlmstools_aggregated_crosslinks_pdb_file,
                                    )
                                else:
                                    if (
                                        xlmstools_aggregated_crosslinks_pdb_code
                                        is not None
                                    ):
                                        if (
                                            len(
                                                xlmstools_aggregated_crosslinks_pdb_code.strip()
                                            )
                                            != 4
                                        ):
                                            raise ValueError(
                                                "Specified PDB code is not a valid 4-letter PDB identification code!"
                                            )
                                        st.session_state[
                                            "export_aggregated_crosslinks_xlmstools"
                                        ] = exporter.to_xlmstools(
                                            aggregated_crosslinks,
                                            xlmstools_aggregated_crosslinks_pdb_code.strip(),
                                            filename_prefix=None,
                                        )
                                    else:
                                        raise RuntimeError(
                                            "Can't export to xlms-tools when neither PDB code nor file are given!"
                                        )
                            except Exception as e:
                                _ = st.error(
                                    "Something went wrong! This is most likely due to missing information in the results!",
                                    icon="⚠️",
                                )
                                with st.expander("Show exception"):
                                    _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_xlmstools" in st.session_state
                and st.session_state["export_aggregated_crosslinks_xlmstools"]
                is not None
            ):
                export_aggregated_crosslinks_xlmstools_download_info = st.markdown(
                    "Your exported aggregated crosslinks in xlms-tools format are ready for download:"
                )
                export_aggregated_crosslinks_xlmstools_download = st.download_button(
                    label="Download in xlms-tools format!",
                    data=to_text(
                        st.session_state["export_aggregated_crosslinks_xlmstools"][
                            "xlms-tools"
                        ]
                    ),
                    file_name="aggregated_crosslinks_xlmstools.txt",
                    on_click="ignore",
                    type="primary",
                    mime="text/plain",
                    icon=":material/download:",
                    width="stretch",
                    help="Downloads the exported aggregated crosslinks in xlms-tools format.",
                    key="export_aggregated_crosslinks_xlmstools_download",
                )
                with st.expander("Download Meta-data"):
                    export_aggregated_crosslinks_xlmstools_download_meta_nr_xl = st.markdown(
                        "**Number of mapped crosslinks:** "
                        + f"{st.session_state['export_aggregated_crosslinks_xlmstools']['Number of mapped crosslinks']}"
                    )
                    export_aggregated_crosslinks_xlmstools_download_meta_mapping = st.download_button(
                        label="Download crosslink mapping!",
                        data=to_text(
                            st.session_state["export_aggregated_crosslinks_xlmstools"][
                                "Mapping"
                            ]
                        ),
                        file_name="aggregated_crosslinks_xlmstools_mapping.txt",
                        on_click="ignore",
                        type="primary",
                        mime="text/plain",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the mapping of aggregated crosslinks to the PDB structure.",
                        key="export_aggregated_crosslinks_xlmstools_download_meta_mapping",
                    )
                    export_aggregated_crosslinks_xlmstools_download_meta_pdb_sequence = st.download_button(
                        label="Download parsed PDB sequence!",
                        data=to_text(
                            pyxlinkviewer_get_fasta(
                                st.session_state[
                                    "export_aggregated_crosslinks_xlmstools"
                                ]["Parsed PDB sequence"]
                            )
                        ),
                        file_name="aggregated_crosslinks_xlmstools_pdb_sequence.fasta",
                        on_click="ignore",
                        type="primary",
                        mime="chemical/seq-aa-fasta",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the parsed PDB sequence.",
                        key="export_aggregated_crosslinks_xlmstools_download_meta_pdb_sequence",
                    )
                    export_aggregated_crosslinks_xlmstools_download_meta_pdb_annotation = st.download_button(
                        label="Download parsed PDB annotation!",
                        data=dataframe_to_csv_stream(
                            pyxlinkviewer_get_annotation(
                                st.session_state[
                                    "export_aggregated_crosslinks_xlmstools"
                                ]["Parsed PDB sequence"],
                                st.session_state[
                                    "export_aggregated_crosslinks_xlmstools"
                                ]["Parsed PDB chains"],
                                st.session_state[
                                    "export_aggregated_crosslinks_xlmstools"
                                ]["Parsed PDB residue numbers"],
                            ),
                            sep=",",
                            index=False,
                        ),
                        file_name="aggregated_crosslinks_xlmstools_pdb_annotation.csv",
                        on_click="ignore",
                        type="primary",
                        mime="text/csv",
                        icon=":material/download:",
                        width="stretch",
                        help="Downloads the parsed PDB annotation.",
                        key="export_aggregated_crosslinks_xlmstools_download_meta_pdb_annotation",
                    )
                export_aggregated_crosslinks_xlmstools_download_goto_tool = (
                    st.link_button(
                        "Go to xlms-tools!",
                        url="https://gitlab.com/topf-lab/xlms-tools",
                        help="Go to the xlms-tools project page.",
                        type="primary",
                        icon="🔗",
                        width="stretch",
                    )
                )
        # XMAS
        elif export_aggregated_crosslinks_picker == "XMAS":
            export_aggregated_crosslinks_xmas_info = st.info(
                "To export to XMAS your aggregated crosslinks should be **unique** and **not** "
                + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                + "You can check this in the "
                + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
            )
            export_aggregated_crosslinks_xmas_button = st.button(
                "Export to XMAS format!",
                type="primary",
                width="stretch",
                key="export_aggregated_crosslinks_xmas_button",
            )
            if export_aggregated_crosslinks_xmas_button:
                with st.spinner(
                    "Exporting aggregated crosslinks to XMAS format...",
                    show_time=True,
                ):
                    try:
                        st.session_state["export_aggregated_crosslinks_xmas"] = (
                            exporter.to_xmas(aggregated_crosslinks, filename=None)
                        )
                    except Exception as e:
                        _ = st.error(
                            "Something went wrong! This is most likely due to missing information in the results!",
                            icon="⚠️",
                        )
                        with st.expander("Show exception"):
                            _ = st.exception(e)
            if (
                "export_aggregated_crosslinks_xmas" in st.session_state
                and st.session_state["export_aggregated_crosslinks_xmas"] is not None
            ):
                export_aggregated_crosslinks_xmas_download_info = st.markdown(
                    "Your exported aggregated crosslinks in XMAS format are ready for download:"
                )
                export_aggregated_crosslinks_xmas_download = st.download_button(
                    label="Download in XMAS format!",
                    data=dataframe_to_xlsx_stream(
                        st.session_state["export_aggregated_crosslinks_xmas"],
                        sheet_name="xmas",
                        index=False,
                    ),
                    file_name="aggregated_crosslinks_xmas.xlsx",
                    on_click="ignore",
                    type="primary",
                    mime="application/vnd.ms-excel",
                    icon=":material/download:",
                    width="stretch",
                    help="Downloads the exported aggregated crosslinks in XMAS format.",
                    key="export_aggregated_crosslinks_xmas_download",
                )
                export_aggregated_crosslinks_xmas_download_goto_tool = st.link_button(
                    "Go to XMAS!",
                    url="https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle",
                    help="Go to the XMAS project page.",
                    type="primary",
                    icon="🔗",
                    width="stretch",
                )
        else:
            pass


# about tab
def about_tab():
    general_description = """
        **pyXLMS** is a python package and web application with graphical user interface that aims to simplify and streamline the intermediate step of
        connecting crosslink search engine results with down-stream analysis tools, enabling researchers even without bioinformatics knowledge to
        conduct in-depth crosslink analyses and shifting the focus from data transformation to data interpretation and therefore gaining biological
        insight.

        Currently pyXLMS supports input from seven different crosslink search engines:
        [MaxLynx (part of MaxQuant)](https://www.maxquant.org/),
        [MeroX](https://www.stavrox.com/),
        [MS Annika](https://github.com/hgb-bin-proteomics/MSAnnika),
        [pLink 2 and pLink 3](http://pfind.ict.ac.cn/se/plink/),
        [Scout](https://github.com/diogobor/Scout),
        [xiSearch](https://www.rappsilberlab.org/software/xisearch/) and [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
        [XlinkX](https://docs.thermofisher.com/r/XlinkX-3.2-Quick-Start-Guide/),
        as well as the [mzIdentML format](https://www.psidev.info/mzidentml) of the HUPO Proteomics Standards Initiative,
        and a well-documented and [human-readable custom tabular format](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md).

        Down-stream analysis is facilitated by functionality that is directly available within pyXLMS such as validation, annotation, aggregation,
        filtering, and visualization - and [much more](https://hgb-bin-proteomics.github.io/pyXLMS/modules.html) - of crosslink-spectrum-matches and crosslinks.

        In addition, the data can easily be exported to the required data format of the various available down-stream analysis tools such as
        [AlphaLink2](https://github.com/Rappsilber-Laboratory/AlphaLink2),
        [ProXL](https://www.yeastrc.org/proxl_public/),
        [xiNET](https://crosslinkviewer.org/index.php),
        [xiVIEW](https://www.xiview.org/index.php),
        [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
        [XlinkDB](https://xlinkdb.gs.washington.edu/xlinkdb/),
        [xlms-tools](https://gitlab.com/topf-lab/xlms-tools),
        PyMOL (via [PyXlinkViewer](https://github.com/BobSchiffrin/PyXlinkViewer)),
        ChimeraX (via [XMAS](https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle)),
        or [IMP-X-FDR](https://github.com/vbc-proteomics-org/imp-x-fdr).
        """
    description = st.markdown(general_description)

    header_2 = st.subheader("Citing", divider="grey")
    citation_str = """
        If you are using pyXLMS please cite the following publication:

        - Manuscript in preparation
          ```
          (wip)
          ```
        """
    citation = st.markdown(citation_str)

    header_3 = st.subheader("Contact", divider="grey")
    contact_str = """
        - [proteomics@fh-hagenberg.at](mailto:proteomics@fh-hagenberg.at)
        - [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at) (primary developer)
        """
    contact = st.markdown(contact_str)

    header_4 = st.subheader("Further Links", divider="grey")
    further_info = st.markdown("Read more about pyXLMS at the links below:")

    l1, center_1, r1 = st.columns(3)

    with l1:
        link_button_1 = st.link_button(
            "GitHub",
            url="https://github.com/hgb-bin-proteomics/pyXLMS",
            type="primary",
            help="Link to the pyXLMS GitHub page.",
            width="stretch",
        )

    with center_1:
        link_button_2 = st.link_button(
            "User Guide",
            url="https://hgb-bin-proteomics.github.io/pyXLMS-docs",
            type="primary",
            help="Link to the pyXLMS user guide page.",
            width="stretch",
        )

    with r1:
        link_button_3 = st.link_button(
            "Documentation",
            url="https://hgb-bin-proteomics.github.io/pyXLMS",
            type="primary",
            help="Link to the pyXLMS documentation page.",
            width="stretch",
        )


# main page content
def main_page():
    title = st.title("pyXLMS")

    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        ["Load Data", "Filter", "Visualize", "Export", "About"]
    )

    with tab1:
        input_tab()

    with tab2:
        filter_tab()

    with tab3:
        visualize_tab()

    with tab4:
        export_tab()

    with tab5:
        about_tab()


# side bar and main page loader
def main():
    about_str = """
    **pyXLMS** is a python package and web application with graphical user interface that aims to simplify and streamline the intermediate step of
    connecting crosslink search engine results with down-stream analysis tools, enabling researchers even without bioinformatics knowledge to
    conduct in-depth crosslink analyses and shifting the focus from data transformation to data interpretation and therefore gaining biological
    insight.
    """

    st.set_page_config(
        page_title="pyXLMS",
        page_icon=":dna:",
        layout="wide",
        initial_sidebar_state="auto",
        menu_items={
            "Get Help": "https://github.com/hgb-bin-proteomics/pyXLMS/discussions",
            "Report a bug": "https://github.com/hgb-bin-proteomics/pyXLMS/issues",
            "About": about_str,
        },
    )

    title = st.sidebar.title("pyXLMS")

    logo = st.sidebar.image(
        ".streamlit/icon/logo.png",
        caption="A python package to process protein cross-linking data.",
    )

    div_1 = st.sidebar.divider()

    doc = st.sidebar.markdown(about_str)

    div_2 = st.sidebar.divider()

    info_str = ""
    info_str += "- **Documentation:**  \n [hgb-bin-proteomics.github.io/pyXLMS-docs/](https://hgb-bin-proteomics.github.io/pyXLMS-docs/)\n"
    info_str += "- **Contact:**  \n  [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)\n"
    info_str += "- **License:**  \n  [MIT License](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/LICENSE)\n"
    info_str += "- **Project Page:**  \n  [GitHub](https://github.com/hgb-bin-proteomics/pyXLMS/)"
    info = st.sidebar.markdown(info_str)

    div_3 = st.sidebar.divider()

    get_help = st.sidebar.link_button(
        "Help",
        url=HELP_URL,
        type="secondary",
        # help="Link to the pyXLMS webapp documentation page.",
        width="stretch",
    )

    version_info = st.sidebar.markdown(
        f"Server is running web app version {__version__} and pyXLMS version {__pyxlms_version__}."
    )

    main_page()

    if "display_help_toast" not in st.session_state:
        st.session_state["display_help_toast"] = False
        help_msg = st.toast(
            f"Need help? Read the [docs]({HELP_URL})!",
            icon="🤔",
            duration="infinite",
        )


if __name__ == "__main__":
    main()
