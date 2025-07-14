#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "streamlit",
#   "pyxlms>=1.1",
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
import json
import pandas as pd
from tempfile import NamedTemporaryFile

from pyXLMS import parser
from pyXLMS import transform
from pyXLMS import constants
from pyXLMS import plotting
from pyXLMS import exporter

import streamlit as st

from typing import Optional
from typing import Dict
from typing import List
from typing import Any


@st.cache_data
def to_text(data: str) -> bytes:
    return data.encode("utf-8")


@st.cache_data
def to_json(data: Dict[str, Any]) -> bytes:
    return json.dumps(data).encode("utf-8")


@st.cache_data
def dataframe_to_csv_stream(dataframe: pd.DataFrame, sep: str, index: bool) -> bytes:
    return dataframe.to_csv(sep=sep, index=index).encode("utf-8")


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
def read_file(
    uploaded_file: io.BytesIO,
    engine: str,
    crosslinker: str,
    parse_modifications: bool,
    crosslinker_mass: Optional[float],
) -> Dict[str, Any]:
    #
    with NamedTemporaryFile(
        suffix=os.path.splitext(uploaded_file.name)[1], delete_on_close=False
    ) as f:  # pyright: ignore[reportCallIssue]
        f.write(uploaded_file.getbuffer())
        f.close()
        if crosslinker_mass is not None:
            try:
                return parser.read(
                    f.name,
                    engine=engine,  # pyright: ignore[reportArgumentType]
                    crosslinker=crosslinker,
                    parse_modifications=parse_modifications,
                    crosslinker_mass=crosslinker_mass,
                )
            except Exception as _e:
                parser.read(
                    f.name,
                    engine=engine,  # pyright: ignore[reportArgumentType]
                    crosslinker=crosslinker,
                    parse_modifications=parse_modifications,
                )
        return parser.read(
            f.name,
            engine=engine,  # pyright: ignore[reportArgumentType]
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


# input tab
def input_tab():
    general_description = """
    _a python package to process protein cross-linking data_

    **pyXLMS** is a python package and web application with graphical user interface that aims to simplify and streamline the intermediate step of
    connecting crosslink search engine results with down-stream analysis tools, enabling researchers even without bioinformatics knowledge to
    conduct in-depth crosslink analyses and shifting the focus from data transformation to data interpretation and therefore gaining biological
    insight.

    Currently pyXLMS supports input from six different crosslink search engines:
    [MaxLynx (part of MaxQuant)](https://www.maxquant.org/),
    [MS Annika](https://github.com/hgb-bin-proteomics/MSAnnika),
    [pLink 2 and pLink 3](http://pfind.ict.ac.cn/se/plink/),
    [Scout](https://github.com/diogobor/Scout),
    [xiSearch](https://www.rappsilberlab.org/software/xisearch/) and [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
    [XlinkX](https://docs.thermofisher.com/r/XlinkX-3.2-Quick-Start-Guide/),
    as well as the [mzIdentML format](https://www.psidev.info/mzidentml) of the HUPO Proteomics Standards Initiative,
    and a well-documented and [human-readable custom tabular format](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md).

    Down-stream analysis is facilitated by functionality that is directly available within pyXLMS such as validation, annotation, aggregation,
    and filtering of crosslink-spectrum-matches and crosslinks.

    In addition, the data can easily be exported to the required data format of the various available down-stream analysis tools such as
    [xiNET](https://crosslinkviewer.org/index.php),
    [xiVIEW](https://www.xiview.org/index.php),
    [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
    [XlinkDB](https://xlinkdb.gs.washington.edu/xlinkdb/),
    [xlms-tools](https://gitlab.com/topf-lab/xlms-tools),
    pyMOL (via [pyXlinkViewer](https://github.com/BobSchiffrin/PyXlinkViewer)),
    ChimeraX (via [XMAS](https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle)),
    or [IMP-X-FDR](https://github.com/vbc-proteomics-org/imp-x-fdr).

    **Try it yourself below!** 😉
    """
    description = st.markdown(general_description)

    header_1 = st.subheader("File Upload", divider="grey")

    uploaded_file = st.file_uploader(
        "Upload a cross-linking result file from any of the supported search engines or formats:",
        type=None,
        accept_multiple_files=False,
        key="uploaded_file",
        help="Upload a cross-linking result file from any of the supported search engines or formats.",
    )

    l1, r1 = st.columns(2)

    with l1:
        search_engine = st.selectbox(
            "Select a crosslink search engine or file format:",
            options=[
                "Custom",
                "MaxQuant",
                "MaxLynx",
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
        read_file_button = st.button(
            "Read file!", type="primary", use_container_width=True
        )

    if read_file_button:
        # reset any exported files
        st.session_state["export_csms_impxfdr"] = None
        st.session_state["export_csms_msannika"] = None
        st.session_state["export_csms_xifdr"] = None
        st.session_state["export_crosslinks_impxfdr"] = None
        st.session_state["export_crosslinks_msannika"] = None
        st.session_state["export_crosslinks_pyxlinkviewer"] = None
        st.session_state["export_crosslinks_xinet"] = None
        # check what is uploaded and set
        if uploaded_file is None:
            _ = st.error("You need to upload a result file first!")
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
            uploaded_file is not None
            and search_engine is not None
            and crosslinker is not None
            and crosslinker != "Custom"
        ):
            with st.spinner("Parsing file...", show_time=True):
                try:
                    st.session_state["pr"] = read_file(
                        uploaded_file,
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
                                if group_by  # pyright: ignore[reportPossiblyUnboundVariable]
                                == "Peptide sequence and peptide crosslink position"
                                else "protein",
                                score="higher_better"
                                if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable]
                                else "lower_better",
                            )
                    if validate_aggregated:
                        if st.session_state["aggregated"] is not None:
                            st.session_state["aggregated"] = transform.validate(
                                st.session_state["aggregated"],
                                fdr=fdr,  # pyright: ignore[reportPossiblyUnboundVariable]
                                formula=formula,  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType]
                                score="higher_better"
                                if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable]
                                else "lower_better",
                                separate_intra_inter=separate  # pyright: ignore[reportPossiblyUnboundVariable]
                                == "Separate FDR for intra and inter matches",
                            )
                    if unique:
                        st.session_state["pr"] = transform.unique(
                            st.session_state["pr"],
                            by="peptide"
                            if group_by  # pyright: ignore[reportPossiblyUnboundVariable]
                            == "Peptide sequence and peptide crosslink position"
                            else "protein",
                            score="higher_better"
                            if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable]
                            else "lower_better",
                        )
                    if validate:
                        st.session_state["pr"] = transform.validate(
                            st.session_state["pr"],
                            fdr=fdr,  # pyright: ignore[reportPossiblyUnboundVariable]
                            formula=formula,  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType]
                            score="higher_better"
                            if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable]
                            else "lower_better",
                            separate_intra_inter=separate  # pyright: ignore[reportPossiblyUnboundVariable]
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
                        if uploaded_fasta is None:  # pyright: ignore[reportPossiblyUnboundVariable]
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
                            uploaded_fasta,  # pyright: ignore[reportPossiblyUnboundVariable]
                        )
                        if (
                            "aggregated" in st.session_state
                            and st.session_state["aggregated"] is not None
                        ):
                            st.session_state["aggregated"] = reannotating_positions(
                                st.session_state["aggregated"],
                                uploaded_fasta,  # pyright: ignore[reportPossiblyUnboundVariable]
                            )
                except Exception as e:
                    _ = st.error(
                        "Something went wrong! This is most likely due to missing information in the results!",
                        icon="⚠️",
                    )
                    with st.expander("Show exception"):
                        _ = st.exception(e)
        elif (
            uploaded_file is not None
            and search_engine is not None
            and crosslinker is not None
            and crosslinker == "Custom"
            and crosslinker_name is not None
            and crosslinker_mass is not None
        ):
            with st.spinner("Parsing file...", show_time=True):
                try:
                    st.session_state["pr"] = read_file(
                        uploaded_file,
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
                                if group_by  # pyright: ignore[reportPossiblyUnboundVariable]
                                == "Peptide sequence and peptide crosslink position"
                                else "protein",
                                score="higher_better"
                                if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable]
                                else "lower_better",
                            )
                    if validate_aggregated:
                        if st.session_state["aggregated"] is not None:
                            st.session_state["aggregated"] = transform.validate(
                                st.session_state["aggregated"],
                                fdr=fdr,  # pyright: ignore[reportPossiblyUnboundVariable]
                                formula=formula,  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType]
                                score="higher_better"
                                if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable]
                                else "lower_better",
                                separate_intra_inter=separate  # pyright: ignore[reportPossiblyUnboundVariable]
                                == "Separate FDR for intra and inter matches",
                            )
                    if unique:
                        st.session_state["pr"] = transform.unique(
                            st.session_state["pr"],
                            by="peptide"
                            if group_by  # pyright: ignore[reportPossiblyUnboundVariable]
                            == "Peptide sequence and peptide crosslink position"
                            else "protein",
                            score="higher_better"
                            if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable]
                            else "lower_better",
                        )
                    if validate:
                        st.session_state["pr"] = transform.validate(
                            st.session_state["pr"],
                            fdr=fdr,  # pyright: ignore[reportPossiblyUnboundVariable]
                            formula=formula,  # pyright: ignore[reportPossiblyUnboundVariable, reportArgumentType]
                            score="higher_better"
                            if score == "Higher better"  # pyright: ignore[reportPossiblyUnboundVariable]
                            else "lower_better",
                            separate_intra_inter=separate  # pyright: ignore[reportPossiblyUnboundVariable]
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
                        if uploaded_fasta is None:  # pyright: ignore[reportPossiblyUnboundVariable]
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
                            uploaded_fasta,  # pyright: ignore[reportPossiblyUnboundVariable]
                        )
                        if (
                            "aggregated" in st.session_state
                            and st.session_state["aggregated"] is not None
                        ):
                            st.session_state["aggregated"] = reannotating_positions(
                                st.session_state["aggregated"],
                                uploaded_fasta,  # pyright: ignore[reportPossiblyUnboundVariable]
                            )
                except Exception as e:
                    _ = st.error(
                        "Something went wrong! This is most likely due to missing information in the results!",
                        icon="⚠️",
                    )
                    with st.expander("Show exception"):
                        _ = st.exception(e)

    if "pr" in st.session_state:
        if st.session_state["pr"]["crosslink-spectrum-matches"] is not None:
            csms_header = st.subheader(
                "Read Crosslink-Spectrum-Matches", divider="grey"
            )
            csms = st.session_state["pr"]["crosslink-spectrum-matches"]
            csms_info = st.markdown(f"**Read {len(csms)} crosslink-spectrum-matches:**")
            csms_df = st.dataframe(
                transform.to_dataframe(csms), use_container_width=True
            )
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
                    use_container_width=True,
                    help="Download crosslink-spectrum-matches in comma-separated format.",
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
                    use_container_width=True,
                    help="Download crosslink-spectrum-matches in Microsoft Excel format.",
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
                    use_container_width=True,
                    help="Download crosslink-spectrum-matches in JavaScript Object Notation (JSON) format.",
                )

        if st.session_state["pr"]["crosslinks"] is not None:
            crosslinks_header = st.subheader("Read Crosslinks", divider="grey")
            crosslinks = st.session_state["pr"]["crosslinks"]
            crosslinks_info = st.markdown(f"**Read {len(crosslinks)} crosslinks:**")
            crosslinks_df = st.dataframe(
                transform.to_dataframe(crosslinks), use_container_width=True
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
                    use_container_width=True,
                    help="Download crosslinks in comma-separated format.",
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
                    use_container_width=True,
                    help="Download crosslinks in Microsoft Excel format.",
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
                    use_container_width=True,
                    help="Download crosslinks in JavaScript Object Notation (JSON) format.",
                )

    if "aggregated" in st.session_state and st.session_state["aggregated"] is not None:
        aggregated_crosslinks_header = st.subheader(
            "Aggregated Crosslinks", divider="grey"
        )
        aggregated_crosslinks = st.session_state["aggregated"]
        aggregated_crosslinks_info = st.markdown(
            f"**Aggregated {len(aggregated_crosslinks)} crosslinks:**"
        )
        aggregated_crosslinks_df = st.dataframe(
            transform.to_dataframe(aggregated_crosslinks), use_container_width=True
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
                use_container_width=True,
                help="Download aggregated crosslinks in comma-separated format.",
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
                use_container_width=True,
                help="Download aggregated crosslinks in Microsoft Excel format.",
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
                use_container_width=True,
                help="Download aggregated crosslinks in JavaScript Object Notation (JSON) format.",
            )


def visualize_tab():
    if "pr" not in st.session_state and "aggregated" not in st.session_state:
        no_data = st.info("You need to upload a result file first!")
    if "pr" in st.session_state and st.session_state["pr"] is not None:
        if st.session_state["pr"]["crosslink-spectrum-matches"] is not None:
            csms = st.session_state["pr"]["crosslink-spectrum-matches"]
            available_keys = transform.get_available_keys(csms)
            if (
                available_keys["score"]
                and available_keys["alpha_decoy"]
                and available_keys["beta_decoy"]
            ):
                csms_score_dist_header = st.subheader(
                    "Crosslink-Spectrum-Match Score Distribution", divider="grey"
                )
                fig, ax = plotting.plot_score_distribution(csms)
                csms_score_dist = st.pyplot(fig, use_container_width=True)
            else:
                csms_not_enough_data = st.info(
                    "Not enough data to plot score distribution for crosslink-spectrum-matches!"
                )
        if st.session_state["pr"]["crosslinks"] is not None:
            crosslinks = st.session_state["pr"]["crosslinks"]
            available_keys = transform.get_available_keys(crosslinks)
            if (
                available_keys["score"]
                and available_keys["alpha_decoy"]
                and available_keys["beta_decoy"]
            ):
                crosslinks_score_dist_header = st.subheader(
                    "Crosslink Score Distribution", divider="grey"
                )
                fig, ax = plotting.plot_score_distribution(crosslinks)
                crosslinks_score_dist = st.pyplot(fig, use_container_width=True)
            else:
                crosslinks_not_enough_data = st.info(
                    "Not enough data to plot score distribution for crosslinks!"
                )
    if "aggregated" in st.session_state and st.session_state["aggregated"] is not None:
        aggregated_crosslinks = st.session_state["aggregated"]
        available_keys = transform.get_available_keys(aggregated_crosslinks)
        if (
            available_keys["score"]
            and available_keys["alpha_decoy"]
            and available_keys["beta_decoy"]
        ):
            aggregated_crosslinks_score_dist_header = st.subheader(
                "Aggregated Crosslink Score Distribution", divider="grey"
            )
            fig, ax = plotting.plot_score_distribution(aggregated_crosslinks)
            aggregated_crosslinks_score_dist = st.pyplot(fig, use_container_width=True)
        else:
            aggregated_crosslinks_not_enough_data = st.info(
                "Not enough data to plot score distribution for aggregated crosslinks!"
            )


def export_tab():
    if "pr" not in st.session_state and "aggregated" not in st.session_state:
        no_data = st.info("You need to upload a result file first!")
    if "pr" in st.session_state and st.session_state["pr"] is not None:
        if st.session_state["pr"]["crosslink-spectrum-matches"] is not None:
            csms = st.session_state["pr"]["crosslink-spectrum-matches"]
            export_csms_header = st.subheader(
                "Export Crosslink-Spectrum-Matches", divider="grey"
            )
            export_csms_options = ["IMP-X-FDR", "MS Annika", "xiFDR"]
            export_csms_picker = st.selectbox(
                "Export crosslink-spectrum-matches to:",
                options=export_csms_options,
                index=None,
                help="Chose a format to export the crosslink-spectrum-matches to.",
            )
            if export_csms_picker is None:
                pass
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
                    use_container_width=True,
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
                        use_container_width=True,
                        help="Downloads the exported crosslink-spectrum-matches in IMP-X-FDR format.",
                        key="export_csms_impxfdr_download",
                    )
            elif export_csms_picker == "MS Annika":
                export_csms_msannika_button = st.button(
                    "Export to MS Annika format!",
                    type="primary",
                    use_container_width=True,
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
                            use_container_width=True,
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
                            use_container_width=True,
                            help="Downloads the exported crosslink-spectrum-matches in MS Annika Microsoft Excel (.xlsx) format.",
                            key="export_csms_msannika_download_xlsx",
                        )
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
                    use_container_width=True,
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
                        use_container_width=True,
                        help="Downloads the exported crosslink-spectrum-matches in xiFDR format.",
                        key="export_csms_xifdr_download",
                    )
            else:
                pass

        if st.session_state["pr"]["crosslinks"] is not None:
            crosslinks = st.session_state["pr"]["crosslinks"]
            export_crosslinks_header = st.subheader("Export Crosslinks", divider="grey")
            export_crosslinks_options = [
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
                help="Chose a format to export the crosslinks to.",
            )
            if export_crosslinks_picker is None:
                pass
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
                    use_container_width=True,
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
                        use_container_width=True,
                        help="Downloads the exported crosslinks in IMP-X-FDR format.",
                        key="export_crosslinks_impxfdr_download",
                    )
            elif export_crosslinks_picker == "MS Annika":
                export_crosslinks_msannika_button = st.button(
                    "Export to MS Annika format!",
                    type="primary",
                    use_container_width=True,
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
                            use_container_width=True,
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
                            use_container_width=True,
                            help="Downloads the exported crosslinks in MS Annika Microsoft Excel (.xlsx) format.",
                            key="export_crosslinks_msannika_download_xlsx",
                        )
            elif export_crosslinks_picker == "PyXlinkViewer":
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
                export_crosslinks_pyxlinkviewer_info = st.info(
                    "To export to PyXlinkViewer your crosslinks should be **unique** and **not** "
                    + "**contain any decoy matches**! Usually you would also want to filter for high-confidence crosslinks! "
                    + "You can check this in the "
                    + "**'Load Data'** tab in the **'Summary Statistics'** of your loaded result!"
                )
                export_crosslinks_pyxlinkviewer_button = st.button(
                    "Export to PyXlinkViewer format!",
                    type="primary",
                    use_container_width=True,
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
                    export_crosslinks_msannika_download_info = st.markdown(
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
                        use_container_width=True,
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
                            use_container_width=True,
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
                            use_container_width=True,
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
                            use_container_width=True,
                            help="Downloads the parsed PDB annotation.",
                            key="export_crosslinks_pyxlinkviewer_download_meta_pdb_annotation",
                        )
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
                    use_container_width=True,
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
                        use_container_width=True,
                        help="Downloads the exported crosslinks in xiNET format.",
                        key="export_crosslinks_xinet_download",
                    )
            elif export_crosslinks_picker == "xiVIEW":
                pass
            elif export_crosslinks_picker == "XlinkDB":
                pass
            elif export_crosslinks_picker == "xlms-tools":
                pass
            elif export_crosslinks_picker == "XMAS":
                pass
            else:
                pass

    if "aggregated" in st.session_state and st.session_state["aggregated"] is not None:
        aggregated_crosslinks = st.session_state["aggregated"]
        export_aggregated_crosslinks_header = st.subheader(
            "Export Aggregated Crosslinks", divider="grey"
        )


def about_tab():
    general_description = """
        **pyXLMS** is a python package and web application with graphical user interface that aims to simplify and streamline the intermediate step of
        connecting crosslink search engine results with down-stream analysis tools, enabling researchers even without bioinformatics knowledge to
        conduct in-depth crosslink analyses and shifting the focus from data transformation to data interpretation and therefore gaining biological
        insight.

        Currently pyXLMS supports input from six different crosslink search engines:
        [MaxLynx (part of MaxQuant)](https://www.maxquant.org/),
        [MS Annika](https://github.com/hgb-bin-proteomics/MSAnnika),
        [pLink 2 and pLink 3](http://pfind.ict.ac.cn/se/plink/),
        [Scout](https://github.com/diogobor/Scout),
        [xiSearch](https://www.rappsilberlab.org/software/xisearch/) and [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
        [XlinkX](https://docs.thermofisher.com/r/XlinkX-3.2-Quick-Start-Guide/),
        as well as the [mzIdentML format](https://www.psidev.info/mzidentml) of the HUPO Proteomics Standards Initiative,
        and a well-documented and [human-readable custom tabular format](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md).

        Down-stream analysis is facilitated by functionality that is directly available within pyXLMS such as validation, annotation, aggregation,
        and filtering of crosslink-spectrum-matches and crosslinks.

        In addition, the data can easily be exported to the required data format of the various available down-stream analysis tools such as
        [xiNET](https://crosslinkviewer.org/index.php),
        [xiVIEW](https://www.xiview.org/index.php),
        [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
        [XlinkDB](https://xlinkdb.gs.washington.edu/xlinkdb/),
        [xlms-tools](https://gitlab.com/topf-lab/xlms-tools),
        pyMOL (via [pyXlinkViewer](https://github.com/BobSchiffrin/PyXlinkViewer)),
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

    header_3 = st.header("Contact", divider="grey")
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
            use_container_width=True,
        )

    with center_1:
        link_button_2 = st.link_button(
            "User Guide",
            url="https://hgb-bin-proteomics.github.io/pyXLMS-docs",
            type="primary",
            help="Link to the pyXLMS user guide page.",
            use_container_width=True,
        )

    with r1:
        link_button_3 = st.link_button(
            "Documentation",
            url="https://hgb-bin-proteomics.github.io/pyXLMS",
            type="primary",
            help="Link to the pyXLMS documentaion page.",
            use_container_width=True,
        )


# main page content
def main_page():
    title = st.title("pyXLMS")

    tab1, tab2, tab3, tab4 = st.tabs(["Load Data", "Visualize", "Export", "About"])

    with tab1:
        input_tab()

    with tab2:
        visualize_tab()

    with tab3:
        export_tab()

    with tab4:
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
        initial_sidebar_state="expanded",
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
    info_str += "- **Documentation:**  \n [hgb-bin-proteomics.github.io/pyXLMS-docs/](hgb-bin-proteomics.github.io/pyXLMS-docs/)\n"
    info_str += "- **Contact:**  \n  [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)\n"
    info_str += "- **License:**  \n  [MIT License](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/LICENSE)\n"
    info_str += "- **Project Page:**  \n  [GitHub](https://github.com/hgb-bin-proteomics/pyXLMS/)"
    info = st.sidebar.markdown(info_str)

    main_page()


if __name__ == "__main__":
    main()
