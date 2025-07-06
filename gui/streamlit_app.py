#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "streamlit",
#   "pyxlms>=1",
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

import io
import os
from tempfile import NamedTemporaryFile

from pyXLMS import parser
from pyXLMS import transform
from pyXLMS import constants

import streamlit as st

from typing import Optional
from typing import Dict
from typing import Any


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


# main page content
def main_page():
    title = st.title("pyXLMS")

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

    parse_modifications = st.toggle(
        "Parse modifications",
        key="parse_modifications",
        help="If post-translational-modifications should be parsed or not.",
    )

    crosslinker_name = None
    crosslinker_mass = None
    if crosslinker == "Custom":
        l2, r2 = st.columns(2)

        with l2:
            crosslinker_name = st.text_input(
                "Name of the crosslinker:",
                value=None,
                max_chars=50,
                placeholder="DSSO",
                key="crosslinker_name",
                help="Name of the crosslinker used in the experiment of the uploaded result file.",
            )

        with r2:
            crosslinker_mass = st.number_input(
                "Mass of the crosslinker:",
                value=None,
                step=0.00001,
                format="%0.5f",
                placeholder="158.00376",
                key="crosslinker_mass",
                help="Monoisotopic delta mass of the crosslinker used in the experiment of the uploaded result file.",
            )

    l3, center_3, r3 = st.columns(3)

    with center_3:
        read_file_button = st.button(
            "Read file!", type="primary", use_container_width=True
        )

    if read_file_button:
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
                st.session_state["pr"] = read_file(
                    uploaded_file,
                    search_engine,
                    crosslinker,
                    parse_modifications,
                    crosslinker_mass,
                )
        elif (
            uploaded_file is not None
            and search_engine is not None
            and crosslinker is not None
            and crosslinker == "Custom"
            and crosslinker_name is not None
            and crosslinker_mass is not None
        ):
            with st.spinner("Parsing file...", show_time=True):
                st.session_state["pr"] = read_file(
                    uploaded_file,
                    search_engine,
                    crosslinker_name,
                    parse_modifications,
                    crosslinker_mass,
                )

    if "pr" in st.session_state:
        if st.session_state["pr"]["crosslink-spectrum-matches"] is not None:
            csms = st.session_state["pr"]["crosslink-spectrum-matches"]
            csms_info = st.markdown(f"**Read {len(csms)} crosslink-spectrum-matches:**")
            csms_df = st.dataframe(
                transform.to_dataframe(csms), use_container_width=True
            )
        if st.session_state["pr"]["crosslinks"] is not None:
            crosslinks = st.session_state["pr"]["crosslinks"]
            crosslinks_info = st.markdown(f"**Read {len(crosslinks)} crosslinks:**")
            crosslinks_df = st.dataframe(
                transform.to_dataframe(crosslinks), use_container_width=True
            )


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
