#!/usr/bin/env python3

# SCRIPT NAME GUI
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

import streamlit as st
from main import my_product

# main page content
def main_page():

    title = st.title("TITLE")

    general_description = \
    """
    DESCRIPTION
    """
    description = st.markdown(general_description)

    header_1 = st.subheader("Product", divider = "rainbow")

    f1 = st.number_input("Number 1:",
                         min_value = 1,
                         max_value = 10000,
                         value = 1,
                         step = 1,
                         help = "First factor of multiplication.")

    f2 = st.number_input("Number 2:",
                         min_value = 1,
                         max_value = 10000,
                         value = 2,
                         step = 1,
                         help = "Second factor of multiplication.")

    l1, l2, center_button, r1, r2 = st.columns(5)

    with center_button:
        compute = st.button("Compute!", type = "primary")

    if compute:
        st.success(f"The product of {f1} and {f2} is {my_product(f1, f2)}")

# side bar and main page loader
def main():

    about_str = \
    """
    ABOUT
    """

    st.set_page_config(page_title = "TITLE",
                       page_icon = ":test_tube:",
                       layout = "wide",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": "https://github.com/YOUR_REPO/discussions",
                                     "Report a bug": "https://github.com/YOUR_REPO/issues",
                                     "About": about_str}
                       )

    title = st.sidebar.title("TITLE")

    logo = st.sidebar.image("https://user-images.githubusercontent.com/7164864/217935870-c0bc60a3-6fc0-4047-b011-7b4c59488c91.png",
                            caption = "CAPTION")

    doc = st.sidebar.markdown(about_str)

    citing_str = "**Citing:** CITATION INFO"
    citing = st.sidebar.markdown(citing_str)

    contact_str = "**Contact:** CONTACT INFO"
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** [MIT License](https://github.com/YOUR_REPO/blob/master/LICENSE.md)"
    license = st.sidebar.markdown(license_str)

    project_str = "**Project Page:** [GitHub](https://github.com/YOUR_REPO/)"
    project = st.sidebar.markdown(project_str)

    main_page()

if __name__ == "__main__":

    main()
