#!/usr/bin/env python3

# pyXLMS GUI - TEST
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file("streamlit_app.py", default_timeout=30.0)
    at.run(timeout=60.0)
    assert not at.exception
