#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

from os.path import splitext


def __get_filename(filename: str, ext: str) -> str:
    r"""Adds the file extension to the filename if it does not have one already.

    Parameters
    ----------
    filename : str
        The filename.
    ext : str
        The file extension without a dot, for example "csv".

    Returns
    -------
    str
        The filename with guaranteed file extension.
    """
    if splitext(filename)[1].lower() == f".{ext}":
        return filename
    return f"{filename}.{ext}"
