#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

AMINO_ACIDS = {
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "U",
    "W",
    "Y",
    "V",
    "J",
    "O",
}
"""List of valid amino acids.

List of one-letter codes for all valid amino acids.

Examples
--------
>>> from pyXLMS.constants import AMINO_ACIDS
>>> "A" in AMINO_ACIDS
True
>>> "B" in AMINO_ACIDS
False
"""

CROSSLINKERS = {
    "BS3": 138.06808,
    "DSS": 138.06808,
    "DSSO": 158.00376,
    "ADH": 138.09054635,
    "DSBSO": 308.03883,
    "PhoX": 209.97181,
}
"""Dictionary of crosslinkers.

Dictionary of pre-defined crosslinkers that maps crosslinker names to crosslinker delta masses.
Currently contains `"BS3"`, `"DSS"`, `"DSSO"`, `"ADH"`, `"DSBSO"`, `"PhoX"`.

Examples
--------
>>> from pyXLMS.constants import CROSSLINKERS
>>> CROSSLINKERS["BS3"]
138.06808
"""

MODIFICATIONS = {
    "Carbamidomethyl": 57.021464,
    "Oxidation": 15.994915,
    "Phospho": 79.966331,
    "Acetyl": 42.010565,
}
"""Dictionary of post-translational-modifications.

Dictionary of pre-defined post-translational-modifications that maps modification names to modification delta masses.
Currently contains `"Carbamidomethyl"`, `"Oxidation"`, `"Phospho"`, `"Acetyl"` and all crosslinkers.

Examples
--------
>>> from pyXLMS.constants import MODIFICATIONS
>>> MODIFICATIONS["Carbamidomethyl"]
57.021464
>>> MODIFICATIONS["BS3"]
138.06808
"""

MODIFICATIONS.update(CROSSLINKERS)
