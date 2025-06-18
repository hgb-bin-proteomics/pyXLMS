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
r"""List of valid amino acids.

List of one-letter codes for all valid amino acids.

Examples
--------
>>> from pyXLMS.constants import AMINO_ACIDS
>>> "A" in AMINO_ACIDS
True
>>> "B" in AMINO_ACIDS
False
"""

AMINO_ACIDS_3TO1 = {
    "GLY": "G",
    "PRO": "P",
    "ALA": "A",
    "VAL": "V",
    "LEU": "L",
    "ILE": "I",
    "MET": "M",
    "CYS": "C",
    "PHE": "F",
    "TYR": "Y",
    "TRP": "W",
    "HIS": "H",
    "LYS": "K",
    "ARG": "R",
    "GLN": "Q",
    "ASN": "N",
    "GLU": "E",
    "ASP": "D",
    "SER": "S",
    "THR": "T",
}
r"""Mapping of amino acid 3-letter codes to their 1-letter codes.

Mapping of all amino acid 3-letter codes to their corresponding 1-letter codes.

Examples
--------
>>> from pyXLMS.constants import AMINO_ACIDS_3TO1
>>> AMINO_ACIDS_3TO1["GLY"]
'G'
"""

AMINO_ACIDS_1TO3 = {v: k for k, v in AMINO_ACIDS_3TO1.items()}
r"""Mapping of amino acid 1-letter codes to their 3-letter codes.

Mapping of all amino acid 1-letter codes to their corresponding 3-letter codes.

Examples
--------
>>> from pyXLMS.constants import AMINO_ACIDS_1TO3
>>> AMINO_ACIDS_1TO3["G"]
'GLY'
"""

CROSSLINKERS = {
    "BS3": 138.06808,
    "DSS": 138.06808,
    "DSSO": 158.00376,
    "DSBU": 196.08479231,
    "ADH": 138.09054635,
    "DSBSO": 308.03883,
    "PhoX": 209.97181,
}
r"""Dictionary of crosslinkers.

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
r"""Dictionary of post-translational-modifications.

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

XI_MODIFICATION_MAPPING = {
    "->": ("Substitution", float("nan")),
    "cm": ("Carbamidomethyl", 57.021464),
    "ox": ("Oxidation", 15.994915),
    "bs3oh": ("BS3 Hydrolized", 156.0786347),
    "bs3nh2": ("BS3 Amidated", 155.094619105),
    "bs3loop": ("BS3 Looplink", 138.06808),
    "bs3_hyd": ("BS3 Hydrolized", 156.0786347),
    "bs3_ami": ("BS3 Amidated", 155.094619105),
    "bs3_tris": ("BS3 Tris", 259.141973),
    "dssoloop": ("DSSO Looplink", 158.00376),
    "dsso_loop": ("DSSO Looplink", 158.00376),
    "dsso_hyd": ("DSSO Hydrolized", 176.0143295),
    "dsso_ami": ("DSSO Amidated", 175.030313905),
    "dsso_tris": ("DSSO Tris", 279.077658),
    "dsbuloop": ("DSBU Looplink", 196.08479231),
    "dsbu_loop": ("DSBU Looplink", 196.08479231),
    "dsbu_hyd": ("DSBU Hydrolized", 214.095357),
    "dsbu_ami": ("DSBU Amidated", 213.111341),
    "dsbu_tris": ("DSBU Tris", 317.158685),
}
r"""Dictionary that maps sequence elements from xiSearch and xiFDR to their corresponding post-translational-modifications.

Dictionary that maps sequence elements (e.g. `"cm"`) from xiSearch and xiFDR to their corresponding
post-translational-modifications (e.g. `("Carbamidomethyl", 57.021464)`).

Examples
--------
>>> from pyXLMS.constants import XI_MODIFICATION_MAPPING
>>> XI_MODIFICATION_MAPPING["cm"]
('Carbamidomethyl', 57.021464)
>>> XI_MODIFICATION_MAPPING["ox"]
('Oxidation', 15.994915)
"""

SCOUT_MODIFICATION_MAPPING = {
    "+57.021460": ("Carbamidomethyl", 57.021464),
    "+15.994900": ("Oxidation", 15.994915),
    "Oxidation of Methionine": ("Oxidation", 15.994915),
    "Carbamidomethyl": ("Carbamidomethyl", 57.021464),
}
r"""Dictionary that maps sequence elements and modifications from Scout to their corresponding post-translational-modifications.

Dictionary that maps sequence elements (e.g. `"+57.021460"`) and modifications (e.g. `"Carbamidomethyl"`) from Scout to their
corresponding post-translational-modifications (e.g. `("Carbamidomethyl", 57.021464)`).

Examples
--------
>>> from pyXLMS.constants import SCOUT_MODIFICATION_MAPPING
>>> SCOUT_MODIFICATION_MAPPING["+57.021460"]
('Carbamidomethyl', 57.021464)
>>> SCOUT_MODIFICATION_MAPPING["Carbamidomethyl"]
('Carbamidomethyl', 57.021464)
>>> SCOUT_MODIFICATION_MAPPING["Oxidation of Methionine"]
('Oxidation', 15.994915)
"""

SCOUT_CROSSLINKER_MAPPING = {k: (k, v) for k, v in CROSSLINKERS.items()}
SCOUT_MODIFICATION_MAPPING.update(SCOUT_CROSSLINKER_MAPPING)
