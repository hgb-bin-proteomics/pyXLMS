#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

PRIDE_URL = "https://www.ebi.ac.uk/pride/archive/projects/PXD014337"
PROTEOMEXCHANGE_URL = (
    "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014337"
)
PUBLICATION_URL = "https://doi.org/10.1038/s41467-020-14608-2"
DATASETS = [
    {
        "files": [
            {
                "filename": "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw",
            },
            {
                "filename": "XLpeplib_Beveridge_QEx-HFX_DSS_R2.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSS_R2.raw",
            },
            {
                "filename": "XLpeplib_Beveridge_QEx-HFX_DSS_R3.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSS_R3.raw",
            },
        ],
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "5 ppm",
            "MS2 Tolerance": "20 ppm",
            "Enzyme": "Trypsin",
            "Cleavage": "K,R blocked by P",
            "Maximum Missed Cleavages": "3",
            "Crosslink Modification": "DSS[K]+138.06808",
            "Fixed Modifications": "Carbamidomethyl[C]+57.021464",
            "Variable Modifications": "Oxidation[M]+15.994915",
            "Minimum Peptide Length": "5",
            "Maximum Peptide Length": "60",
        },
    },
    {
        "files": [
            {
                "filename": "XLpeplib_Beveridge_QEx-HFX_DSSO_stHCD.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSSO_stHCD.raw",
            },
        ],
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "",
            "MS2 Tolerance": "",
            "Enzyme": "",
            "Cleavage": "",
            "Maximum Missed Cleavages": "",
            "Crosslink Modification": "",
            "Fixed Modifications": "",
            "Variable Modifications": "",
            "Minimum Peptide Length": "",
            "Maximum Peptide Length": "",
        },
    },
    {
        "files": [
            {
                "filename": "XLpeplib_Beveridge_QEx-HFX_DSBU_stHCD.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSBU_stHCD.raw",
            },
        ],
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "",
            "MS2 Tolerance": "",
            "Enzyme": "",
            "Cleavage": "",
            "Maximum Missed Cleavages": "",
            "Crosslink Modification": "",
            "Fixed Modifications": "",
            "Variable Modifications": "",
            "Minimum Peptide Length": "",
            "Maximum Peptide Length": "",
        },
    },
    {
        "files": [
            {
                "filename": "XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw",
            },
        ],
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "",
            "MS2 Tolerance": "",
            "Enzyme": "",
            "Cleavage": "",
            "Maximum Missed Cleavages": "",
            "Crosslink Modification": "",
            "Fixed Modifications": "",
            "Variable Modifications": "",
            "Minimum Peptide Length": "",
            "Maximum Peptide Length": "",
        },
    },
    {
        "files": [
            {
                "filename": "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_Lumos_DSSO_MS3.raw",
            },
        ],
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "",
            "MS2 Tolerance": "",
            "MS3 Tolerance": "",
            "Enzyme": "",
            "Cleavage": "",
            "Maximum Missed Cleavages": "",
            "Crosslink Modification": "",
            "Fixed Modifications": "",
            "Variable Modifications": "",
            "Minimum Peptide Length": "",
            "Maximum Peptide Length": "",
        },
    },
    {
        "files": [
            {
                "filename": "XLpeplib_Beveridge_Lumos_DSSO_MS3-EThcD.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_Lumos_DSSO_MS3-EThcD.raw",
            },
        ],
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "",
            "MS2 Tolerance": "",
            "MS3 Tolerance": "",
            "Enzyme": "",
            "Cleavage": "",
            "Maximum Missed Cleavages": "",
            "Crosslink Modification": "",
            "Fixed Modifications": "",
            "Variable Modifications": "",
            "Minimum Peptide Length": "",
            "Maximum Peptide Length": "",
        },
    },
    {
        "files": [
            {
                "filename": "XLpeplib_Beveridge_Lumos_DSSO_CID-ETD.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_Lumos_DSSO_CID-ETD.raw",
            },
        ],
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "",
            "MS2 Tolerance": "",
            "Enzyme": "",
            "Cleavage": "",
            "Maximum Missed Cleavages": "",
            "Crosslink Modification": "",
            "Fixed Modifications": "",
            "Variable Modifications": "",
            "Minimum Peptide Length": "",
            "Maximum Peptide Length": "",
        },
    },
]
LIBRARIES = {"standard": {"PEPTIDE1": {"PEPTIDE1", "PEPTIDE2", "PEPTIDE3"}}}
