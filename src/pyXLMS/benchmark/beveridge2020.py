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
        "description": "Triplicate measurement of Cas9 crosslinked with DSS.",
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
        "fasta": {
            "filename": "Cas9_plus10.fasta",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14608-2/MediaObjects/41467_2020_14608_MOESM4_ESM.zip",
            "alt_url": "https://github.com/hgb-bin-proteomics/MSAnnika_NC_Results/raw/master/Peplib_Beveridge/Cas9_plus10.fasta",
        },
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
        "description": "",
        "files": [
            {
                "filename": "XLpeplib_Beveridge_QEx-HFX_DSSO_stHCD.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSSO_stHCD.raw",
            },
        ],
        "fasta": {
            "filename": "cas9_crapome.fasta",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14608-2/MediaObjects/41467_2020_14608_MOESM4_ESM.zip",
            "alt_url": "https://github.com/hgb-bin-proteomics/MSAnnika_MS3_Results/raw/master/PXD014337/cas9_crapome.fasta",
        },
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "5 ppm",
            "MS2 Tolerance": "20 ppm",
            "Enzyme": "Trypsin",
            "Cleavage": "K,R blocked by P",
            "Maximum Missed Cleavages": "3",
            "Crosslink Modification": "DSSO/+158.004",
            "Fixed Modifications": "Carbamidomethyl/+57.021",
            "Variable Modifications": "Oxidation/+15.995",
            "Minimum Peptide Length": "5",
            "Maximum Peptide Length": "60",
        },
    },
    {
        "description": "",
        "files": [
            {
                "filename": "XLpeplib_Beveridge_QEx-HFX_DSBU_stHCD.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSBU_stHCD.raw",
            },
        ],
        "fasta": {
            "filename": "cas9_crapome.fasta",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14608-2/MediaObjects/41467_2020_14608_MOESM4_ESM.zip",
            "alt_url": "https://github.com/hgb-bin-proteomics/MSAnnika_MS3_Results/raw/master/PXD014337/cas9_crapome.fasta",
        },
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "5 ppm",
            "MS2 Tolerance": "20 ppm",
            "Enzyme": "Trypsin",
            "Cleavage": "K,R blocked by P",
            "Maximum Missed Cleavages": "3",
            "Crosslink Modification": "DSSO/+158.004",
            "Fixed Modifications": "Carbamidomethyl/+57.021",
            "Variable Modifications": "Oxidation/+15.995",
            "Minimum Peptide Length": "5",
            "Maximum Peptide Length": "60",
        },
    },
    {
        "description": "",
        "files": [
            {
                "filename": "XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_Lumos_DSSO_stHCD-MS2.raw",
            },
        ],
        "fasta": {
            "filename": "cas9_crapome.fasta",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14608-2/MediaObjects/41467_2020_14608_MOESM4_ESM.zip",
            "alt_url": "https://github.com/hgb-bin-proteomics/MSAnnika_MS3_Results/raw/master/PXD014337/cas9_crapome.fasta",
        },
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "5 ppm",
            "MS2 Tolerance": "20 ppm",
            "Enzyme": "Trypsin",
            "Cleavage": "K,R blocked by P",
            "Maximum Missed Cleavages": "3",
            "Crosslink Modification": "DSSO/+158.004",
            "Fixed Modifications": "Carbamidomethyl/+57.021",
            "Variable Modifications": "Oxidation/+15.995",
            "Minimum Peptide Length": "5",
            "Maximum Peptide Length": "60",
        },
    },
    {
        "description": "",
        "files": [
            {
                "filename": "XLpeplib_Beveridge_Lumos_DSSO_MS3.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_Lumos_DSSO_MS3.raw",
            },
        ],
        "fasta": {
            "filename": "cas9_crapome.fasta",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14608-2/MediaObjects/41467_2020_14608_MOESM4_ESM.zip",
            "alt_url": "https://github.com/hgb-bin-proteomics/MSAnnika_MS3_Results/raw/master/PXD014337/cas9_crapome.fasta",
        },
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "5 ppm",
            "MS2 Tolerance": "20 ppm",
            "MS3 Tolerance": "",
            "Enzyme": "Trypsin",
            "Cleavage": "K,R blocked by P",
            "Maximum Missed Cleavages": "3",
            "Crosslink Modification": "DSSO/+158.004",
            "Fixed Modifications": "Carbamidomethyl/+57.021",
            "Variable Modifications": "Oxidation/+15.995",
            "Minimum Peptide Length": "5",
            "Maximum Peptide Length": "60",
        },
    },
    {
        "description": "",
        "files": [
            {
                "filename": "XLpeplib_Beveridge_Lumos_DSSO_MS3-EThcD.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_Lumos_DSSO_MS3-EThcD.raw",
            },
        ],
        "fasta": {
            "filename": "cas9_crapome.fasta",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14608-2/MediaObjects/41467_2020_14608_MOESM4_ESM.zip",
            "alt_url": "https://github.com/hgb-bin-proteomics/MSAnnika_MS3_Results/raw/master/PXD014337/cas9_crapome.fasta",
        },
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "5 ppm",
            "MS2 Tolerance": "20 ppm",
            "MS3 Tolerance": "",
            "Enzyme": "Trypsin",
            "Cleavage": "K,R blocked by P",
            "Maximum Missed Cleavages": "3",
            "Crosslink Modification": "DSSO/+158.004",
            "Fixed Modifications": "Carbamidomethyl/+57.021",
            "Variable Modifications": "Oxidation/+15.995",
            "Minimum Peptide Length": "5",
            "Maximum Peptide Length": "60",
        },
    },
    {
        "description": "",
        "files": [
            {
                "filename": "XLpeplib_Beveridge_Lumos_DSSO_CID-ETD.raw",
                "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD014337/XLpeplib_Beveridge_Lumos_DSSO_CID-ETD.raw",
            },
        ],
        "fasta": {
            "filename": "cas9_crapome.fasta",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14608-2/MediaObjects/41467_2020_14608_MOESM4_ESM.zip",
            "alt_url": "https://github.com/hgb-bin-proteomics/MSAnnika_MS3_Results/raw/master/PXD014337/cas9_crapome.fasta",
        },
        "library": "standard",
        "settings": {
            "MS1 Tolerance": "5 ppm",
            "MS2 Tolerance": "20 ppm",
            "Enzyme": "Trypsin",
            "Cleavage": "K,R blocked by P",
            "Maximum Missed Cleavages": "3",
            "Crosslink Modification": "DSSO/+158.004",
            "Fixed Modifications": "Carbamidomethyl/+57.021",
            "Variable Modifications": "Oxidation/+15.995",
            "Minimum Peptide Length": "5",
            "Maximum Peptide Length": "60",
        },
    },
]
LIBRARIES = {"standard": {"PEPTIDE1": {"PEPTIDE1", "PEPTIDE2", "PEPTIDE3"}}}