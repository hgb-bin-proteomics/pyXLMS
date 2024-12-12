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
CROSSLINKERS = {"BS3": 138.06808, "DSS": 138.06808, "DSSO": 158.00376,"ADH": 138.09054635, "DSBSO": 308.03883, "PhoX": 209.97181}
MODIFICATIONS = {"Carbamidomethyl": 57.021464, "Oxidation": 15.994915, "Phospho": 79.966331, "Acetyl": 42.010565}
MODIFICATIONS.update(CROSSLINKERS)
