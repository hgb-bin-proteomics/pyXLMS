#!/usr/bin/env python3

# SCRIPT NAME - TESTS
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

def test1():

    import pandas as pd
    df = pd.read_csv("example_data.csv")

    from main import my_product

    for i, row in df.iterrows():
        assert my_product(row["x"], row["y"]) == row["x"] * row["y"]

def test2():

    from main import main

    assert main(["-f1", "20"]) == 40
    assert main(["-f1", "2", "-f2", "3"]) == 6
