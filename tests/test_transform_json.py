#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from pyXLMS.data import CrosslinkSpectrumMatch, Crosslink, ParserResult
    from pyXLMS.parser import read
    from pyXLMS.transform import to_json, from_json

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    json_data_pr = to_json(pr)
    json_data_csms = to_json(pr.csms())
    _json_data_xls = to_json(pr.xls(), output_file="xls.json")
    pr = from_json(json_data_pr)
    assert isinstance(pr, ParserResult)
    csms = from_json(json_data_csms)
    assert len(csms) == 826
    assert isinstance(csms[0], CrosslinkSpectrumMatch)
    xls = from_json("xls.json")
    assert len(xls) == 300
    assert isinstance(xls[0], Crosslink)


def test2():
    from pyXLMS.data import CrosslinkSpectrumMatch
    from pyXLMS.parser import read
    from pyXLMS.transform import to_json, from_json

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    with open("csms.json", "w", encoding="utf-8") as f:
        _json_data_csms = to_json(pr.csms(), output_file=f)
    with open("csms.json", "r", encoding="utf-8") as f:
        csms = from_json(f)
    assert len(csms) == 826
    assert isinstance(csms[0], CrosslinkSpectrumMatch)
