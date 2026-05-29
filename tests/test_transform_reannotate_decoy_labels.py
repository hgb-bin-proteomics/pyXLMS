#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

CSMS = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt"
XLS = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt"
FASTA = "data/_fasta/Cas9_plus10.fasta"


def test1():
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.transform import reannotate_decoy_labels

    xls = [create_crosslink_min("ADANLDK", 7, "GNTDRHSIK", 9)]
    xls = reannotate_decoy_labels(xls, by_mapping={None: False})
    assert not xls[0]["alpha_decoy"]
    assert not xls[0]["beta_decoy"]


def test2():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    csms_2 = reannotate_decoy_labels(
        csms, by_mapping={True: None, False: None, None: None}
    )
    xls_2 = reannotate_decoy_labels(
        xls, by_mapping={True: None, False: None, None: None}
    )

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is not None
        assert csms[i]["beta_decoy"] is not None
        assert csms_2[i]["alpha_decoy"] is None
        assert csms_2[i]["beta_decoy"] is None

    for i in range(len(xls)):
        assert xls[i]["alpha_decoy"] is not None
        assert xls[i]["beta_decoy"] is not None
        assert xls_2[i]["alpha_decoy"] is None
        assert xls_2[i]["beta_decoy"] is None


def test3():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_2 = reannotate_decoy_labels(pr, by_mapping={True: None, False: None, None: None})
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is not None
        assert csms[i]["beta_decoy"] is not None
        assert csms_2[i]["alpha_decoy"] is None
        assert csms_2[i]["beta_decoy"] is None

    for i in range(len(xls)):
        assert xls[i]["alpha_decoy"] is not None
        assert xls[i]["beta_decoy"] is not None
        assert xls_2[i]["alpha_decoy"] is None
        assert xls_2[i]["beta_decoy"] is None


def test4():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_2 = reannotate_decoy_labels(pr, by_mapping={True: None, False: None})
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is not None
        assert csms[i]["beta_decoy"] is not None
        assert csms_2[i]["alpha_decoy"] is None
        assert csms_2[i]["beta_decoy"] is None

    for i in range(len(xls)):
        assert xls[i]["alpha_decoy"] is not None
        assert xls[i]["beta_decoy"] is not None
        assert xls_2[i]["alpha_decoy"] is None
        assert xls_2[i]["beta_decoy"] is None


def test5():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")

    with pytest.raises(
        TypeError,
        match=r"Parameter 'by_mapping' has to be a dictionary that maps bool | None -> bool | None!",
    ):
        _ = reannotate_decoy_labels(
            pr, by_mapping={True: None, False: None, "None": None}
        )


def test6():
    import copy
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = copy.deepcopy(pr)
    csms = list()
    xls = list()
    for csm in pr_none["crosslink-spectrum-matches"]:
        if csm["alpha_decoy"]:
            csm = csm.copy_with_update(
                {
                    "alpha_proteins": [
                        f"REV__{protein}" for protein in csm["alpha_proteins"]
                    ]
                }
            )
        if csm["beta_decoy"]:
            csm = csm.copy_with_update(
                {
                    "beta_proteins": [
                        f"REV__{protein}" for protein in csm["beta_proteins"]
                    ]
                }
            )
        csms.append(csm)
    for xl in pr_none["crosslinks"]:
        if xl["alpha_decoy"]:
            xl = xl.copy_with_update(
                {
                    "alpha_proteins": [
                        f"REV__{protein}" for protein in xl["alpha_proteins"]
                    ]
                }
            )
        if xl["beta_decoy"]:
            xl = xl.copy_with_update(
                {
                    "beta_proteins": [
                        f"REV__{protein}" for protein in xl["beta_proteins"]
                    ]
                }
            )
        xls.append(xl)
    pr_none = pr_none.copy_with_update(
        {"crosslink-spectrum-matches": csms, "crosslinks": xls}
    )
    pr_none = reannotate_decoy_labels(
        pr_none, by_mapping={True: None, False: None, None: None}
    )
    pr_2 = reannotate_decoy_labels(pr_none, by_decoy_protein_prefix="REV__")
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        assert xls[i]["alpha_decoy"] is xls_2[i]["alpha_decoy"]
        assert xls[i]["beta_decoy"] is xls_2[i]["beta_decoy"]


def test7():
    import copy
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = copy.deepcopy(pr)
    csms = list()
    xls = list()
    # previously this test (non pydantic) used whitespace in front of the REV__
    # this doesn't work with pydantic anymore because whitespace is automatically stripped
    for csm in pr_none["crosslink-spectrum-matches"]:
        if csm["alpha_decoy"]:
            csm = csm.copy_with_update(
                {
                    "alpha_proteins": [
                        f"S REV__{protein}" for protein in csm["alpha_proteins"]
                    ]
                }
            )
        if csm["beta_decoy"]:
            csm = csm.copy_with_update(
                {
                    "beta_proteins": [
                        f"S REV__{protein}" for protein in csm["beta_proteins"]
                    ]
                }
            )
        csms.append(csm)
    for xl in pr_none["crosslinks"]:
        if xl["alpha_decoy"]:
            xl = xl.copy_with_update(
                {
                    "alpha_proteins": [
                        f"S REV__{protein}" for protein in xl["alpha_proteins"]
                    ]
                }
            )
        if xl["beta_decoy"]:
            xl = xl.copy_with_update(
                {
                    "beta_proteins": [
                        f"S REV__{protein}" for protein in xl["beta_proteins"]
                    ]
                }
            )
        xls.append(xl)
    pr_none = pr_none.copy_with_update(
        {"crosslink-spectrum-matches": csms, "crosslinks": xls}
    )
    pr_none = reannotate_decoy_labels(
        pr_none, by_mapping={True: None, False: None, None: None}
    )
    pr_2 = reannotate_decoy_labels(pr_none, by_decoy_protein_prefix="REV__")
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert not csms_2[i]["alpha_decoy"]
        assert not csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        assert not xls_2[i]["alpha_decoy"]
        assert not xls_2[i]["beta_decoy"]


def test8():
    import copy
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = copy.deepcopy(pr)
    csms = list()
    xls = list()
    for csm in pr_none["crosslink-spectrum-matches"]:
        if csm["alpha_decoy"]:
            csm = csm.copy_with_update(
                {
                    "alpha_proteins": [
                        f"REV__{protein}" for protein in csm["alpha_proteins"]
                    ]
                }
            )
        if csm["beta_decoy"]:
            csm = csm.copy_with_update(
                {
                    "beta_proteins": [
                        f"REV__{protein}" for protein in csm["beta_proteins"]
                    ]
                }
            )
        csms.append(csm)
    for xl in pr_none["crosslinks"]:
        if xl["alpha_decoy"]:
            xl = xl.copy_with_update(
                {
                    "alpha_proteins": [
                        f"REV__{protein}" for protein in xl["alpha_proteins"]
                    ]
                }
            )
        if xl["beta_decoy"]:
            xl = xl.copy_with_update(
                {
                    "beta_proteins": [
                        f"REV__{protein}" for protein in xl["beta_proteins"]
                    ]
                }
            )
        xls.append(xl)
    pr_none = pr_none.copy_with_update(
        {"crosslink-spectrum-matches": csms, "crosslinks": xls}
    )
    pr_none = reannotate_decoy_labels(
        pr_none, by_mapping={True: None, False: None, None: None}
    )
    pr_2 = reannotate_decoy_labels(pr_none, by_decoy_protein_substring="REV__")
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        assert xls[i]["alpha_decoy"] is xls_2[i]["alpha_decoy"]
        assert xls[i]["beta_decoy"] is xls_2[i]["beta_decoy"]


def test9():
    import copy
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = copy.deepcopy(pr)
    csms = list()
    xls = list()
    for csm in pr_none["crosslink-spectrum-matches"]:
        if csm["alpha_decoy"]:
            csm = csm.copy_with_update(
                {
                    "alpha_proteins": [
                        f"REV__{protein}" for protein in csm["alpha_proteins"]
                    ]
                }
            )
        if csm["beta_decoy"]:
            csm = csm.copy_with_update(
                {
                    "beta_proteins": [
                        f"REV__{protein}" for protein in csm["beta_proteins"]
                    ]
                }
            )
        csms.append(csm)
    for xl in pr_none["crosslinks"]:
        if xl["alpha_decoy"]:
            xl = xl.copy_with_update(
                {
                    "alpha_proteins": [
                        f"REV__{protein}" for protein in xl["alpha_proteins"]
                    ]
                }
            )
        if xl["beta_decoy"]:
            xl = xl.copy_with_update(
                {
                    "beta_proteins": [
                        f"REV__{protein}" for protein in xl["beta_proteins"]
                    ]
                }
            )
        xls.append(xl)
    pr_none = pr_none.copy_with_update(
        {"crosslink-spectrum-matches": csms, "crosslinks": xls}
    )
    pr_none = reannotate_decoy_labels(
        pr_none, by_mapping={True: None, False: None, None: None}
    )
    pr_2 = reannotate_decoy_labels(pr_none, by_decoy_protein_substring="REV__")
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        assert xls[i]["alpha_decoy"] is xls_2[i]["alpha_decoy"]
        assert xls[i]["beta_decoy"] is xls_2[i]["beta_decoy"]


def test10():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = reannotate_decoy_labels(
        pr, by_mapping={True: None, False: None, None: None}
    )
    pr_2 = reannotate_decoy_labels(pr_none, by_target_fasta=FASTA)
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        # because MS Annika does not report decoy labels for individual peptides at XL level
        # we have to do like this
        new_alpha = xls_2[i]["alpha_decoy"] or xls_2[i]["beta_decoy"]
        assert xls[i]["alpha_decoy"] is new_alpha
        new_beta = xls_2[i]["alpha_decoy"] or xls_2[i]["beta_decoy"]
        assert xls[i]["beta_decoy"] is new_beta


def test11():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = reannotate_decoy_labels(
        pr, by_mapping={True: None, False: None, None: None}
    )
    with open(FASTA, "r", encoding="utf-8") as f:
        pr_2 = reannotate_decoy_labels(pr_none, by_target_fasta=f)
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        # because MS Annika does not report decoy labels for individual peptides at XL level
        # we have to do like this
        new_alpha = xls_2[i]["alpha_decoy"] or xls_2[i]["beta_decoy"]
        assert xls[i]["alpha_decoy"] is new_alpha
        new_beta = xls_2[i]["alpha_decoy"] or xls_2[i]["beta_decoy"]
        assert xls[i]["beta_decoy"] is new_beta


def test12():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = reannotate_decoy_labels(
        pr, by_mapping={True: None, False: None, None: None}
    )
    pr_2 = reannotate_decoy_labels(pr_none, by_decoy_fasta=FASTA)
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is not csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is not csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        # because MS Annika does not report decoy labels for individual peptides at XL level
        # we have to do like this
        new_alpha = xls_2[i]["alpha_decoy"] and xls_2[i]["beta_decoy"]
        assert xls[i]["alpha_decoy"] is not new_alpha
        new_beta = xls_2[i]["alpha_decoy"] and xls_2[i]["beta_decoy"]
        assert xls[i]["beta_decoy"] is not new_beta


def test13():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = reannotate_decoy_labels(
        pr, by_mapping={True: None, False: None, None: None}
    )
    with open(FASTA, "r", encoding="utf-8") as f:
        pr_2 = reannotate_decoy_labels(pr_none, by_decoy_fasta=f)
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is not csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is not csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        # because MS Annika does not report decoy labels for individual peptides at XL level
        # we have to do like this
        new_alpha = xls_2[i]["alpha_decoy"] and xls_2[i]["beta_decoy"]
        assert xls[i]["alpha_decoy"] is not new_alpha
        new_beta = xls_2[i]["alpha_decoy"] and xls_2[i]["beta_decoy"]
        assert xls[i]["beta_decoy"] is not new_beta


def test14():
    import copy
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_none = copy.deepcopy(pr)
    csms = list()
    xls = list()
    for csm in pr_none["crosslink-spectrum-matches"]:
        if csm["alpha_decoy"]:
            csm = csm.copy_with_update(
                {
                    "alpha_proteins": [
                        f"REV__{protein}" for protein in csm["alpha_proteins"]
                    ]
                }
            )
        if csm["beta_decoy"]:
            csm = csm.copy_with_update(
                {
                    "beta_proteins": [
                        f"REV__{protein}" for protein in csm["beta_proteins"]
                    ]
                }
            )
        csms.append(csm)
    for xl in pr_none["crosslinks"]:
        if xl["alpha_decoy"]:
            xl = xl.copy_with_update(
                {
                    "alpha_proteins": [
                        f"REV__{protein}" for protein in xl["alpha_proteins"]
                    ]
                }
            )
        if xl["beta_decoy"]:
            xl = xl.copy_with_update(
                {
                    "beta_proteins": [
                        f"REV__{protein}" for protein in xl["beta_proteins"]
                    ]
                }
            )
        xls.append(xl)
    pr_none = pr_none.copy_with_update(
        {"crosslink-spectrum-matches": csms, "crosslinks": xls}
    )
    pr_none = reannotate_decoy_labels(
        pr_none, by_mapping={True: None, False: None, None: None}
    )

    def test_annotate(item: dict[str, any]) -> tuple[bool, bool]:
        alpha = False
        beta = False
        if all([protein.startswith("REV__") for protein in item["alpha_proteins"]]):
            alpha = True
        if all([protein.startswith("REV__") for protein in item["beta_proteins"]]):
            beta = True
        return alpha, beta

    pr_2 = reannotate_decoy_labels(pr_none, by_function=test_annotate)
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        assert xls[i]["alpha_decoy"] is xls_2[i]["alpha_decoy"]
        assert xls[i]["beta_decoy"] is xls_2[i]["beta_decoy"]


def test15():
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.transform import reannotate_decoy_labels

    xls = [
        create_crosslink_min(
            "ADANLDK", 7, "GNTDRHSIK", 9, proteins_a=None, proteins_b=["B"]
        )
    ]
    with pytest.warns(
        RuntimeWarning,
        match=r"Could not annotate alpha decoy label at index=0 because alpha proteins is 'None'!",
    ):
        _ = reannotate_decoy_labels(xls, by_decoy_protein_prefix="C")


def test16():
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.transform import reannotate_decoy_labels

    xls = [
        create_crosslink_min(
            "ADANLDK", 7, "GNTDRHSIK", 9, proteins_a=["A"], proteins_b=None
        )
    ]
    with pytest.warns(
        RuntimeWarning,
        match=r"Could not annotate beta decoy label at index=0 because beta proteins is 'None'!",
    ):
        _ = reannotate_decoy_labels(xls, by_decoy_protein_prefix="C")


def test17():
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.transform import reannotate_decoy_labels

    xls = [
        create_crosslink_min(
            "ADANLDK", 7, "GNTDRHSIK", 9, proteins_a=None, proteins_b=["B"]
        )
    ]
    with pytest.warns(
        RuntimeWarning,
        match=r"Could not annotate alpha decoy label at index=0 because alpha proteins is 'None'!",
    ):
        _ = reannotate_decoy_labels(xls, by_decoy_protein_substring="C")


def test18():
    from pyXLMS.data import create_crosslink_min
    from pyXLMS.transform import reannotate_decoy_labels

    xls = [
        create_crosslink_min(
            "ADANLDK", 7, "GNTDRHSIK", 9, proteins_a=["A"], proteins_b=None
        )
    ]
    with pytest.warns(
        RuntimeWarning,
        match=r"Could not annotate beta decoy label at index=0 because beta proteins is 'None'!",
    ):
        _ = reannotate_decoy_labels(xls, by_decoy_protein_substring="C")


def test19():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    xls = pr["crosslinks"]

    pr_2 = reannotate_decoy_labels(pr)
    csms_2 = pr_2["crosslink-spectrum-matches"]
    xls_2 = pr_2["crosslinks"]

    for i in range(len(csms)):
        assert csms[i]["alpha_decoy"] is csms_2[i]["alpha_decoy"]
        assert csms[i]["beta_decoy"] is csms_2[i]["beta_decoy"]

    for i in range(len(xls)):
        assert xls[i]["alpha_decoy"] is xls_2[i]["alpha_decoy"]
        assert xls[i]["beta_decoy"] is xls_2[i]["beta_decoy"]


def test20():
    from pyXLMS.parser import read
    from pyXLMS.transform import reannotate_decoy_labels

    pr = read([CSMS, XLS], engine="MS Annika", crosslinker="DSS")

    with pytest.raises(
        RuntimeError,
        match=r"Please only specify one option for reannotation, e.g. 'by_mapping' or 'by_target_fasta' but not both!",
    ):
        _ = reannotate_decoy_labels(
            pr, by_mapping={None: False}, by_decoy_protein_prefix="REV__"
        )
