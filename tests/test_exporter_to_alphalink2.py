#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

CAS9 = """
GAASMDKKYSIGLAIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIG
ALLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHRLE
ESFLVEEDKKHERHPIFGNIVDEVAYHEKYPTIYHLRKKLVDSTDKADLRLIYL
ALAHMIKFRGHFLIEGDLNPDNSDVDKLFIQLVQTYNQLFEENPINASGVDA
KAILSARLSKSRRLENLIAQLPGEKKNGLFGNLIALSLGLTPNFKSNFDLAE
DAKLQLSKDTYDDDLDNLLAQIGDQYADLFLAAKNLSDAILLSDILRVNTEITKA
PLSASMIKRYDEHHQDLTLLKALVRQQLPEKYKEIFFDQSKNGYAGYIDGG
ASQEEFYKFIKPILEKMDGTEELLVKLNREDLLRKQRTFDNGSIPHQIHLGE
LHAILRRQEDFYPFLKDNREKIEKILTFRIPYYVGPLARGNSRFAWMTRKSEETI
TPWNFEEVVDKGASAQSFIERMTNFDKNLPNEKVLPKHSLLYEYFTVYNELT
KVKYVTEGMRKPAFLSGEQKKAIVDLLFKTNRKVTVKQLKEDYFKKIECFD
SVEISGVEDRFNASLGTYHDLLKIIKDKDFLDNEENEDILEDIVLTLTLFEDREM
IEERLKTYAHLFDDKVMKQLKRRRYTGWGRLSRKLINGIRDKQSGKTILDFL
KSDGFANRNFMQLIHDDSLTFKEDIQKAQVSGQGDSLHEHIANLAGSPAIKK
GILQTVKVVDELVKVMGRHKPENIVIEMARENQTTQKGQKNSRERMKRIEEGIK
ELGSQILKEHPVENTQLQNEKLYLYYLQNGRDMYVDQELDINRLSDYDVDAI
VPQSFLKDDSIDNKVLTRSDKNRGKSDNVPSEEVVKKMKNYWRQLLNAKLIT
QRKFDNLTKAERGGLSELDKAGFIKRQLVETRQITKHVAQILDSRMNTKYDEND
KLIREVKVITLKSKLVSDFRKDFQFYKVREINNYHHAHDAYLNAVVGTALIK
KYPKLESEFVYGDYKVYDVRKMIAKSEQEIGKATAKYFFYSNIMNFFKTEI
TLANGEIRKRPLIETNGETGEIVWDKGRDFATVRKVLSMPQVNIVKKTEVQTGGF
SKESILPKRNSDKLIARKKDWDPKKYGGFDSPTVAYSVLVVAKVEKGKSKK
LKSVKELLGITIMERSSFEKNPIDFLEAKGYKEVKKDLIIKLPKYSLFEL
ENGRKRMLASAGELQKGNELALPSKYVNFLYLASHYEKLKGSPEDNEQKQLFVEQ
HKHYLDEIIEQISEFSKRVILADANLDKVLSAYNKHRDKPIREQAENIIHL
FTLTNLGAPAAFKYFDTTIDRKQYRSTKEVLDATLIHQSITGLYETRIDLS
QLGGD
""".replace("\n", "").replace(" ", "")


def test1():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import filter_proteins
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    cas9 = filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]
    export = to_alphalink2(
        cas9, fasta="data/_fasta/Cas9_plus10.fasta", filename_prefix="Cas9"
    )
    assert export["Exported files"][0] == "Cas9_AlphaLink2.txt"
    assert export["Exported files"][1] == "Cas9_AlphaLink2.fasta"


def test2():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import filter_proteins
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    cas9 = filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]
    export = to_alphalink2(cas9, fasta="data/_fasta/Cas9_plus10.fasta")
    assert len(
        [
            line
            for line in export["AlphaLink2 crosslinks"].split("\n")
            if line.strip() != ""
        ]
    ) == len(cas9)
    assert export["AlphaLink2 FASTA"].split("\n")[0].split(">")[1].split("|")[0] == "A"
    assert export["AlphaLink2 FASTA"].split("\n")[1].strip() == CAS9
    assert export["AlphaLink2 DataFrame"].shape == (223, 5)
    for i, row in export["AlphaLink2 DataFrame"].iterrows():
        assert (
            int(row["residueFrom"]) == cas9[i]["alpha_proteins_crosslink_positions"][0]
        )
        assert str(row["chain1"]) == "A"
        assert int(row["residueTo"]) == cas9[i]["beta_proteins_crosslink_positions"][0]
        assert str(row["chain2"]) == "A"
        assert float(row["FDR"]) == pytest.approx(0.01)
    assert len(export["Exported files"]) == 0


def test3():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import filter_proteins
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    cas9 = filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]
    export = to_alphalink2(
        cas9, fasta="data/_fasta/Cas9_plus10.fasta", annotated_fdr=0.02
    )
    assert len(
        [
            line
            for line in export["AlphaLink2 crosslinks"].split("\n")
            if line.strip() != ""
        ]
    ) == len(cas9)
    assert export["AlphaLink2 FASTA"].split("\n")[0].split(">")[1].split("|")[0] == "A"
    assert export["AlphaLink2 FASTA"].split("\n")[1].strip() == CAS9
    assert export["AlphaLink2 DataFrame"].shape == (223, 5)
    for i, row in export["AlphaLink2 DataFrame"].iterrows():
        assert (
            int(row["residueFrom"]) == cas9[i]["alpha_proteins_crosslink_positions"][0]
        )
        assert str(row["chain1"]) == "A"
        assert int(row["residueTo"]) == cas9[i]["beta_proteins_crosslink_positions"][0]
        assert str(row["chain2"]) == "A"
        assert float(row["FDR"]) == pytest.approx(0.02)
    assert len(export["Exported files"]) == 0


def test4():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.transform import filter_proteins
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    cas9 = filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]
    fdr = [r / len(cas9) for r in range(1, len(cas9) + 1)]
    export = to_alphalink2(
        cas9, fasta="data/_fasta/Cas9_plus10.fasta", annotated_fdr=fdr
    )
    assert len(
        [
            line
            for line in export["AlphaLink2 crosslinks"].split("\n")
            if line.strip() != ""
        ]
    ) == len(cas9)
    assert export["AlphaLink2 FASTA"].split("\n")[0].split(">")[1].split("|")[0] == "A"
    assert export["AlphaLink2 FASTA"].split("\n")[1].strip() == CAS9
    assert export["AlphaLink2 DataFrame"].shape == (223, 5)
    for i, row in export["AlphaLink2 DataFrame"].iterrows():
        assert (
            int(row["residueFrom"]) == cas9[i]["alpha_proteins_crosslink_positions"][0]
        )
        assert str(row["chain1"]) == "A"
        assert int(row["residueTo"]) == cas9[i]["beta_proteins_crosslink_positions"][0]
        assert str(row["chain2"]) == "A"
        assert float(row["FDR"]) == pytest.approx(fdr[i])
    assert len(export["Exported files"]) == 0


def test5():
    from pyXLMS.parser import read
    from pyXLMS import transform
    from pyXLMS.exporter import to_alphalink2

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    pr = transform.annotate_fdr(pr)
    pr = transform.unique(pr, by="protein")
    pr = transform.validate(pr)
    pr = transform.targets_only(pr)
    cas9 = transform.filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]

    fdr = [r / len(cas9) for r in range(1, len(cas9) + 1)]
    annotated_fdr = [
        xl["additional_information"]["pyXLMS_annotated_FDR"] for xl in cas9
    ]
    export = to_alphalink2(
        cas9, fasta="data/_fasta/Cas9_plus10.fasta", annotated_fdr=fdr
    )
    assert len(
        [
            line
            for line in export["AlphaLink2 crosslinks"].split("\n")
            if line.strip() != ""
        ]
    ) == len(cas9)
    assert export["AlphaLink2 FASTA"].split("\n")[0].split(">")[1].split("|")[0] == "A"
    assert export["AlphaLink2 FASTA"].split("\n")[1].strip() == CAS9
    assert export["AlphaLink2 DataFrame"].shape == (223, 5)
    for i, row in export["AlphaLink2 DataFrame"].iterrows():
        assert (
            int(row["residueFrom"]) == cas9[i]["alpha_proteins_crosslink_positions"][0]
        )
        assert str(row["chain1"]) == "A"
        assert int(row["residueTo"]) == cas9[i]["beta_proteins_crosslink_positions"][0]
        assert str(row["chain2"]) == "A"
        assert float(row["FDR"]) == pytest.approx(annotated_fdr[i])
    assert len(export["Exported files"]) == 0


def test6():
    from pyXLMS.parser import read
    from pyXLMS import transform
    from pyXLMS.exporter import to_alphalink2

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    pr = transform.annotate_fdr(pr)
    pr = transform.unique(pr, by="protein")
    pr = transform.validate(pr)
    pr = transform.targets_only(pr)
    cas9 = transform.filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]

    fdr = [r / len(cas9) for r in range(1, len(cas9) + 1)]
    export = to_alphalink2(
        cas9,
        fasta="data/_fasta/Cas9_plus10.fasta",
        annotated_fdr=fdr,
        try_use_annotated_fdr=False,
    )
    assert len(
        [
            line
            for line in export["AlphaLink2 crosslinks"].split("\n")
            if line.strip() != ""
        ]
    ) == len(cas9)
    assert export["AlphaLink2 FASTA"].split("\n")[0].split(">")[1].split("|")[0] == "A"
    assert export["AlphaLink2 FASTA"].split("\n")[1].strip() == CAS9
    assert export["AlphaLink2 DataFrame"].shape == (223, 5)
    for i, row in export["AlphaLink2 DataFrame"].iterrows():
        assert (
            int(row["residueFrom"]) == cas9[i]["alpha_proteins_crosslink_positions"][0]
        )
        assert str(row["chain1"]) == "A"
        assert int(row["residueTo"]) == cas9[i]["beta_proteins_crosslink_positions"][0]
        assert str(row["chain2"]) == "A"
        assert float(row["FDR"]) == pytest.approx(fdr[i])
    assert len(export["Exported files"]) == 0


def test7():
    from pyXLMS.pipelines import pipeline
    from pyXLMS import transform
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    cas9 = transform.filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]
    fdr = [r / len(cas9) for r in range(1, len(cas9) + 1)][: len(cas9) - 10]
    with pytest.raises(
        ValueError,
        match=(
            "Length of annotated_fdr does not match length of crosslinks! "
            "When providing a list it needs to contain FDR values for every crosslink and therefore be of equal length!"
        ),
    ):
        _export = to_alphalink2(
            cas9, fasta="data/_fasta/Cas9_plus10.fasta", annotated_fdr=fdr
        )


def test8():
    from pyXLMS.pipelines import pipeline
    from pyXLMS import transform
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    cas9 = transform.filter_proteins(pr["crosslinks"], proteins=["Cas9"])["Both"]
    with pytest.raises(
        IndexError,
        match="Found more than the supported 62 proteins/chains in the fasta file! Please trim fasta file to a maximum of 62 sequences!",
    ):
        _export = to_alphalink2(
            cas9,
            fasta="data/_fasta/Ribosome_for_XL_all_contained_Proteins_from_shotgun.fasta",
        )


def test9():
    from pyXLMS.exporter.to_alphalink2 import __get_proteins_and_positions

    proteins_and_positions = __get_proteins_and_positions(
        "GSQKDR", {"A": {"header": "Cas9", "sequence": CAS9}}
    )
    assert len(proteins_and_positions[0]) == 0
    assert len(proteins_and_positions[1]) == 0


def test10():
    from pyXLMS.exporter.to_alphalink2 import __get_proteins_and_positions

    with pytest.raises(RuntimeError, match="No match found for peptide GSQKDR!"):
        _proteins_and_positions = __get_proteins_and_positions(
            "GSQKDR", {"A": {"header": "Cas9", "sequence": CAS9}}, True
        )


def test11():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    export = to_alphalink2(pr["crosslinks"], fasta="data/_fasta/Cas9.fasta", verbose=0)
    assert export["AlphaLink2 DataFrame"].shape[0] == 223


def test12():
    from pyXLMS.parser import read
    from pyXLMS.transform import targets_only
    from pyXLMS.transform import filter_proteins
    from pyXLMS.exporter import to_alphalink2

    pr = read(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    should = len(
        filter_proteins(targets_only(pr["crosslinks"]), proteins=["Cas9"])["Both"]
    )
    assert len(pr["crosslinks"]) == 300
    export = to_alphalink2(pr["crosslinks"], fasta="data/_fasta/Cas9.fasta", verbose=0)
    assert export["AlphaLink2 DataFrame"].shape[0] == should


def test13():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    with pytest.warns(
        RuntimeWarning,
        match=r"Could not find matching proteins in FASTA file for crosslink id 148:SDNVPSEEVVKKMK-VASMASEKMK! This warning can be ignored if this is to be expected\.",
    ):
        _export = to_alphalink2(pr["crosslinks"], fasta="data/_fasta/Cas9.fasta")


def test14():
    from pyXLMS.pipelines import pipeline
    from pyXLMS.exporter import to_alphalink2

    pr = pipeline(
        "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult",
        engine="MS Annika",
        crosslinker="DSS",
    )
    with pytest.raises(
        RuntimeError,
        match=r"Could not find matching proteins in FASTA file for crosslink id 148:SDNVPSEEVVKKMK-VASMASEKMK! If this is to be expected please set verbose level to either 1 or 0!",
    ):
        _export = to_alphalink2(
            pr["crosslinks"], fasta="data/_fasta/Cas9.fasta", verbose=2
        )
