import pandas as pd

def load_acidic_library(filename: str) -> dict[str, set[str]]:
    df = pd.read_excel(filename, header=1)
    df.columns = [col.strip() for col in df.columns]

    peptide_to_groups = {}
    for _, row in df.iterrows():
        peptide = row["tryptic sequence (miscleavage at XL site)"]
        group_cell = row["group annotation in acidic library"]
        if pd.isna(group_cell) or str(group_cell).lower() == "no":
            peptide_to_groups[peptide] = set()
        else:
            peptide_to_groups[peptide] = {g.strip() for g in str(group_cell).split(";")}
    return peptide_to_groups

def is_correct_crosslink(peptide1: str, peptide2: str, peptide_to_groups: dict[str, set[str]]) -> bool:
    return not peptide_to_groups.get(peptide1, set()).isdisjoint(
        peptide_to_groups.get(peptide2, set())
    )

def match_acidic_library(input_file: str, peptide_to_groups: dict[str, set[str]]) -> pd.DataFrame:
    sequences_df = pd.read_excel(input_file, usecols=["Sequence A", "Sequence B"])
    sequences_df = sequences_df.astype(str).replace(r"[\[\]]", "", regex=True)
    sequences_df["Same Group (acidic library)"] = sequences_df.apply(
        lambda row: is_correct_crosslink(row["Sequence A"], row["Sequence B"], peptide_to_groups),
        axis=1
    )
    return sequences_df

if __name__ == "__main__":
    peptide_file = "41467_2022_31701_MOESM4_ESM.xlsx"
    input_file = "R1_standard_deiso.xlsx"

    peptide_groups = load_acidic_library(peptide_file)
    result = match_acidic_library(input_file, peptide_groups)

    result.to_csv("out_acidic_library.csv", sep="\t", index=False)
    print(result)

    counts = result["Same Group (acidic library)"].value_counts()
    print(f"True: {counts.get(True, 0)}")
    print(f"False: {counts.get(False, 0)}")
