import pandas as pd

def load_groups(filename: str) -> dict[str, dict[str, set[str]]]:
    df = pd.read_excel(filename, header=1)
    df.columns = [col.strip() for col in df.columns]

    libraries = ["group annotation in main library",
                 "group annotation in enrichable library",
                 "group annotation in acidic library"]
    
    result = {lib: {} for lib in libraries}

    for _, row in df.iterrows():
        peptide = row["tryptic sequence (miscleavage at XL site)"]
        for lib in libraries:
            group_cell = row[lib]
            if pd.isna(group_cell) or str(group_cell).lower() == "no":
                result[lib][peptide] = set() 
            else:
                groups = {g.strip() for g in str(group_cell).split(";")}
                result[lib][peptide] = groups
    
    return result


def is_correct_crosslink(peptide1: str, peptide2: str, peptide_to_groups: dict[str, set[str]]) -> bool:
    return not peptide_to_groups.get(peptide1, set()).isdisjoint(
        peptide_to_groups.get(peptide2, set())
    )

def match_peptides_multi(input_file: str, peptide_groups_by_library: dict[str, dict[str, set[str]]]) -> pd.DataFrame:
    sequences_df = pd.read_excel(input_file, usecols=["Sequence A", "Sequence B"])
    sequences_df = sequences_df.astype(str).replace(r"[\[\]]", "", regex=True)

    for lib_name, peptide_to_groups in peptide_groups_by_library.items():
        sequences_df[f"Same Group ({lib_name})"] = sequences_df.apply(
            lambda row: is_correct_crosslink(row["Sequence A"], row["Sequence B"], peptide_to_groups),
            axis=1
        )

    return sequences_df


if __name__ == "__main__":
    peptide_xlsx_file = "41467_2022_31701_MOESM4_ESM.xlsx"
    input_file = "R1_standard_deiso.xlsx"

    peptide_groups = load_groups(peptide_xlsx_file)
    result = match_peptides_multi(input_file, peptide_groups)

    result.to_csv("out2.csv", sep="\t", index=False)
    print(result)

    # True/False Counts pro Library
    for lib_name in peptide_groups.keys():
        counts = result[f"Same Group ({lib_name})"].value_counts()
        print(f"\nLibrary: {lib_name}")
        print(f"True: {counts.get(True, 0)}")
        print(f"False: {counts.get(False, 0)}")



