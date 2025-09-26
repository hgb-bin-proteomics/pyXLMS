import pandas as pd


def load_groups(filename: str) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    current_group: str = ""

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_group = line[1:].strip()
                groups[current_group] = []
            else:
                groups[current_group].append(line.rstrip())
    return groups

def invert_groups(groups: dict[str, list[str]]) -> dict[str, set[str]]:
    peptide_to_groups: dict[str, set[str]] = {}
    for group, peptides in groups.items():
        for peptide in peptides:
            peptide_to_groups.setdefault(peptide, set()).add(group)
    return peptide_to_groups

def is_correct_crosslink(peptide1: str, peptide2: str, peptide_to_groups: dict[str, set[str]]) -> bool:
    return not peptide_to_groups.get(peptide1, set()).isdisjoint(
        peptide_to_groups.get(peptide2, set())
    )

def match_peptides(input_file: str, peptide_to_groups: dict[str, set[str]]) -> pd.DataFrame:
    sequences_df: pd.DataFrame = pd.read_excel(input_file, usecols=["Sequence A", "Sequence B"])
    sequences_df = sequences_df.astype(str).replace(r"[\[\]]", "", regex=True)
    sequences_df["Same Group"] = sequences_df.apply(
        lambda row: is_correct_crosslink(row["Sequence A"], row["Sequence B"], peptide_to_groups),
        axis=1
    )
    return sequences_df

if __name__ == "__main__":
    peptide_groups_file = "peptide_groups.txt"
    input_file = "r3_standard.xlsx"

    groups: dict[str, list[str]] = load_groups(peptide_groups_file)
    peptide_to_groups: dict[str, set[str]] = invert_groups(groups)

    result: pd.DataFrame = match_peptides(input_file, peptide_to_groups)
    result.to_csv("out.csv", sep="\t", index=False)
    print(result)

    counts = result["Same Group"].value_counts()
    print(f"True: {counts.get(True, 0)}")
    print(f"False: {counts.get(False, 0)}")
