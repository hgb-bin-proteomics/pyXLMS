Description of the pyXLMS file format [TODO]

## Crosslink-Spectrum-Matches

| Column Name                        | Required | Data Type | Example 1 | Example 2     | Description |
|:-----------------------------------|:--------:|:----------|:----------|:--------------|:------------|
| Alpha Peptide                      | ✅       | str       | PEPKTIDE  | KPEPMTIDE     | Unmodified amino acid sequence of the alpha peptide in uppercase letters |
| Alpha Peptide Modifications        | ❌       | str       | (4:[DSS|138.06808]) | (1:[DSS|138.06808]);(5:[Oxidation|15.994915]) | See [Modification Encoding](#Modification_Encoding) |
| Alpha Peptide Crosslink Position   | ✅       | int       | 4         | 1             | Position of the crosslinker in the alpha peptide (1-based) |
| Alpha Proteins                     | ❌       | str       | G3ECR1    | G3ECR1;J7RUA5 | Accession of the associated protein(s) of the alpha peptide, if multiple proteins are given they should be delimited by a semicolon |
| Alpha Proteins Crosslink Positions | ❌       | int, str  | 13        | 13;15         | Position of the crosslinker in the associated alpha proteins, positions in multiple proteins should be delimited by a semicolon (1-based) |
| Alpha Decoy                        | ❌       | bool, str | False     | True          | Whether the alpha peptide is from the target (False) or decoy (True) database |
| Beta Peptide                       | ✅       | str       | PEPKTIDE  | KPEPTIDE      | Unmodified amino acid sequence of the beta peptide in uppercase letters |
| Beta Peptide Crosslink Position    | ✅       | int       | 4         | 1             | Position of the crosslinker in the beta peptide (1-based) |
| Beta Proteins                      | ❌       | str       | G3ECR1    | G3ECR1;J7RUA5 | Accession of the associated protein(s) of the beta peptide, if multiple proteins are given they should be delimited by a semicolon |
| Beta Proteins Crosslink Positions  | ❌       | int, str  | 13        | 13;15         | Position of the crosslinker in the associated beta proteins, positions in multiple proteins should be delimited by a semicolon (1-based) |
| Beta Decoy                         | ❌       | bool, str | False     | True          | Whether the beta peptide is from the target (False) or decoy (True) database |
| Score                              | ❌       | float     | 0.99513   | 170.3         | Score of the crosslink |

## Modification Encoding

## Crosslinks

| Column Name                        | Required | Data Type | Example 1 | Example 2     | Description |
|:-----------------------------------|:--------:|:----------|:----------|:--------------|:------------|
| Alpha Peptide                      | ✅       | str       | PEPKTIDE  | KPEPTIDE      | Unmodified amino acid sequence of the alpha peptide in uppercase letters |
| Alpha Peptide Crosslink Position   | ✅       | int       | 4         | 1             | Position of the crosslinker in the alpha peptide (1-based) |
| Alpha Proteins                     | ❌       | str       | G3ECR1    | G3ECR1;J7RUA5 | Accession of the associated protein(s) of the alpha peptide, if multiple proteins are given they should be delimited by a semicolon |
| Alpha Proteins Crosslink Positions | ❌       | int, str  | 13        | 13;15         | Position of the crosslinker in the associated alpha proteins, positions in multiple proteins should be delimited by a semicolon (1-based) |
| Alpha Decoy                        | ❌       | bool, str | False     | True          | Whether the alpha peptide is from the target (False) or decoy (True) database |
| Beta Peptide                       | ✅       | str       | PEPKTIDE  | KPEPTIDE      | Unmodified amino acid sequence of the beta peptide in uppercase letters |
| Beta Peptide Crosslink Position    | ✅       | int       | 4         | 1             | Position of the crosslinker in the beta peptide (1-based) |
| Beta Proteins                      | ❌       | str       | G3ECR1    | G3ECR1;J7RUA5 | Accession of the associated protein(s) of the beta peptide, if multiple proteins are given they should be delimited by a semicolon |
| Beta Proteins Crosslink Positions  | ❌       | int, str  | 13        | 13;15         | Position of the crosslinker in the associated beta proteins, positions in multiple proteins should be delimited by a semicolon (1-based) |
| Beta Decoy                         | ❌       | bool, str | False     | True          | Whether the beta peptide is from the target (False) or decoy (True) database |
| Score                              | ❌       | float     | 0.99513   | 170.3         | Score of the crosslink |
