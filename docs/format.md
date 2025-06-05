Description of the pyXLMS file format [TODO]

Crosslink-Spectrum-Matches

| Column Name                        | Required | Data Type | Description | Example 1 | Example 2     |
|:-----------------------------------|:--------:|:----------|:------------|:----------|:--------------|
| Alpha Peptide                      | ✅       | str       | desc        | PEPKTIDE  | KPEPTIDE      |
| Alpha Peptide Crosslink Position   | ✅       | int       | desc        | 4         | 1             |
| Alpha Proteins                     | ❌       | str       | desc        | G3ECR1    | G3ECR1;J7RUA5 |
| Alpha Proteins Crosslink Positions | ❌       | int, str  | desc        | 13        | 13;15         |
| Alpha Decoy                        | ❌       | bool, str | desc        | False     | True          |
| Beta Peptide                       | ✅       | str       | desc        | PEPKTIDE  | KPEPTIDE      |
| Beta Peptide Crosslink Position    | ✅       | int       | desc        | 4         | 1             |
| Beta Proteins                      | ❌       | str       | desc        | G3ECR1    | G3ECR1;J7RUA5 |
| Beta Proteins Crosslink Positions  | ❌       | int, str  | desc        | 13        | 13;15         |
| Beta Decoy                         | ❌       | bool, str | desc        | False     | True          |
| Score                              | ❌       | float     | desc        | 0.99513   | 170.3         |
