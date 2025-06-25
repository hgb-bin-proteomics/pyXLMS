# Description of the pyXLMS file format

Reading files with `parser.read(*, engine="Custom")` requires the following data format for crosslinks and crosslink-spectrum-matches. While the column names can be adjusted via the parameter
`column_mapping`, the format of the columns needs to stay the same for successful parsing. Any column that is not required can be safely omitted. This format is also output by
`transform.to_dataframe()`.

## Crosslink-Spectrum-Matches

Data required for parsing crosslink-spectrum-matches:

| Column Name                        | Required | Data Type | Example 1 | Example 2     | Description |
|:-----------------------------------|:--------:|:----------|:----------|:--------------|:------------|
| Alpha Peptide                      | ✅       | str       | PEPKTIDE  | KPEPMTIDE     | Unmodified amino acid sequence of the alpha peptide in uppercase letters |
| Alpha Peptide Modifications        | ❌       | str       | (4:[DSS\|138.06808]) | (1:[DSS\|138.06808]);(5:[Oxidation\|15.994915]) | Modifications of the alpha peptide, see ➡️ [Modification Encoding](#modification-encoding) |
| Alpha Peptide Crosslink Position   | ✅       | int       | 4         | 1             | Position of the crosslinker in the alpha peptide (1-based) |
| Alpha Proteins                     | ❌       | str       | G3ECR1    | G3ECR1;J7RUA5 | Accession of the associated protein(s) of the alpha peptide, if multiple proteins are given they should be delimited by a semicolon |
| Alpha Proteins Crosslink Positions | ❌       | int, str  | 13        | 13;15         | Position of the crosslinker in the associated alpha protein(s), positions in multiple proteins should be delimited by a semicolon (1-based) |
| Alpha Proteins Peptide Positions   | ❌       | int, str  | 10        | 13;15         | Position of the alpha peptide in the associated alpha protein(s), positions in multiple proteins should be delimited by a semicolon (1-based) |
| Alpha Score                        | ❌       | float     | 0.837     | 45.73         | Score of the alpha peptide |
| Alpha Decoy                        | ❌       | bool, str | False     | True          | Whether the alpha peptide is from the target (False) or decoy (True) database |
| Beta Peptide                       | ✅       | str       | PEPKTIDE  | KPEPTIDE      | Unmodified amino acid sequence of the beta peptide in uppercase letters |
| Beta Peptide Modifications         | ❌       | str       | (4:[DSS\|138.06808]) | (1:[DSS\|138.06808]);(5:[Oxidation\|15.994915]) | Modifications of the beta peptide, see ➡️ [Modification Encoding](#modification-encoding) |
| Beta Peptide Crosslink Position    | ✅       | int       | 4         | 1             | Position of the crosslinker in the beta peptide (1-based) |
| Beta Proteins                      | ❌       | str       | G3ECR1    | G3ECR1;J7RUA5 | Accession of the associated protein(s) of the beta peptide, if multiple proteins are given they should be delimited by a semicolon |
| Beta Proteins Crosslink Positions  | ❌       | int, str  | 13        | 13;15         | Position of the crosslinker in the associated beta protein(s), positions in multiple proteins should be delimited by a semicolon (1-based) |
| Beta Proteins Peptide Positions    | ❌       | int, str  | 10        | 13;15         | Position of the beta peptide in the associated beta protein(s), positions in multiple proteins should be delimited by a semicolon (1-based) |
| Beta Score                         | ❌       | float     | 0.837     | 45.73         | Score of the beta peptide |
| Beta Decoy                         | ❌       | bool, str | False     | True          | Whether the beta peptide is from the target (False) or decoy (True) database |
| CSM Score                          | ❌       | float     | 0.99513   | 170.3         | Score of the crosslink-spectrum-match |
| Spectrum File                      | ✅       | str       | 2025_03_17_EXP1_RUN3_R1.raw | 2025_03_17_EXP1_RUN3_R1 | File name of the spectrum file |
| Scan Nr                            | ✅       | int       | 1703      | 38901         | The scan number of the spectrum the match was identified in |
| Precursor Charge                   | ❌       | int       | 3         | 4             | Precursor charge of the crosslink spectrum |
| Retention Time                     | ❌       | float     | 530.17    | 1701.7        | Retention time of the crosslink spectrum in seconds |
| Ion Mobility                       | ❌       | float     | 170.41    | -50.0         | Ion mobility, CCS, or compensation voltage of the crosslink spectrum |

Additional resources:
- [API Documentation of the parser]()
- [API Documentation of the crosslink-spectrum-match creator]()

## Modification Encoding

Modifications are encoded with the following values:

- _position_: The 1-based position of the modification in the peptide sequence
  - should be parse-able as `int` data type
- _name_: The name of the modification, for example `Oxidation`
  - should be parse-able as `str` data type
- _mass_: The monoisotopic delta mass of the modification, for example `15.994915`
  - should be parse-able as `float` data type

Any modification is then encoded as `(position:[name|mass])`, multiple modifications should be delimited by a semicolon `;`. In the rare case that there is more than one modification on the
same position, their names should be delimited by a comma `,`. See examples below:

- `(4:[DSS|138.06808])`
- `(1:[DSS|138.06808]);(5:[Oxidation|15.994915])`
- `(5:[Substitution, Oxidation|13.541798])`

## Crosslinks

Data required for parsing crosslinks:

| Column Name                        | Required | Data Type | Example 1 | Example 2     | Description |
|:-----------------------------------|:--------:|:----------|:----------|:--------------|:------------|
| Alpha Peptide                      | ✅       | str       | PEPKTIDE  | KPEPTIDE      | Unmodified amino acid sequence of the alpha peptide in uppercase letters |
| Alpha Peptide Crosslink Position   | ✅       | int       | 4         | 1             | Position of the crosslinker in the alpha peptide (1-based) |
| Alpha Proteins                     | ❌       | str       | G3ECR1    | G3ECR1;J7RUA5 | Accession of the associated protein(s) of the alpha peptide, if multiple proteins are given they should be delimited by a semicolon |
| Alpha Proteins Crosslink Positions | ❌       | int, str  | 13        | 13;15         | Position of the crosslinker in the associated alpha protein(s), positions in multiple proteins should be delimited by a semicolon (1-based) |
| Alpha Decoy                        | ❌       | bool, str | False     | True          | Whether the alpha peptide is from the target (False) or decoy (True) database |
| Beta Peptide                       | ✅       | str       | PEPKTIDE  | KPEPTIDE      | Unmodified amino acid sequence of the beta peptide in uppercase letters |
| Beta Peptide Crosslink Position    | ✅       | int       | 4         | 1             | Position of the crosslinker in the beta peptide (1-based) |
| Beta Proteins                      | ❌       | str       | G3ECR1    | G3ECR1;J7RUA5 | Accession of the associated protein(s) of the beta peptide, if multiple proteins are given they should be delimited by a semicolon |
| Beta Proteins Crosslink Positions  | ❌       | int, str  | 13        | 13;15         | Position of the crosslinker in the associated beta protein(s), positions in multiple proteins should be delimited by a semicolon (1-based) |
| Beta Decoy                         | ❌       | bool, str | False     | True          | Whether the beta peptide is from the target (False) or decoy (True) database |
| Crosslink Score                    | ❌       | float     | 0.99513   | 170.3         | Score of the crosslink |

Additional resources:
- [API Documentation of the parser]()
- [API Documentation of the crosslink creator]()
