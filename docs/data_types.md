# Data Types in pyXLMS

Most functions in pyXLMS return dictionaries with very specific structures, the
most notable are `crosslink dictionaries` storing crosslink information,
`crosslink-spectrum-match dictionaries` storing crosslink-spectrum-match
information, and `parser result dictionaries` storing crosslinks and/or
crosslink-spectrum-matches that have been parsed from a file. Each of these
dictionaries is explained below.

## Crosslink-Spectrum-Matches

Any method that returns or expects a crosslink-spectrum-match refers to a
dictionary with the following keys and values, values that are not guaranteed
may be `None`:

| Key                                  | Value Type                            | Value Guaranteed | Value Description |
|:-------------------------------------|:--------------------------------------|:----------------:|:------------------|
| "data_type"                          | `Literal["crosslink-spectrum-match"]` | ✅ | Data type of the dictionary. |
| "completeness"                       | `Literal["full", "partial"]`          | ✅ | Full crosslink-spectrum-matches have all their values set, partial crosslink-spectrum-matches have one or more missing values. |
| "alpha_peptide"                      | `str`                                 | ✅ | The unmodified amino acid sequence of the alpha peptide. |
| "alpha_modifications"                | `Dict[int, Tuple[str, float]]`        | ❌ | The modifications of the alpha peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass. N-terminal modifications are denoted with position 0. C-terminal modifications are denoted with position `len(peptide) + 1`. If the peptide is not modified an empty dictionary is given. |
| "alpha_peptide_crosslink_position"   | `int`                                 | ✅ | The position of the crosslinker in the sequence of the alpha peptide (1-based). |
| "alpha_proteins"                     | `List[str]`                           | ❌ | The accessions of proteins that the alpha peptide is associated with. |
| "alpha_proteins_crosslink_positions" | `List[int]`                           | ❌ | Positions of the crosslink in the proteins of the alpha peptide (1-based). |
| "alpha_proteins_peptide_positions"   | `List[int]`                           | ❌ | Positions of the alpha peptide in the corresponding proteins (1-based). |
| "alpha_score"                        | `float`                               | ❌ | Identification score of the alpha peptide. |
| "alpha_decoy"                        | `bool`                                | ❌ | Whether the alpha peptide is from the decoy database or not. |
| "beta_peptide"                       | `str`                                 | ✅ | The unmodified amino acid sequence of the beta peptide. |
| "beta_modifications"                 | `Dict[int, Tuple[str, float]]`        | ❌ | The modifications of the beta peptide given as a dictionary that maps peptide position (1-based) to modification given as a tuple of modification name and modification delta mass. N-terminal modifications are denoted with position 0. C-terminal modifications are denoted with position `len(peptide) + 1`. If the peptide is not modified an empty dictionary is given. |
| "beta_peptide_crosslink_position"    | `int`                                 | ✅ | The position of the crosslinker in the sequence of the beta peptide (1-based). |
| "beta_proteins"                      | `List[str]`                           | ❌ | The accessions of proteins that the beta peptide is associated with. |
| "beta_proteins_crosslink_positions"  | `List[int]`                           | ❌ | Positions of the crosslink in the proteins of the beta peptide (1-based). |
| "beta_proteins_peptide_positions"    | `List[int]`                           | ❌ | Positions of the beta peptide in the corresponding proteins (1-based). |
| "beta_score"                         | `float`                               | ❌ | Identification score of the beta peptide. |
| "beta_decoy"                         | `bool`                                | ❌ | Whether the beta peptide is from the decoy database or not. |
| "crosslink_type"                     | `Literal["inter", "intra"]`           | ✅ | If it is an inter- or intra-crosslink. |
| "score"                              | `float`                               | ❌ | Score of the crosslink-spectrum-match. |
| "spectrum_file"                      | `str`                                 | ✅ | Name of the spectrum file the crosslink-spectrum-match was identified in. |
| "scan_nr"                            | `int`                                 | ✅ | The corresponding scan number of the crosslink-spectrum-match. |
| "charge"                             | `int`                                 | ❌ | The precursor charge of the corresponding mass spectrum of the crosslink-spectrum-match. |
| "retention_time"                     | `float`                               | ❌ | The retention time of the corresponding mass spectrum of the crosslink-spectrum-match in seconds. |
| "ion_mobility"                       | `float`                               | ❌ | The ion mobility or compensation voltage of the corresponding mass spectrum of the crosslink-spectrum-match. |
| "additional_information"             | `Dict[str, Any]`                      | ❌ | A dictionary with additional information associated with the crosslink-spectrum-match. |

Additional resources:
- [API Documentation of the crosslink-spectrum-match creator](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.html#pyXLMS.data.create_csm)

## Crosslinks

Any method that returns or expects a crosslink refers to a
dictionary with the following keys and values, values that are not guaranteed
may be `None`:

| Key                                  | Value Type                            | Value Guaranteed | Value Description |
|:-------------------------------------|:--------------------------------------|:----------------:|:------------------|
| "data_type"                          | `Literal["crosslink"]`                | ✅ | Data type of the dictionary. |
| "completeness"                       | `Literal["full", "partial"]`          | ✅ | Full crosslinks have all their values set, partial crosslinks have one or more missing values. |
| "alpha_peptide"                      | `str`                                 | ✅ | The unmodified amino acid sequence of the alpha peptide. |
| "alpha_peptide_crosslink_position"   | `int`                                 | ✅ | The position of the crosslinker in the sequence of the alpha peptide (1-based). |
| "alpha_proteins"                     | `List[str]`                           | ❌ | The accessions of proteins that the alpha peptide is associated with. |
| "alpha_proteins_crosslink_positions" | `List[int]`                           | ❌ | Positions of the crosslink in the proteins of the alpha peptide (1-based). |
| "alpha_decoy"                        | `bool`                                | ❌ | Whether the alpha peptide is from the decoy database or not. |
| "beta_peptide"                       | `str`                                 | ✅ | The unmodified amino acid sequence of the beta peptide. |
| "beta_peptide_crosslink_position"    | `int`                                 | ✅ | The position of the crosslinker in the sequence of the beta peptide (1-based). |
| "beta_proteins"                      | `List[str]`                           | ❌ | The accessions of proteins that the beta peptide is associated with. |
| "beta_proteins_crosslink_positions"  | `List[int]`                           | ❌ | Positions of the crosslink in the proteins of the beta peptide (1-based). |
| "beta_decoy"                         | `bool`                                | ❌ | Whether the beta peptide is from the decoy database or not. |
| "crosslink_type"                     | `Literal["inter", "intra"]`           | ✅ | If it is an inter- or intra-crosslink. |
| "score"                              | `float`                               | ❌ | Score of the crosslink. |
| "additional_information"             | `Dict[str, Any]`                      | ❌ | A dictionary with additional information associated with the crosslink. |

Additional resources:
- [API Documentation of the crosslink creator](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.html#pyXLMS.data.create_crosslink)

## Parser Results

Any method that returns or expects a parser result (or `parser_result`) refers to a
dictionary with the following keys and values, values that are not guaranteed
may be `None`:

| Key                                  | Value Type                            | Value Guaranteed | Value Description |
|:-------------------------------------|:--------------------------------------|:----------------:|:------------------|
| "data_type"                          | `Literal["parser_result"]`            | ✅ | Data type of the dictionary. |
| "completeness"                       | `Literal["full", "partial", "empty"]` | ✅ | Full parser results have all their values set, partial parser results have one missing value, empty parser results contain no crosslink-spectrum-matches or crosslinks. |
| "search_engine"                      | `str`                                 | ✅ | Name of the identifying crosslink search engine. |
| "crosslink-spectrum-matches"         | `List[Dict[str, Any]]`                | ❌ | List of parsed crosslink-spectrum-matches. |
| "crosslinks"                         | `List[Dict[str, Any]]`                | ❌ | List of parsed crosslinks. |

Additional resources:
- [API Documentation of the parser result creator](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.html#pyXLMS.data.create_parser_result)
