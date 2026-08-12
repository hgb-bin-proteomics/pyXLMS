# pyXLMS Agents Instructions

pyXLMS is a python package and web application with graphical user interface that
aims to simplify and streamline the intermediate step of connecting crosslink search
engine results with down-stream analysis tools, enabling researchers even without
bioinformatics knowledge to conduct in-depth crosslink analyses and shifting the
focus from data transformation to data interpretation and therefore gaining biological
insight.

## Contributing Code

- pyXLMS is built using a functional approach rather than an object-oriented one,
  preferably your code follows this functional-oriented style.
- Code needs to be type-hinted and type-checked.
- We use [ty](https://docs.astral.sh/ty/) for type-checking.
- Type-checking is done with `ty check`.
- If you modify package code you should run: `ty check --config-file ty.toml`.
- If you modify GUI code you should run: `cd gui && ty check --config-file ty-gui.toml`.
- `ignore` flags may be used when there is a clear issue with the type checker.
- Code needs to pass linting and adhere to the formatting style.
- We use [ruff](https://docs.astral.sh/ruff/) for linting and formatting.
- Linting is done with `ruff check --config ruff.toml`.
- Formatting is done with `ruff format`.
- Every function and object needs to be sufficiently documented.
- Every function and object needs to have sufficient tests.
- Every public function and object needs to have examples.
- Changes that break backward compatibility should be avoided! Please discuss with a
  maintainer first when it is absolutely necessary!

## Branching

- Please base your work on the 'develop' branch!
- **Pull requests should be based on and target the 'develop' branch. Merging directly**
  **into 'master' is not allowed.**

## Linting

- We use [ruff](https://docs.astral.sh/ruff/) for linting and formatting.
- Make sure to use an up-to-date version of `ruff`.
- Linting is done with `ruff check --config ruff.toml`.
- Make sure that `ruff check --config ruff.toml` passes!

## Type Checking

- We use [ty](https://docs.astral.sh/ty/) for type-checking.
- Make sure to use an up-to-date version of `ty`.
- Type-checking is done with `ty check`.
- If you modify package code you should run: `ty check --config-file ty.toml`.
- If you modify GUI code you should run: `cd gui && ty check --config-file ty-gui.toml`.
- `ignore` flags may be used when there is a clear issue with the type checker.
- Make sure that `ty check --config-file ty.toml` passes!
- Make sure that `cd gui && ty check --config-file ty-gui.toml` passes!

## Formatting

- We use [ruff](https://docs.astral.sh/ruff/) for linting and formatting.
- Make sure to use an up-to-date version of `ruff`.
- Formatting is done with `ruff format`.
- Make sure that `ruff format --check` passes!

## Testing

- Please run [pytest](https://docs.pytest.org/en/stable/) before opening a pull request:
  ```bash
  pytest -c pytest.ini --runslow --runext tests/
  ```

## Important Concepts

- Please bump the version info in `pyproject.toml`, `src/pyXLMS/__init__.py`, and
  `sphinx/conf.py` according to [semantic versioning rules](https://semver.org/)
  before opening a pull request.
- Pull requests must pass all GitHub Actions checks.
- GitHub Actions workflows can be found in `.github/workflows`.
- If a new parser is implemented, make sure to update `pyXLMS.parser.read()` and
  `pyXLMS.pipelines.pipeline()`.

### Fail Fast

- Functions should be designed to fail fast meaning that inputs should be excessively
  checked and fail early if either a wrong data type is provided for any of the
  parameters or if the input crosslink-spectrum-matches, crosslinks, or parser results
  do not contain sufficient information to successfully run the function. For example,
  grouping by residue pairs requires that all crosslink-spectrum-matches have associated
  proteins and protein crosslink positions. The function will check this preemptively
  and throw an exception before performing the grouping if any of the information
  is missing!
- Use `pyXLMS.data.check_input()` and `pyXLMS.data.check_input_multi()` for function input checks!
- Use `pyXLMS.transform.assert_csms()` for functions that use crosslink-spectrum-matches as input!
- Use `pyXLMS.transform.assert_xls()` for functions that use crosslinks as input!
- Use `pyXLMS.transform.assert_csms_or_xls()` for functions that use crosslink-spectrum-matches or crosslinks as input!
- Use `pyXLMS.transform.check_available_keys()` to check if input data has all required attributes, e.g. attributes are not `None`!
- Use `pyXLMS.parser.format_sequence()` for reading amino acid sequences!
- Use `pXLMS.parser.get_bool_from_value()` for converting any value to a boolean!

### Peptide Order

- Upon creation of crosslink-spectrum-matches and crosslinks pyXLMS will re-order
  alpha and beta based on the peptide sequences and peptide crosslink positions.
  The peptide sequence and peptide crosslink position are fused for each peptide
  and the alpha (peptide) will always be the one where that fusion is alphabetically
  first, e.g. even if `peptide_a="TIDE"` and `peptide_b="PEP"` is specified, the
  resulting crosslink-spectrum-match or crosslink will have `"alpha_peptide"="PEP"`
  and `"beta_peptide"="TIDE"`. This ensures consistency and allows easy filtering
  of redundant/non-unique crosslinks.
- Ordering has to be carefully considered when retrieving data via `"additional_information"`:
  data associated with the alpha peptide in the additional information might map
  to the beta peptide instead due to the re-ordering! This needs to be manually
  checked via the peptide sequences!

### Modifying Data

- Attributes of crosslink-spectrum-matches, crosslinks, and parser results are frozen.
- To modify a crosslink-spectrum-match, crosslink, or parser result use its `copy_with_update()`
  method!

## Additional Information

- A user guide that documents all available functionality is available via [hgb-bin-proteomics.github.io/pyXLMS-docs](https://hgb-bin-proteomics.github.io/pyXLMS-docs).
- Example jupyter notebooks can be found in `/examples`.
- A full documentation of the python package can be accessed via [hgb-bin-proteomics.github.io/pyXLMS](https://hgb-bin-proteomics.github.io/pyXLMS).

## LLM/Gen-AI Policy

- Please label LLM/AI generated content as such, e.g. via commit attribution to
  the generating model!
