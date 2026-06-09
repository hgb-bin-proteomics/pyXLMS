# Contributing to pyXLMS

Contributions to pyXLMS are very welcome! Contributions do not need to be in the
form of code, you can also contribute by letting us know of any bugs, suggesting
a new feature, writing more documentation, or generally giving us any kind of
feedback!

## Issues and Discussions

Please use [issues](https://github.com/hgb-bin-proteomics/pyXLMS/issues) for any
kind of bug reports or feature requests and select a suitable issue template.
For questions either open up an issue or - preferably - start a
[discussion](https://github.com/hgb-bin-proteomics/pyXLMS/discussions). If you
are unsure, you can always do any of the two or also contact us directly via
[proteomics@fh-hagenberg.at](mailto:proteomics@fh-hagenberg.at).

## Contributing Code

We are happy if you are deciding to directly contribute to the pyXLMS codebase.
If you are looking for something to get started on, check if there is any issues
labelled with `good first issue` that might need attention. Please be aware that
there are a few guidelines for code contributions:

- pyXLMS is built using a functional approach rather than an object-oriented one
  for two reasons:
  - Performance is usually better.
  - All data remains easy to serialize.
- Preferably your code follows this functional-oriented style.
- Please read through the [user guide](https://hgb-bin-proteomics.github.io/pyXLMS-docs/),
  specifically the section: `Documentation` ➡️ `Working with pyXLMS` ➡️ `Important Concepts`.
- Code needs to be type-hinted and type-checked.
  - We use [ty](https://docs.astral.sh/ty/) for type-checking.
  - Type-checking is done with `ty check`.
    - If you modify package code you should run: `ty check --config-file ty.toml`.
    - If you modify GUI code you should run: `cd gui && ty check --config-file ty-gui.toml`.
  - `ignore` flags may be used when there is a clear issue with the type checker.
- Code needs to pass linting and adhere to the formatting style.
  - We use [ruff](https://docs.astral.sh/ruff/) for linting and formatting.
  - Linting is done with `ruff check`.
  - Formatting is done with `ruff format`.
- Every function and object needs to be sufficiently documented.
- Every function and object needs to have sufficient tests.
- Every public function and object needs to have examples.
- Changes that break backward compatibility should be avoided! Please discuss with a
  maintainer first when it is absolutely necessary!
- **Pull requests should be based on and target the 'develop' branch. Merging directly**
  **into 'master' is not allowed.**
- Please run [pytest](https://docs.pytest.org/en/stable/) before opening a pull request:
  ```bash
  pytest -c pytest.ini --runslow --runext tests/
  ```
- Please bump the version info in `pyproject.toml`, `src/pyXLMS/__init__.py`, and
  `sphinx/conf.py` according to [semantic versioning rules](https://semver.org/)
  before opening a pull request.
- Pull requests must pass all GitHub actions checks.
- If a new parser is implemented, make sure to update `parser.read()` and
  `pipelines.pipeline()`.

Feel free to open a pull request even if your code does not conform with these
guidelines, we will still have a look and then adapt ourselves. Thank you!

## Responsible Use of LLMs/Gen-AI

pyXLMS is currently completely free of AI generated code or content and we aim to keep
it this way! While we do not prohibit outside contributors from using generative AI for
any kind of task, please ask yourself if it is absolutely necessary considering both
environmental and ethical concerns. Please also note that when contributing AI generated
code, issues, discussions, etc. you are still responsible for its content! AI generated
contributions have to be labelled as such! We retain the right to reject any kind of
contribution to ensure the quality of pyXLMS.

_Does this affect me if I just use pyXLMS?_ No! You may use pyXLMS in any way you want,
in fact pyXLMS integrates well with coding assistants and agents. This policy only affects
contributions to pyXLMS. You can also fork the code at any time and do whatever you want
with it!

*****

**Thank you for using and contributing to pyXLMS!**
