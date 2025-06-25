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
  - Performance is usually better
  - All data remains easy to serialize
  Preferably your code follows this style.
- Code needs to be type-hinted and type-checked.
  - Type-checking is done with `pyright`.
  - `ignore` flags may be used when there is a clear issue with the type checker.
- Code needs to pass linting and adhere to the formatting style.
  - Linting is done with `ruff check`.
  - Formatting is done with `ruff format`.
- Every function and object needs to be sufficiently documented.
- Every function and object needs to have sufficient tests.
- Every public function and object needs to have examples.
- Pull requests should be opened for branch 'develop' first. Merging directly
  into 'master' is not allowed.
- Pull requests must pass all GitHub actions checks.
- If a new parser is implemented, make sure to update `parser.read()` and
  `pipelines.pipeline()`.

Feel free to open a pull request even if your code does not conform with these
guidelines, we will still have a look and then adapt ourselves. Thank you!

*****

Thank you for using and contributing to pyXLMS!
