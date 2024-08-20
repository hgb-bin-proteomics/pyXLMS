![workflow_state](https://github.com/michabirklbauer/python_template/workflows/example/badge.svg)

# Template Repository for Python Scripts

A template repository for linting, testing, GUI building and dockerizing python scripts.

## Checklist

- Replace `YOURUSERNAME` and `IMAGENAME` in `.github/workflows/docker-image.yml` [or delete file].
- Replace `example` [line 5] with an appropriate name in `.github/workflows/python-app.yml`.
- Setup the appropriate copy statements [line >29] in `.github/workflows/python-app.yml`.
- Replace test data in `data` with your own data [or delete if you don't have test data].
- Adjust `.gitattributes` according to your needs.
- Adjust `.gitignore` according to your needs.
- Setup your `CITATION.cff` according to your needs [or delete file].
- Replace dummy values in `Dockerfile` and write image instructions.
- Replace copyright name in `LICENSE`.
- Replace lines 3 - 6 and write your script in `main.py`.
- Replace lines 3 - 6 and write your gui in `gui/streamlit_app.py`.
- Replace lines 3 - 6 and write tests in `tests/tests.py`.
- Add your requirements to `requirements.txt`.
- Document your script with [Sphinx](https://www.sphinx-doc.org/):
  - Install Sphinx and the [PyData](https://github.com/pydata/pydata-sphinx-theme) theme: `pip install sphinx pydata-sphinx-theme`.
  - Adjust the configuration to your needs in `docs_src/conf.py`.
  - Write documentation!
  - Build documentation with:
    ```
    sphinx-apidoc -f -o docs_src .
    sphinx-build -b html docs_src docs
    ```
  - Publish documentation [optional]!
    - Serving with GitHub pages needs the addition of an empty `.nojekyll` file to your `/docs`.
- Adjust this `README.md` to your needs!

## Known Issues

[List of known issues](https://github.com/michabirklbauer/python_template/issues)

## Citing

If you are using PLACEHOLDER script please cite:
```
Very important title
Important Author, and Another Important Author
Journal of Cool Stuff 2023 12 (3), 4567-4589
DOI: 12.3456/cool-stuff
```

## License

- [MIT](https://github.com/michabirklbauer/python_template/blob/master/LICENSE)

## Contact

- [micha.birklbauer@gmail.com](mailto:micha.birklbauer@gmail.com)
