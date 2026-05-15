# pyXLMS-pydantic
_a python package to process protein cross-linking data_

<img src="https://github.com/hgb-bin-proteomics/pyXLMS/raw/master/docs/logo/logo_padded_shadow.png" class="dark-light" align="left" width="200px" style="padding: 5px 20px 10px 20px;"/>

**pyXLMS-pydantic** is the [Pydantic](https://pydantic.dev/) version of [**pyXLMS**](https://github.com/hgb-bin-proteomics/pyXLMS):
a python package that aims to simplify and streamline the intermediate step of
connecting crosslink search engine results with down-stream analysis tools, enabling researchers even without bioinformatics knowledge to
conduct in-depth crosslink analyses and shifting the focus from data transformation to data interpretation and therefore gaining biological
insight. Currently pyXLMS supports input from several different crosslink search engines including:
[MaxLynx (part of MaxQuant)](https://www.maxquant.org/),
[MeroX](https://www.stavrox.com/),
[MS Annika](https://github.com/hgb-bin-proteomics/MSAnnika),
[pLink 2 and pLink 3](http://pfind.ict.ac.cn/se/plink/),
[Scout](https://github.com/diogobor/Scout),
[xiSearch](https://www.rappsilberlab.org/software/xisearch/) and [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
[XlinkX](https://docs.thermofisher.com/r/XlinkX-3.2-Quick-Start-Guide/),
as well as the [mzIdentML format](https://www.psidev.info/mzidentml)
of the HUPO Proteomics Standards Initiative, and a well-documented and
[human-readable custom tabular format](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md).
Down-stream analysis is facilitated by functionality that is directly available within pyXLMS such as validation, annotation, aggregation, filtering, and visualization - and [much more](https://hgb-bin-proteomics.github.io/pyXLMS/modules.html) - of crosslink-spectrum-matches and crosslinks. In addition, the data can easily be exported to the required data format of the various available down-stream analysis tools such as
[AlphaLink2](https://github.com/Rappsilber-Laboratory/AlphaLink2),
[ProXL](https://www.yeastrc.org/proxl_public/),
[xiNET](https://crosslinkviewer.org/index.php),
[xiVIEW](https://www.xiview.org/index.php),
[xiFDR](https://www.rappsilberlab.org/software/xifdr/),
[XlinkDB](https://xlinkdb.gs.washington.edu/xlinkdb/),
[xlms-tools](https://gitlab.com/topf-lab/xlms-tools),
PyMOL (via [PyXlinkViewer](https://github.com/BobSchiffrin/PyXlinkViewer)),
ChimeraX (via [XMAS](https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle)),
or [IMP-X-FDR](https://github.com/vbc-proteomics-org/imp-x-fdr).

## Why pyXLMS-pydantic and what are the differences to the classic pyXLMS?

pyXLMS-pydantic is the [Pydantic](https://pydantic.dev/) version of the [classic pyXLMS](https://github.com/hgb-bin-proteomics/pyXLMS). Using Pydantic offers
increased type safety and validation, programmatically easier handling of crosslink-spectrum-matches and crosslinks. However, this comes at the slight cost of
performance and compatibility - pyXLMS-pydantic requires python version 3.14 or newer!

For writing code there should not be any differences between the Pydantic and classic version, e.g. all the code written for the classic pyXLMS should also work for
pyXLMS-pydantic. Crosslink-spectrum-matches, crosslinks, and parser results are now standalone classes based on the Pydantic
[BaseModel](https://pydantic.dev/docs/validation/latest/api/pydantic/base_model/)
instead of dictionaries though. Some private functions have also been moved. The terms pyXLMS and pyXLMS-pydantic are used synonymously throughout this repository, unless
specifically stated!

## Installation

**Installing pyXLMS-pydantic will override the [classic pyXLMS](https://github.com/hgb-bin-proteomics/pyXLMS) version if it was previously installed!** 

pyXLMS-pydantic supports python **version 3.14 and greater** and can easily be installed via pip:
```
pip install git+https://github.com/hgb-bin-proteomics/pyXLMS-pydantic.git
```

## Quick Start

After installation you can use pyXLMS in python like this:

_This example shows reading of MS Annika crosslink-spectrum-matches and exporting_
_them to xiFDR format for external validation._

```python
>>> import pyXLMS
>>> pr = pyXLMS.parser.read(
...     "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx",
...     engine="MS Annika",
...     crosslinker="DSS"
... )
Reading MS Annika CSMs...: 100%|████████████████| 826/826 [00:00<00:00, 20731.70it/s]
>>> _ = pyXLMS.transform.summary(pr)
Number of CSMs: 826.0
Number of unique CSMs: 826.0
Number of intra CSMs: 803.0
Number of inter CSMs: 23.0
Number of target-target CSMs: 786.0
Number of target-decoy CSMs: 39.0
Number of decoy-decoy CSMs: 1.0
Minimum CSM score: 1.11
Maximum CSM score: 452.99
>>> _ = pyXLMS.exporter.to_xifdr(
...     pr["crosslink-spectrum-matches"],
...     filename="msannika_CSMs_for_xiFDR.csv"
... )
```

## Examples and Documentation

- Example jupyter notebooks can be found in `/examples`.
- A full documentation of the python package can be accessed via [hgb-bin-proteomics.github.io/pyXLMS-pydantic](https://hgb-bin-proteomics.github.io/pyXLMS-pydantic).

## FAQ

Not sure if pyXLMS is what you are looking for? You can find a collection of common questions and answers about pyXLMS
[here](https://github.com/hgb-bin-proteomics/pyXLMS/discussions/categories/q-a) or in the
[user guide](https://hgb-bin-proteomics.github.io/pyXLMS-docs) under `Documentation` ➡️ `FAQ`.

## Limitations

Despite our best efforts pyXLMS still comes with some limitations that are a direct result of the differences in the output formats of crosslink search engines.
Many crosslink search engines do not report any kind of decoy matches which makes validation in pyXLMS or export to xiFDR impossible which is why **we recommend**
**using validated results for pyXLMS**. Validation within pyXLMS is currently only supported for MaxLynx, MS Annika, and xiSearch.

Furthermore, the different down-stream
analysis tools require varying input information which might not be consistently available from all crosslink search engines. Some of this can be mitigated by functionality
in pyXLMS such as annotation or by additional information that needs to be passed to pyXLMS for a successful export. Generally, the export to all downstream analysis tools
should work for all crosslink search engines and input formats, with the exception of the export to xiFDR which is limited to MaxLynx, MS Annika, and xiSearch for above reasons.
For safety pyXLMS makes sure before the export that all the required information is available and will otherwise throw an error. For more information please check the specific
export pages in the [user guide](https://hgb-bin-proteomics.github.io/pyXLMS-docs) and the [documentation](https://hgb-bin-proteomics.github.io/pyXLMS).

The web app supports most of the features of the python package, features that are not supported in the web app are listed in the [user guide](https://hgb-bin-proteomics.github.io/pyXLMS-docs)
under `Documentation` ➡️ `Web Application` ➡️ `Feature Support`.

Interacting with [STRING](https://string-db.org/) requires an active internet connection and depends on the service availability of STRING.

## Citing

If you are using pyXLMS please cite as described [here](https://github.com/hgb-bin-proteomics/pyXLMS).

## Acknowledgements

We thank Melanie Birklbauer for designing the logo.

## Contact

- [proteomics@fh-hagenberg.at](mailto:proteomics@fh-hagenberg.at)
- [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at) (primary developer)
