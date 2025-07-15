# pyXLMS
_a python package to process protein cross-linking data_

**pyXLMS** is a python package and web application with graphical user interface that aims to simplify and streamline the intermediate step of
connecting crosslink search engine results with down-stream analysis tools, enabling researchers even without bioinformatics knowledge to
conduct in-depth crosslink analyses and shifting the focus from data transformation to data interpretation and therefore gaining biological
insight. Currently pyXLMS supports input from six different crosslink search engines: [MaxLynx (part of MaxQuant)](https://www.maxquant.org/),
[MS Annika](https://github.com/hgb-bin-proteomics/MSAnnika),
[pLink 2 and pLink 3](http://pfind.ict.ac.cn/se/plink/),
[Scout](https://github.com/diogobor/Scout),
[xiSearch](https://www.rappsilberlab.org/software/xisearch/) and [xiFDR](https://www.rappsilberlab.org/software/xifdr/),
[XlinkX](https://docs.thermofisher.com/r/XlinkX-3.2-Quick-Start-Guide/),
as well as the [mzIdentML format](https://www.psidev.info/mzidentml)
of the HUPO Proteomics Standards Initiative, and a well-documented and
[human-readable custom tabular format](https://github.com/hgb-bin-proteomics/pyXLMS/blob/master/docs/format.md).
Down-stream analysis is facilitated by functionality that is directly available within pyXLMS such as validation, annotation, aggregation, and filtering of crosslink-spectrum-matches and crosslinks. In addition, the data can easily be exported to the required data format of the various available down-stream analysis tools such as
[xiNET](https://crosslinkviewer.org/index.php),
[xiVIEW](https://www.xiview.org/index.php),
[xiFDR](https://www.rappsilberlab.org/software/xifdr/),
[XlinkDB](https://xlinkdb.gs.washington.edu/xlinkdb/),
[xlms-tools](https://gitlab.com/topf-lab/xlms-tools),
pyMOL (via [pyXlinkViewer](https://github.com/BobSchiffrin/PyXlinkViewer)),
ChimeraX (via [XMAS](https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle)),
or [IMP-X-FDR](https://github.com/vbc-proteomics-org/imp-x-fdr).

## Installation

pyXLMS supports python version 3.7 and greater!

pyXLMS can easily be installed via pip:
```
pip install pyxlms
```

## Quick Start

After installation you can use pyXLMS in python like this:

_This example shows reading of MS Annika crosslink-spectrum-matches and exporting_
_them to xiFDR format for external validation._

```python
>>> import pyXLMS
>>> pr = pyXLMS.parser.read("data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx", engine="MS Annika", crosslinker="DSS")
Reading MS Annika CSMs...: 100%|████████████████████████████████████████████████████████████████████████████████| 826/826 [00:00<00:00, 20731.70it/s]
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
>>> _ = pyXLMS.exporter.to_xifdr(pr["crosslink-spectrum-matches"], filename="msannika_CSMs_for_xiFDR.csv")
```

## Web App

The web app is publicly accessible for free via [hgb-bin-proteomics.github.io/pyXLMS-app](https://hgb-bin-proteomics.github.io/pyXLMS-app).

Additionally, it can be run locally or self-hosted via:

- Running the [streamlit](https://streamlit.io/) app locally:
  - `cd gui`
  - `pip install -r requirements.txt`
  - `streamlit run streamlit_app.py`
- Building the [Docker](https://docs.docker.com/get-started/get-docker/) container yourself:
  - `docker build . -f Dockerfile -t pyxlms`
  - `docker run -p 8501:8501 pyxlms`
- _or_ by pulling the provided image from Docker Hub:
  - `docker run -p 8501:8501 michabirklbauer/pyxlms:latest`

## User Guide, Examples and Documentation

- A user guide that documents all available functionality is available via [hgb-bin-proteomics.github.io/pyXLMS-docs](https://hgb-bin-proteomics.github.io/pyXLMS-docs).
- Example jupyter notebooks can be found in `/examples`.
- A full documentation of the python package can be accessed via [hgb-bin-proteomics.github.io/pyXLMS](https://hgb-bin-proteomics.github.io/pyXLMS).

## Citing

If you are using pyXLMS please cite the following publication:

- Manuscript in preparation
  ```
  (wip)
  ```

## Contact

- [proteomics@fh-hagenberg.at](mailto:proteomics@fh-hagenberg.at)
- [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at) (primary developer)
