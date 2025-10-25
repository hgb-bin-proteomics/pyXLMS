# Supported pyXLMS I/O

The following section gives an overview of the supported crosslink search engines
and input file formats for reading crosslink results, as well as the supported
export options.

The following abbreviations are used:
- **CSMs:** crosslink-spectrum-matches
- **XLs:** crosslinks/residue pairs

## Supported Search Engines and Formats

The following crosslink search engines and formats are currently supported:

| Search Engine/Format | Parser            | CSMs | XLs | Supported Formats                           | API Docs |
|:---------------------|:------------------|:----:|:---:|:--------------------------------------------|:---------|
| MaxQuant/MaxLynx     | `read_maxquant()` | ✅  | ❌  | `.txt`, any text-like                       | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_maxquant.read_maxquant) |
| MeroX                | `read_merox()`    | ✅  | ❌  | `.csv`, `.zhrm`, any text-like              | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_merox.read_merox) |
| MS Annika            | `read_msannika()` | ✅  | ✅  | `.pdResult`, `.xlsx`, `.txt`, any text-like | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_msannika.read_msannika) |
| pLink 2 & 3          | `read_plink()`    | ✅  | ✅  | `.csv`, any text-like                       | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_plink.read_plink) |
| Scout                | `read_scout()`    | ✅  | ✅  | `.csv`, any text-like                       | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_scout.read_scout) |
| xiSearch & xiFDR     | `read_xi()`       | ✅  | ✅  | `.csv`, any text-like                       | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_xi.read_xi) |
| XlinkX               | `read_xlinkx()`   | ✅  | ✅  | `.pdResult`, `.xlsx`, `.txt`, any text-like | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_xlinkx.read_xlinkx) |
| mzIdentML            | `read_mzid()`     | ✅  | ❌  | `.mzid`                                     | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_mzid.read_mzid) |
| pyXLMS/custom        | `read_custom()`   | ✅  | ✅  | `.xlsx`, `.parquet`, `.csv`, any text-like  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.parser.html#pyXLMS.parser.parser_xldbse_custom.read_custom) |

_\*any text-like denotes all `.csv` like file formats that use a character delimiter e.g. `.csv`, `.tsv` and tab-delimited `.txt`._

## Supported Export Formats

The following export options are currently available:

| Tool/Format                             | Exporter                      | CSMs | XLs | API Docs | Tool Link |
|:----------------------------------------|:------------------------------|:----:|:---:|:---------|:----------|
| pandas.DataFrame/`.csv`/`.xlsx`/`.json` | `transform.to_dataframe()`    | ✅  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.transform.html#pyXLMS.transform.to_dataframe.to_dataframe) |   |
| AlphaLink2                              | `exporter.to_alphalink2()`    | ❌  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_alphalink2.to_alphalink2) | [link](https://github.com/Rappsilber-Laboratory/AlphaLink2) |
| IMP-X-FDR                               | `exporter.to_impxfdr()`       | ✅  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_impxfdr.to_impxfdr) | [link](https://github.com/vbc-proteomics-org/imp-x-fdr) |
| MS Annika                               | `exporter.to_msannika()`      | ✅  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_msannika.to_msannika) | [link](https://github.com/hgb-bin-proteomics/MSAnnika) |
| ProXL                                   | `exporter.to_proxl()`         | ✅  | ❌  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_proxl.to_proxl) | [link](https://www.yeastrc.org/proxl_public/) |
| PyXlinkViewer (pyMOL)                   | `exporter.to_pyxlinkviewer()` | ❌  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_pyxlinkviewer.to_pyxlinkviewer) | [link](https://github.com/BobSchiffrin/PyXlinkViewer) |
| xiFDR                                   | `exporter.to_xifdr()`         | ✅  | ❌  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_xifdr.to_xifdr) | [link](https://www.rappsilberlab.org/software/xifdr/) |
| xiNET                                   | `exporter.to_xinet()`         | ❌  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_xinet.to_xinet) | [link](https://crosslinkviewer.org/index.php) |
| xiVIEW                                  | `exporter.to_xiview()`        | ❌  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_xiview.to_xiview) | [link](https://www.xiview.org/index.php) |
| XlinkDB                                 | `exporter.to_xlinkdb()`       | ❌  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_xlinkdb.to_xlinkdb) | [link](https://xlinkdb.gs.washington.edu/xlinkdb/) |
| xlms-tools                              | `exporter.to_xlmstools()`     | ❌  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_xlmstools.to_xlmstools) | [link](https://gitlab.com/topf-lab/xlms-tools) |
| XMAS (ChimeraX)                         | `exporter.to_xmas()`          | ❌  | ✅  | [link](https://hgb-bin-proteomics.github.io/pyXLMS/pyXLMS.exporter.html#pyXLMS.exporter.to_xmas.to_xmas) | [link](https://github.com/ScheltemaLab/ChimeraX_XMAS_bundle) |
