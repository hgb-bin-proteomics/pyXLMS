# pyXLMS

Supported search engines:
- MS Annika
- xiSearch / xiFDR
- MaxLynx
- ?

General interface with csv input

Packages to include:
- Export to xiNET
  - Req:
  - fasta file
  - Crosslink:
    - Sequence
    - XL position in peptide
    - protein accession
    - XL position in protein
    - Score
- Export to xiVIEW
  - Req: cover by xiNET Req
- Export to xiFDR
  - Req:
  - CSM:
    - Spectrum File
    - Scan Nr
    - Sequence
    - XL position in peptide
    - Precursor charge
    - Score CSM
    - Score peptide
    - protein accession
    - position of peptide in protein
    - decoy peptide
- Export to pyXlinkViewer
  - Req: all covered
- Export to XMAS
  - Req: all covered
- Export to Spectral Library
  - Req:
  - MGF
  - CSM:
    - Modification
    - RT
    - Ion Mobility / Compensation Voltage  
- MS Annika FDR
  - Req: all covered
- MS Annika Combine Results
  - Req: all covered
- CSM Annotation
  - Req: this is probably MS Annika only, as it requires doublet information
