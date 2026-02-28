# tcga-gdc-downloader

![Tests](https://github.com/onedimkurt/tcga-gdc-downloader/actions/workflows/tests.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Linux%20%7C%20Windows-lightgrey)

A command-line and graphical tool to download TCGA RNA-seq gene expression count data from the [NCI GDC portal](https://portal.gdc.cancer.gov), assemble a ready-to-analyse gene × sample count matrix, and automatically annotate every sample with GDC clinical data and harmonised survival and subtype annotations from the [PanCanAtlas Clinical Data Resource](https://www.cell.com/cell/fulltext/S0092-8674(18)30229-0) (Liu et al., *Cell* 2018).

No coding required. Works from a graphical interface or a single terminal command.

---

## Features

- Downloads all open-access STAR-Counts RNA-seq files for any of the 33 TCGA cancer projects
- Assembles a gene × sample count matrix (genes as rows, samples as columns)
- Fetches GDC clinical data (demographics, tumour stage, histology, survival)
- Merges PanCanAtlas CDR annotations: four curated survival endpoints (OS, DSS, DFI, PFI) and molecular subtypes
- Classifies every sample as Tumor or Normal from the TCGA barcode
- Produces separate output files for all samples, tumor-only, and normal-only
- Flags CDR-complete cases (full survival data) for clean analytical cohorts
- Checkpoint system — safely resumes interrupted downloads without re-downloading
- Graphical 9-step wizard (Streamlit) and command-line interface

---

## Installation

Requires Python 3.10 or newer.

```bash
pip install tcga-gdc-downloader
```

For the graphical interface, install the GUI extra:

```bash
pip install "tcga-gdc-downloader[gui]"
```

---

## Quick Start

### Command line

```bash
# Discover files without downloading (recommended first step)
tcga-download --project TCGA-BRCA --dry-run

# Full download and annotation
tcga-download --project TCGA-BRCA --output ~/my_data

# Skip CDR annotation step
tcga-download --project TCGA-BRCA --output ~/my_data --no-cdr

# Resume an interrupted download
tcga-download --project TCGA-BRCA --output ~/my_data

# Start completely fresh
tcga-download --project TCGA-BRCA --output ~/my_data --fresh
```

### Graphical interface

```bash
tcga-download --gui
```

A browser window opens with a 9-step wizard. No further terminal commands needed.

### Python API

```python
from tcga_downloader import GDCClient, build_count_matrix, run_cdr_pipeline

client = GDCClient()
hits = client.discover_star_files("TCGA-BRCA")
```

---

## Output Files

For a project called `TCGA-BRCA`, the following files are written to your output directory:

| File | Contents |
|------|----------|
| `TCGA-BRCA_STAR_unstranded_merged_ALL.tsv` | All samples — counts + GDC clinical |
| `TCGA-BRCA_STAR_unstranded_TUMOR_ONLY.tsv` | Tumor samples only |
| `TCGA-BRCA_STAR_unstranded_NORMAL_ONLY.tsv` | Normal samples only (if present) |
| `TCGA-BRCA_sample_metadata_clinical.tsv` | Metadata + clinical only, no counts |
| `TCGA-BRCA_FULL_merged_with_CDR.tsv` | All samples — counts + GDC clinical + CDR |
| `TCGA-BRCA_CDR_annotations.tsv` | CDR columns only |
| `TCGA-BRCA_CDR_coverage_report.tsv` | Per-field CDR coverage statistics |
| `TCGA-BRCA_FULL_merged_CDR_complete_cases.tsv` | Samples with complete survival data |
| `TCGA-BRCA_CDR_unmatched_cases.txt` | Cases not found in CDR (if any) |

### Opening your files

```python
# Python
import pandas as pd
df = pd.read_csv("TCGA-BRCA_FULL_merged_with_CDR.tsv", sep="\t", index_col=0)
```

```r
# R
df <- read.table("TCGA-BRCA_FULL_merged_with_CDR.tsv", sep="\t", header=TRUE, row.names=1)
```

---

## Supported Projects

All 33 TCGA cancer projects are supported:

`TCGA-ACC` `TCGA-BLCA` `TCGA-BRCA` `TCGA-CESC` `TCGA-CHOL` `TCGA-COAD`
`TCGA-DLBC` `TCGA-ESCA` `TCGA-GBM` `TCGA-HNSC` `TCGA-KICH` `TCGA-KIRC`
`TCGA-KIRP` `TCGA-LAML` `TCGA-LGG` `TCGA-LIHC` `TCGA-LUAD` `TCGA-LUSC`
`TCGA-MESO` `TCGA-OV` `TCGA-PAAD` `TCGA-PCPG` `TCGA-PRAD` `TCGA-READ`
`TCGA-SARC` `TCGA-SKCM` `TCGA-STAD` `TCGA-TGCT` `TCGA-THCA` `TCGA-THYM`
`TCGA-UCEC` `TCGA-UCS` `TCGA-UVM`

---

## GDC Authentication

TCGA STAR-Counts gene expression data is **open access** — no authentication token is required. A GDC token is only needed for controlled-access data (raw BAM files, genotype data), which this tool does not download.

To use a token if you have one:

```bash
tcga-download --project TCGA-BRCA --token ~/gdc-user-token.txt
```

---

## PanCanAtlas CDR Annotations

The tool automatically downloads and merges the TCGA Pan-Cancer Clinical Data Resource (CDR) published in:

> Liu, J. et al. (2018). An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. *Cell*, 173(2), 400–416.e11. https://doi.org/10.1016/j.cell.2018.02.052

The CDR provides four curated survival endpoints that are more reliable than raw GDC fields for survival analysis:

| Column | Description |
|--------|-------------|
| `cdr_OS` | Overall survival event (0/1) |
| `cdr_OS.time` | Overall survival time (days) |
| `cdr_PFI` | Progression-free interval event |
| `cdr_PFI.time` | Progression-free interval time (days) |
| `cdr_DSS` | Disease-specific survival event |
| `cdr_DSS.time` | Disease-specific survival time (days) |
| `cdr_DFI` | Disease-free interval event |
| `cdr_DFI.time` | Disease-free interval time (days) |

Three audit flag columns are added to every sample:

| Column | Description |
|--------|-------------|
| `cdr_matched` | True if case was found in the CDR |
| `cdr_subtype_available` | True if molecular subtype is populated |
| `cdr_survival_complete` | True if OS, OS.time, PFI, PFI.time are all present |

Samples not found in the CDR (e.g. cases added to GDC after the 2018 data freeze) are kept in all output files with `cdr_matched = False` and empty CDR columns — they are never silently dropped.

---

## Requirements

- Python >= 3.10
- pandas >= 2.0
- requests >= 2.31
- openpyxl >= 3.1
- tqdm >= 4.65
- streamlit >= 1.35 *(GUI only)*

---

## Limitations

- Downloads one project at a time (no multi-project batch mode)
- RNA-seq STAR-Counts only (mutations, CNV, methylation not included)
- TCGA projects only (other GDC programs such as CPTAC and TARGET are not supported)
- Open-access files only
- CDR covers the TCGA 2018 data freeze; samples added after 2018 will not have CDR annotations
- Tested on macOS, Linux, and Windows with Python 3.10, 3.11, 3.12, and 3.13

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Author

**Orhan Nedim Kurt**
Independent Researcher

---

## Citation

If you use this tool in your research, please cite:

```
Kurt, O.N. (2026). tcga-gdc-downloader: A tool for downloading and annotating 
TCGA RNA-seq data from the GDC portal. GitHub. 
https://github.com/onedimkurt/tcga-gdc-downloader
```

Please also cite the PanCanAtlas CDR if you use the survival or subtype annotations:

```
Liu, J. et al. (2018). An Integrated TCGA Pan-Cancer Clinical Data Resource 
to Drive High-Quality Survival Outcome Analytics. Cell, 173(2), 400-416.e11.
https://doi.org/10.1016/j.cell.2018.02.052
```

---

## Acknowledgements

Data downloaded via the [NCI Genomic Data Commons (GDC)](https://gdc.cancer.gov) API.
Survival and subtype annotations from the [TCGA PanCanAtlas](https://www.cell.com/pb-assets/consortium/pancanceratlas/pancani3/index.html).
