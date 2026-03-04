"""
constants.py — GDC ID types, endpoints, and download method documentation.

GDC ID TYPE REFERENCE
---------------------
┌────────────────────┬──────────────────────────────────────────────────────────┐
│ ID type            │ Description and example                                  │
├────────────────────┼──────────────────────────────────────────────────────────┤
│ file_id (UUID)     │ GDC UUID for one downloadable file.                      │
│                    │ Format: 8-4-4-4-12 hex chars                             │
│                    │ Example: a3e5f7c2-1234-5678-abcd-ef0123456789           │
│                    │ USE: POST /data  body={"ids":[file_id,...]}              │
├────────────────────┼──────────────────────────────────────────────────────────┤
│ case_id (UUID)     │ GDC UUID for a patient/case. NOT the TCGA barcode.       │
│                    │ Example: b1c2d3e4-0000-1111-2222-333344445555           │
│                    │ USE: clinical data joins                                  │
├────────────────────┼──────────────────────────────────────────────────────────┤
│ case_submitter_id  │ Human TCGA case barcode. 3-field format.                 │
│ (TCGA barcode)     │ Example: TCGA-BH-A0B3                                   │
│                    │ USE: clinical <-> metadata joins                         │
├────────────────────┼──────────────────────────────────────────────────────────┤
│ sample_id (UUID)   │ GDC UUID for a tissue sample.                            │
│                    │ USE: metadata table (informational)                       │
├────────────────────┼──────────────────────────────────────────────────────────┤
│ sample_submitter_id│ Human TCGA sample barcode. 4-field format.               │
│ (TCGA barcode)     │ 4th field encodes tumor/normal: 01=tumor, 11=normal etc. │
│                    │ Example: TCGA-BH-A0B3-01A                                │
│                    │ USE: count matrix column labels; tumor/normal class.      │
└────────────────────┴──────────────────────────────────────────────────────────┘

DOWNLOAD METHOD: REST API /data endpoint (NOT gdc-client manifest)
------------------------------------------------------------------
We POST file_id UUIDs directly to the REST API:
    POST https://api.gdc.cancer.gov/data
    Body: {"ids": ["uuid1", "uuid2", ...]}
    Response: ZIP archive containing one subdirectory per file_id UUID

This is completely different from the gdc-client manifest approach:
    gdc-client download -m manifest.txt
    (requires multi-column TSV: UUID + MD5 + size + filename)

We never create a manifest file and never use gdc-client.
See: https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/

EXTRACTED DIRECTORY STRUCTURE
------------------------------
raw_counts/
    <file_id_uuid>/           <- directory named with file_id UUID
        <filename>.tsv        <- STAR counts file

Primary key for matching files back to metadata: the UUID directory name.
Fallback: the filename.
"""

import re

# ── GDC API endpoints ─────────────────────────────────────────────────────────
GDC_BASE              = "https://api.gdc.cancer.gov"
GDC_FILES_ENDPOINT    = f"{GDC_BASE}/files"
GDC_CASES_ENDPOINT    = f"{GDC_BASE}/cases"
GDC_DATA_ENDPOINT     = f"{GDC_BASE}/data"
GDC_PROJECTS_ENDPOINT = f"{GDC_BASE}/projects"
GDC_STATUS_ENDPOINT   = f"{GDC_BASE}/status"

# ── UUID validation ───────────────────────────────────────────────────────────
_UUID_RE = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$",
    re.IGNORECASE,
)

def is_valid_gdc_uuid(value: str) -> bool:
    """Return True if value matches the GDC UUID format."""
    return bool(_UUID_RE.match(str(value).strip()))


# ── TCGA sample type codes ────────────────────────────────────────────────────
SAMPLE_TYPE_CODES: dict = {
    "01": "Primary Tumor",
    "02": "Recurrent Tumor",
    "03": "Primary Blood Derived Cancer - Peripheral Blood",
    "04": "Recurrent Blood Derived Cancer - Bone Marrow",
    "05": "Additional - New Primary",
    "06": "Metastatic",
    "07": "Additional Metastatic",
    "08": "Human Tumor Original Cells",
    "09": "Primary Blood Derived Cancer - Bone Marrow",
    "10": "Blood Derived Normal",
    "11": "Solid Tissue Normal",
    "12": "Buccal Cell Normal",
    "13": "EBV Immortalized Normal",
    "14": "Bone Marrow Normal",
    "20": "Control Analyte",
    "40": "Recurrent Blood Derived Cancer - Peripheral Blood",
    "50": "Cell Lines",
    "60": "Primary Xenograft Tissue",
    "61": "Cell Line Derived Xenograft Tissue",
}

TUMOR_CODES:  frozenset = frozenset(
    {"01","02","03","04","05","06","07","08","09","40","50","60","61"}
)
NORMAL_CODES: frozenset = frozenset(
    {"10","11","12","13","14","20"}
)

# ── STAR-Counts filter specification ─────────────────────────────────────────
# All four constraints are required together to specifically identify
# STAR unstranded counts TSV files (not FPKM, not miRNA, not WGS, etc.)
STAR_COUNTS_FILTERS = [
    {"op": "=", "content": {"field": "data_type",
                             "value": "Gene Expression Quantification"}},
    {"op": "=", "content": {"field": "experimental_strategy",
                             "value": "RNA-Seq"}},
    {"op": "=", "content": {"field": "analysis.workflow_type",
                             "value": "STAR - Counts"}},
    {"op": "=", "content": {"field": "data_format",
                             "value": "TSV"}},
    {"op": "=", "content": {"field": "access",
                             "value": "open"}},
]

# ── Fields requested from GDC /files ─────────────────────────────────────────
# file_id   = UUID  → used in /data download endpoint
# file_name = str   → used to match extracted TSV files
# file_size = int   → used to estimate total download size
FILE_FIELDS = (
    "file_id,file_name,file_size,"
    "cases.case_id,cases.submitter_id,"
    "cases.samples.sample_id,cases.samples.submitter_id,"
    "cases.samples.sample_type,cases.samples.tissue_type"
)

# ── Fields requested from GDC /cases ─────────────────────────────────────────
CLINICAL_FIELDS = ",".join([
    "case_id","submitter_id",
    "demographic.gender","demographic.age_at_index","demographic.race",
    "demographic.ethnicity","demographic.vital_status",
    "demographic.days_to_death","demographic.year_of_birth",
    "diagnoses.age_at_diagnosis","diagnoses.days_to_last_follow_up",
    "diagnoses.days_to_death","diagnoses.primary_diagnosis",
    "diagnoses.tumor_stage","diagnoses.tumor_grade",
    "diagnoses.morphology","diagnoses.tissue_or_organ_of_origin",
    "diagnoses.site_of_resection_or_biopsy",
    "diagnoses.prior_malignancy","diagnoses.synchronous_malignancy",
    "diagnoses.ajcc_pathologic_stage","diagnoses.ajcc_pathologic_t",
    "diagnoses.ajcc_pathologic_n","diagnoses.ajcc_pathologic_m",
    "exposures.cigarettes_per_day","exposures.years_smoked",
    "exposures.pack_years_smoked","exposures.alcohol_history",
    "exposures.bmi","exposures.height","exposures.weight",
])

# ── STAR column name synonyms (ordered by preference) ────────────────────────
# Current GDC format:  gene_id  gene_name  gene_type  unstranded  stranded_first  stranded_second
# Older GDC format:    gene_id  unstranded  stranded_first  stranded_second
UNSTRANDED_SYNONYMS = ("unstranded","stranded_unspecific","count","expected_count","raw_count")
STRANDED1_SYNONYMS  = ("stranded_first","sense","stranded")
STRANDED2_SYNONYMS  = ("stranded_second","antisense")

# ── GDC pagination limits ─────────────────────────────────────────────────────
GDC_PAGE_SIZE   = 10_000   # GDC hard cap per page
MAX_TOTAL_FILES = 100_000  # safety ceiling

# ── Request settings ──────────────────────────────────────────────────────────
REQUEST_TIMEOUT_SHORT = 60
REQUEST_TIMEOUT_LONG  = 900
MAX_RETRIES           = 3
RETRY_WAIT_SECONDS    = 10
CHUNK_SIZE_BYTES      = 1024 * 1024

# GDC bulk download truncates archives above ~100 MB.
# Batching at 50 files keeps each request well under this limit.
GDC_DOWNLOAD_BATCH_SIZE = 50

# ── Checkpoint ────────────────────────────────────────────────────────────────
CHECKPOINT_FILE = "pipeline_checkpoint.json"

# ── Program scope ─────────────────────────────────────────────────────────────
# This tool supports TCGA projects only.
# CDR integration is only available for TCGA.
SUPPORTED_PROGRAMS:     frozenset = frozenset({"TCGA"})
CDR_SUPPORTED_PROGRAMS: frozenset = frozenset({"TCGA"})

# ── PanCanAtlas CDR ───────────────────────────────────────────────────────────
# The TCGA Pan-Cancer Clinical Data Resource (Liu et al., Cell 2018).
# Stable GDC file UUID — this file does not change.
# Download via: POST /data  body={"ids": [CDR_FILE_UUID]}
CDR_FILE_UUID     = "1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
CDR_CACHE_FILENAME = "TCGA_CDR_SupplementalTableS1.xlsx"

# Join key: CDR uses 'bcr_patient_barcode' = our 'case_submitter_id'
CDR_JOIN_KEY = "bcr_patient_barcode"

# Survival endpoint columns in CDR (used to compute cdr_survival_complete flag)
CDR_SURVIVAL_COLS = ("OS", "OS.time", "DSS", "DSS.time",
                     "DFI", "DFI.time", "PFI", "PFI.time")

# Subtype column in CDR (used to compute cdr_subtype_available flag).
# The GDC-hosted CDR uses "Subtype_Selected"; the original Cell 2018 paper
# supplement used "Subtype_mRNA". We check both, preferring Subtype_Selected.
CDR_SUBTYPE_COL         = "Subtype_Selected"    # primary (GDC-hosted version)
CDR_SUBTYPE_COL_LEGACY  = "Subtype_mRNA"        # fallback (original paper supplement)

# All CDR columns that will be prefixed with 'cdr_' in output
# Covers: survival, subtype, harmonised clinical, treatment
CDR_COLUMNS_TO_KEEP = [
    "bcr_patient_barcode",   # join key — kept unprefixed internally, dropped in output
    "type",                  # TCGA cancer type abbreviation (BRCA, GBM, etc.)
    # Demographics
    "age_at_initial_pathologic_diagnosis",
    "gender",
    "race",
    "birth_days_to",
    "menopause_status",
    # Tumour characteristics
    "ajcc_pathologic_tumor_stage",
    "clinical_stage",
    "histological_type",
    "histological_grade",
    "initial_pathologic_dx_year",
    # Outcome fields
    "vital_status",
    "tumor_status",
    "last_contact_days_to",
    "death_days_to",
    "cause_of_death",
    "treatment_outcome_first_course",
    "margin_status",
    "residual_tumor",
    # New tumour events
    "new_tumor_event_type",
    "new_tumor_event_site",
    "new_tumor_event_site_other",
    "new_tumor_event_dx_days_to",
    # Survival endpoints (four curated endpoints from Liu et al. Cell 2018)
    "OS", "OS.time",
    "DSS", "DSS.time",
    "DFI", "DFI.time",
    "PFI", "PFI.time",
    # Molecular subtypes
    # Note: the GDC-hosted CDR uses Subtype_Selected and Subtype_Integrative.
    # The per-platform columns (Subtype_mRNA, Subtype_miRNA, Subtype_protein,
    # Subtype_CNV_low, Subtype_CNV_high) are present only in the original
    # Cell 2018 paper supplement, not in the GDC-hosted version.
    "Subtype_Selected",      # Harmonised subtype chosen by CDR authors
    "Subtype_Integrative",   # Pan-cancer iCluster subtype
    # Keep legacy names as fallbacks — silently skipped if not present
    "Subtype_mRNA",
    "Subtype_miRNA",
    "Subtype_protein",
    "Subtype_CNV_low",
    "Subtype_CNV_high",
]

# GDC /files fields addition — add program.name so we can gate non-TCGA projects
# (Added here as a reminder — client.py validate_project uses this)
PROJECTS_FIELDS = "project_id,name,primary_site,disease_type,summary,program.name"
