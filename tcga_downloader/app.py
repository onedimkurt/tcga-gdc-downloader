"""
app.py  —  Streamlit graphical interface for tcga_downloader
============================================================
Launch with:
    streamlit run app.py
Or via the CLI:
    tcga-download --gui

Step order
----------
1  Token             — upload GDC token (or skip for open-access)
2  Project           — select TCGA project
3  Discover          — find STAR-Counts files and show tumor/normal breakdown
4  Output Dir        — choose where to save files
5  Download          — bulk download via GDC /data endpoint
6  Clinical          — fetch GDC clinical data
7  CDR               — download and merge PanCanAtlas CDR annotations
8  Build Matrix      — parse TSV files, assemble gene × sample matrix
9  Merge & Save      — join all data, write output files
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import streamlit as st

st.set_page_config(
    page_title="TCGA GDC Downloader",
    page_icon="🧬",
    layout="centered",
    initial_sidebar_state="collapsed",
)

from tcga_downloader.client     import GDCClient
from tcga_downloader.merge      import (flatten_file_metadata, flatten_clinical,
                                        merge_outputs, save_full_merged_with_cdr)
from tcga_downloader.matrix     import build_count_matrix
from tcga_downloader.sample     import summarise_sample_types
from tcga_downloader.checkpoint import Checkpoint
from tcga_downloader.exceptions import GDCError
from tcga_downloader.cdr        import run_cdr_pipeline


# =============================================================================
#   Styling
# =============================================================================

st.markdown("""
<style>
    .main { max-width: 800px; margin: auto; }
    .step-header { font-size: 1.1rem; font-weight: 600; margin-top: 1.2rem; }
    .info-box    { background:#f0f7ff; border-left:4px solid #1f77b4;
                   padding:.8rem 1rem; border-radius:4px; margin:.5rem 0; }
    .warn-box    { background:#fff8e1; border-left:4px solid #f9a825;
                   padding:.8rem 1rem; border-radius:4px; margin:.5rem 0; }
    .success-box { background:#e8f5e9; border-left:4px solid #388e3c;
                   padding:.8rem 1rem; border-radius:4px; margin:.5rem 0; }
    .error-box   { background:#ffebee; border-left:4px solid #c62828;
                   padding:.8rem 1rem; border-radius:4px; margin:.5rem 0; }
    .cdr-box     { background:#f3e5f5; border-left:4px solid #7b1fa2;
                   padding:.8rem 1rem; border-radius:4px; margin:.5rem 0; }
</style>
""", unsafe_allow_html=True)


# =============================================================================
#   Session state
# =============================================================================

def _init_state():
    defaults = {
        "step":             1,
        "client":           None,
        "project_id":       None,
        "project_info":     None,
        "program_name":     "TCGA",
        "file_hits":        None,
        "metadata_df":      None,
        "metadata_with_cdr": None,   # metadata enriched with CDR columns
        "sample_summary":   None,
        "output_dir":       None,
        "checkpoint":       None,
        "clinical_df":      None,
        "cdr_result":       None,    # full dict from run_cdr_pipeline()
        "matrix":           None,
        "output_paths":     None,
        "error":            None,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

_init_state()


def _show_error(err: GDCError | Exception):
    if isinstance(err, GDCError):
        st.markdown(f"""
<div class="error-box">
<b>❌ Error:</b> {err.message}<br><br>
<b>💡 How to fix it:</b><br>
<pre style="margin:0;white-space:pre-wrap">{err.fix}</pre>
</div>""", unsafe_allow_html=True)
    else:
        st.error(f"Unexpected error: {type(err).__name__}: {err}")


# =============================================================================
#   Header and step progress indicator
# =============================================================================

st.title("🧬 TCGA GDC RNA-seq Downloader")
st.markdown(
    "Download STAR unstranded counts from the **NCI GDC portal**, "
    "merged with clinical data and PanCanAtlas CDR annotations — no coding required."
)
st.divider()

STEP_NAMES = [
    "Token", "Project", "Discover", "Output Dir",
    "Download", "Clinical", "CDR", "Build Matrix", "Merge & Save"
]
TOTAL_STEPS = len(STEP_NAMES)
current_step = st.session_state["step"]

cols = st.columns(TOTAL_STEPS)
for i, (col, name) in enumerate(zip(cols, STEP_NAMES)):
    n = i + 1
    if n < current_step:
        col.markdown(
            f"<div style='text-align:center;color:#388e3c;font-size:.8rem'>"
            f"✅<br>{name}</div>", unsafe_allow_html=True)
    elif n == current_step:
        col.markdown(
            f"<div style='text-align:center;color:#1f77b4;font-size:.8rem'>"
            f"<b>▶<br>{name}</b></div>", unsafe_allow_html=True)
    else:
        col.markdown(
            f"<div style='text-align:center;color:#aaa;font-size:.8rem'>"
            f"○<br>{name}</div>", unsafe_allow_html=True)

st.divider()


# =============================================================================
#   STEP 1 — Authentication
# =============================================================================

if current_step == 1:
    st.markdown('<p class="step-header">Step 1 — GDC Authentication Token</p>',
                unsafe_allow_html=True)
    st.markdown("""
<div class="info-box">
<b>How to get your GDC token:</b><br>
1. Go to <a href="https://portal.gdc.cancer.gov" target="_blank">portal.gdc.cancer.gov</a><br>
2. Click <b>Login</b> (top right) and sign in (free account)<br>
3. Click your <b>username → Download Token</b><br>
4. Upload that file below
</div>""", unsafe_allow_html=True)

    token_file = st.file_uploader("Upload your GDC token file (.txt)", type=["txt"])

    col1, col2 = st.columns(2)
    with col1:
        if token_file is not None:
            if st.button("✅ Load Token", use_container_width=True, type="primary"):
                try:
                    token = token_file.read().decode("utf-8").strip()
                    if len(token) < 20:
                        st.error("Token file appears empty or invalid.")
                    else:
                        st.session_state["client"] = GDCClient(token=token)
                        st.session_state["step"]   = 2
                        st.rerun()
                except Exception as e:
                    st.error(f"Could not read token file: {e}")
    with col2:
        if st.button("⏭ Skip (open-access only)", use_container_width=True):
            st.session_state["client"] = GDCClient()
            st.session_state["step"]   = 2
            st.rerun()


# =============================================================================
#   STEP 2 — Project Selection
# =============================================================================

elif current_step == 2:
    st.markdown('<p class="step-header">Step 2 — Select a TCGA Project</p>',
                unsafe_allow_html=True)

    common_projects = [
        "TCGA-BRCA", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-COAD",
        "TCGA-KIRC", "TCGA-LIHC", "TCGA-STAD", "TCGA-BLCA", "TCGA-HNSC",
        "TCGA-GBM",  "TCGA-OV",   "TCGA-UCEC", "TCGA-THCA", "TCGA-SKCM",
        "TCGA-KIRP", "TCGA-SARC", "TCGA-LAML", "TCGA-ESCA", "TCGA-PAAD",
        "TCGA-CHOL", "TCGA-MESO", "TCGA-UVM",  "TCGA-ACC",  "TCGA-PCPG",
        "TCGA-TGCT", "TCGA-THYM", "TCGA-DLBC", "TCGA-UCS",  "TCGA-READ",
        "TCGA-LGG",  "TCGA-CESC", "TCGA-KICH",
        "Other (type below)",
    ]

    selected = st.selectbox("Choose a project:", common_projects)
    custom   = ""
    if selected == "Other (type below)":
        custom = st.text_input("Enter project ID:", placeholder="TCGA-XXXX").strip().upper()

    pid = custom if selected == "Other (type below)" else selected

    if st.button("▶ Validate Project", type="primary", disabled=not pid):
        with st.spinner(f"Checking {pid}..."):
            try:
                info = st.session_state["client"].validate_project(pid)
                program = info.get("program", {}).get("name", "TCGA")
                st.session_state["project_id"]   = pid
                st.session_state["project_info"] = info
                st.session_state["program_name"] = program
                st.session_state["step"]         = 3
                st.rerun()
            except GDCError as e:
                _show_error(e)


# =============================================================================
#   STEP 3 — Discover Files
# =============================================================================

elif current_step == 3:
    pid = st.session_state["project_id"]
    st.markdown(f'<p class="step-header">Step 3 — Discovering STAR-Counts Files for {pid}</p>',
                unsafe_allow_html=True)

    with st.spinner("Querying GDC for STAR-Counts RNA-Seq files..."):
        try:
            hits        = st.session_state["client"].discover_star_files(pid)
            metadata_df = flatten_file_metadata(hits)
            summary     = summarise_sample_types(metadata_df)

            st.session_state["file_hits"]   = hits
            st.session_state["metadata_df"] = metadata_df
            st.session_state["sample_summary"] = summary

            n_tumor  = summary.get("Tumor",  0)
            n_normal = summary.get("Normal", 0)
            n_total  = len(hits)
            total_gb = sum(f.get("file_size", 0) for f in hits) / (1024**3)

            st.markdown(f"""
<div class="success-box">
<b>✅  Found {n_total} files  (~{total_gb:.1f} GB)</b><br>
Tumor samples: <b>{n_tumor}</b> &nbsp;|&nbsp; Normal samples: <b>{n_normal}</b>
</div>""", unsafe_allow_html=True)

            if n_normal > 0:
                st.markdown("""
<div class="warn-box">
⚠️ <b>Normal samples included.</b> Output will contain a
<code>sample_category</code> column (Tumor / Normal) and separate
files for each group. Always verify this column before analysis.
</div>""", unsafe_allow_html=True)

            st.session_state["step"] = 4
            st.rerun()
        except GDCError as e:
            _show_error(e)


# =============================================================================
#   STEP 4 — Output Directory
# =============================================================================

elif current_step == 4:
    pid  = st.session_state["project_id"]
    hits = st.session_state["file_hits"]
    st.markdown('<p class="step-header">Step 4 — Choose Output Directory</p>',
                unsafe_allow_html=True)

    total_gb = sum(f.get("file_size", 0) for f in hits) / (1024**3)
    st.caption(f"Estimated download size: ~{total_gb:.1f} GB")

    default_dir = str(Path.home() / "TCGA_downloads" / pid)
    out_str = st.text_input("Output directory:", value=default_dir)

    if st.button("▶ Start Download", type="primary"):
        out_path = Path(out_str.strip())
        try:
            out_path.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            st.error(f"Cannot create directory: {e}")
            st.stop()
        st.session_state["output_dir"]  = out_path
        st.session_state["checkpoint"]  = Checkpoint(out_path)
        st.session_state["step"]        = 5
        st.rerun()


# =============================================================================
#   STEP 5 — Download
# =============================================================================

elif current_step == 5:
    hits       = st.session_state["file_hits"]
    output_dir = st.session_state["output_dir"]
    cp         = st.session_state["checkpoint"]

    st.markdown('<p class="step-header">Step 5 — Downloading Files</p>',
                unsafe_allow_html=True)

    if cp.is_done("downloaded"):
        st.success("✅ Download already completed (checkpoint found).")
        st.session_state["step"] = 6
        st.rerun()

    import tarfile
    from tcga_downloader.constants import GDC_DOWNLOAD_BATCH_SIZE
    counts_dir = output_dir / "raw_counts"
    counts_dir.mkdir(exist_ok=True)
    file_ids  = [f["file_id"] for f in hits]
    n_batches = (len(file_ids) + GDC_DOWNLOAD_BATCH_SIZE - 1) // GDC_DOWNLOAD_BATCH_SIZE

    try:
        progress_bar = st.progress(0, text=f"Downloading batch 1 of {n_batches}...")
        for batch_num, batch_start in enumerate(
                range(0, len(file_ids), GDC_DOWNLOAD_BATCH_SIZE), start=1):
            batch_ids  = file_ids[batch_start:batch_start + GDC_DOWNLOAD_BATCH_SIZE]
            batch_path = output_dir / f"gdc_batch_{batch_num:03d}.zip"

            progress_bar.progress(
                int((batch_num - 1) / n_batches * 100),
                text=f"Downloading batch {batch_num}/{n_batches} ({len(batch_ids)} files)..."
            )

            if not batch_path.exists():
                st.session_state["client"].stream_download(batch_ids, str(batch_path))

            # Extract batch
            with open(batch_path, "rb") as fh:
                magic = fh.read(4)
            try:
                if magic[:2] == b'\x1f\x8b':
                    archive_path = batch_path.with_suffix(".tar.gz")
                    batch_path.rename(archive_path)
                    with tarfile.open(archive_path, "r:gz") as t:
                        t.extractall(counts_dir)
                    archive_path.unlink(missing_ok=True)
                else:
                    with zipfile.ZipFile(batch_path, "r") as z:
                        z.extractall(counts_dir)
                    batch_path.unlink(missing_ok=True)
            except (EOFError, tarfile.ReadError, zipfile.BadZipFile) as extract_err:
                for p in [batch_path, batch_path.with_suffix(".tar.gz")]:
                    if p.exists():
                        p.unlink()
                st.error(
                    f"⚠️ Batch {batch_num}/{n_batches} was incomplete or corrupted. "
                    f"The partial file has been deleted. "
                    f"Click **Re-run step** to resume from this batch."
                )
                st.stop()

        progress_bar.progress(100, text=f"✅ All {n_batches} batches downloaded and extracted.")
        cp.save("downloaded", {"counts_dir": str(counts_dir)})
        st.session_state["step"] = 6
        st.rerun()
    except GDCError as e:
        _show_error(e)


# =============================================================================
#   STEP 6 — Clinical Data
# =============================================================================

elif current_step == 6:
    pid        = st.session_state["project_id"]
    output_dir = st.session_state["output_dir"]
    cp         = st.session_state["checkpoint"]

    st.markdown('<p class="step-header">Step 6 — Fetching Clinical Data</p>',
                unsafe_allow_html=True)

    clinical_cache = output_dir / "clinical_data.tsv"

    if cp.is_done("clinical") and clinical_cache.exists():
        import pandas as pd
        st.session_state["clinical_df"] = pd.read_csv(clinical_cache, sep="\t", dtype=str)
        st.success(f"✅ Clinical data loaded ({len(st.session_state['clinical_df'])} cases).")
        st.session_state["step"] = 7
        st.rerun()

    with st.spinner("Retrieving demographics, diagnoses, and survival data from GDC..."):
        try:
            case_hits   = st.session_state["client"].fetch_clinical_data(pid)
            clin_df     = flatten_clinical(case_hits)
            clin_df.to_csv(clinical_cache, sep="\t", index=False)
            cp.save("clinical", {"n_cases": len(clin_df)})
            st.session_state["clinical_df"] = clin_df
            st.success(f"✅  {len(clin_df)} clinical records fetched.")
            st.session_state["step"] = 7
            st.rerun()
        except GDCError as e:
            _show_error(e)


# =============================================================================
#   STEP 7 — PanCanAtlas CDR Annotations
# =============================================================================

elif current_step == 7:
    pid          = st.session_state["project_id"]
    program_name = st.session_state["program_name"]
    metadata_df  = st.session_state["metadata_df"]
    output_dir   = st.session_state["output_dir"]
    cp           = st.session_state["checkpoint"]

    st.markdown('<p class="step-header">Step 7 — PanCanAtlas CDR Annotations</p>',
                unsafe_allow_html=True)

    st.markdown("""
<div class="cdr-box">
<b>📚 What is the CDR?</b><br>
The <b>TCGA Pan-Cancer Clinical Data Resource</b> (Liu et al., <i>Cell</i> 2018)
provides harmonised clinical annotations for all 33 TCGA cancer types, including:<br>
• Four curated survival endpoints (OS, DSS, DFI, PFI)<br>
• Molecular subtypes (PAM50 for BRCA, IDH status for GBM, etc.)<br>
• Harmonised tumour stage, histological type, and demographics<br><br>
CDR columns are prefixed <code>cdr_</code> to distinguish them from GDC API fields.
Cases not in the CDR (post-2018 additions) are kept with <code>cdr_matched = False</code>.
</div>""", unsafe_allow_html=True)

    cdr_cache = output_dir / "metadata_with_cdr.tsv"

    # ── Already done ──────────────────────────────────────────────────────────
    if cp.is_done("cdr"):
        import pandas as pd
        saved = cp.get("cdr")
        if saved.get("skipped"):
            st.warning(f"⚠️  CDR step was previously skipped: {saved.get('reason', '')}")
            st.session_state["metadata_with_cdr"] = metadata_df.copy()
        elif cdr_cache.exists():
            n_matched = saved.get("n_matched", "?")
            n_total   = saved.get("n_total",   "?")
            st.success(f"✅ CDR already merged ({n_matched}/{n_total} cases matched).")
            st.session_state["metadata_with_cdr"] = pd.read_csv(cdr_cache, sep="\t", dtype=str)
        else:
            st.success("✅ CDR: already done.")
            st.session_state["metadata_with_cdr"] = metadata_df.copy()
        st.session_state["step"] = 8
        st.rerun()

    # ── Run CDR pipeline ──────────────────────────────────────────────────────
    col1, col2 = st.columns(2)
    run_cdr   = col1.button("▶ Fetch CDR Annotations", type="primary",
                             use_container_width=True)
    skip_cdr  = col2.button("⏭ Skip CDR step", use_container_width=True)

    if skip_cdr:
        cp.save("cdr", {"skipped": True, "reason": "User chose to skip"})
        st.session_state["metadata_with_cdr"] = metadata_df.copy()
        st.session_state["step"] = 8
        st.rerun()

    if run_cdr:
        with st.spinner("Downloading PanCanAtlas CDR file (~2 MB, cached after first use)..."):
            try:
                cdr_result = run_cdr_pipeline(
                    client       = st.session_state["client"],
                    metadata_df  = metadata_df,
                    project_id   = pid,
                    program_name = program_name,
                    output_dir   = output_dir,
                    verbose      = False,
                )
                merged_df = cdr_result["merged_df"]
                merged_df.to_csv(cdr_cache, sep="\t", index=False)

                n_matched  = int(merged_df["cdr_matched"].sum()) if "cdr_matched" in merged_df.columns else 0
                n_total    = len(merged_df)
                n_subtype  = int(merged_df["cdr_subtype_available"].sum()) if "cdr_subtype_available" in merged_df.columns else 0
                n_survival = int(merged_df["cdr_survival_complete"].sum()) if "cdr_survival_complete" in merged_df.columns else 0

                cp.save("cdr", {
                    "n_matched":  n_matched,
                    "n_total":    n_total,
                    "cache":      str(cdr_cache),
                    "output_paths": {k: str(v) for k, v in cdr_result["output_paths"].items()},
                })

                st.session_state["cdr_result"]       = cdr_result
                st.session_state["metadata_with_cdr"] = merged_df

                # ── CDR summary ───────────────────────────────────────────────
                st.markdown(f"""
<div class="success-box">
<b>✅  CDR annotations merged</b><br>
Cases matched to CDR: <b>{n_matched} / {n_total}</b>
({100*n_matched/n_total:.1f}%)<br>
Subtype available: <b>{n_subtype}</b> &nbsp;|&nbsp;
Survival complete (OS + PFI): <b>{n_survival}</b>
</div>""", unsafe_allow_html=True)

                if n_matched < n_total:
                    st.markdown(f"""
<div class="warn-box">
⚠️ <b>{n_total - n_matched} case(s) not found in CDR.</b><br>
These are likely samples added to the GDC after the 2018 TCGA data freeze.
They are kept in all output files with <code>cdr_matched = False</code> and
empty CDR columns. Their count data and GDC clinical data are unaffected.
</div>""", unsafe_allow_html=True)

                # Show coverage report as table
                cov_path = output_dir / f"{pid}_CDR_coverage_report.tsv"
                if cov_path.exists():
                    import pandas as pd
                    cov_df = pd.read_csv(cov_path, sep="\t")
                    with st.expander("📊 View CDR field coverage"):
                        st.dataframe(
                            cov_df.style.format({"pct_present": "{:.1f}%"}),
                            use_container_width=True,
                            height=300,
                        )

                st.session_state["step"] = 8
                st.rerun()

            except GDCError as e:
                _show_error(e)
                st.markdown("""
<div class="warn-box">
You can skip the CDR step and continue with GDC clinical data only.
</div>""", unsafe_allow_html=True)


# =============================================================================
#   STEP 8 — Build Count Matrix
# =============================================================================

elif current_step == 8:
    output_dir   = st.session_state["output_dir"]
    metadata_df  = st.session_state["metadata_df"]
    cp           = st.session_state["checkpoint"]
    counts_dir   = output_dir / "raw_counts"
    matrix_cache = output_dir / "count_matrix_cache.tsv"

    st.markdown('<p class="step-header">Step 8 — Building Count Matrix</p>',
                unsafe_allow_html=True)

    if cp.is_done("matrix") and matrix_cache.exists():
        import pandas as pd
        with st.spinner("Loading cached count matrix..."):
            matrix = pd.read_csv(matrix_cache, sep="\t", index_col=0)
        st.success(f"✅ Matrix loaded: {matrix.shape[0]:,} genes × {matrix.shape[1]:,} samples")
        st.session_state["matrix"] = matrix
        st.session_state["step"]   = 9
        st.rerun()

    with st.spinner("Parsing count files and assembling matrix "
                    "(may take several minutes for large projects)..."):
        try:
            matrix = build_count_matrix(counts_dir, metadata_df, verbose=False)
            matrix.to_csv(matrix_cache, sep="\t")
            cp.save("matrix", {"cache": str(matrix_cache),
                                "shape": list(matrix.shape)})
            st.success(f"✅ Matrix built: {matrix.shape[0]:,} genes × {matrix.shape[1]:,} samples")
            st.session_state["matrix"] = matrix
            st.session_state["step"]   = 9
            st.rerun()
        except GDCError as e:
            _show_error(e)


# =============================================================================
#   STEP 9 — Merge and Save
# =============================================================================

elif current_step == 9:
    pid               = st.session_state["project_id"]
    matrix            = st.session_state["matrix"]
    metadata_df       = st.session_state["metadata_df"]
    _cdr = st.session_state.get("metadata_with_cdr")
    metadata_with_cdr = _cdr if _cdr is not None else metadata_df
    clinical_df       = st.session_state["clinical_df"]
    output_dir        = st.session_state["output_dir"]
    cp                = st.session_state["checkpoint"]

    st.markdown('<p class="step-header">Step 9 — Merging and Saving Output</p>',
                unsafe_allow_html=True)

    if cp.is_done("merged"):
        paths = {k: Path(v) for k, v in cp.get("merged").items()}
        st.markdown('<div class="success-box">🎉 All done! Your files are ready.</div>',
                    unsafe_allow_html=True)
        st.session_state["output_paths"] = paths
    else:
        with st.spinner("Joining all data and writing output files..."):
            try:
                # Standard output: counts + GDC clinical
                paths = merge_outputs(matrix, metadata_df, clinical_df,
                                      pid, output_dir, verbose=False)

                # CDR-enriched output (if CDR was run)
                has_cdr = any(c.startswith("cdr_")
                              for c in metadata_with_cdr.columns)
                if has_cdr:
                    cdr_paths = save_full_merged_with_cdr(
                        count_matrix      = matrix,
                        metadata_with_cdr = metadata_with_cdr,
                        project_id        = pid,
                        output_dir        = output_dir,
                        verbose           = False,
                    )
                    paths.update(cdr_paths)

                cp.save("merged", {k: str(v) for k, v in paths.items()})
                st.session_state["output_paths"] = paths
                st.balloons()
            except GDCError as e:
                _show_error(e)
                st.stop()

    # ── Output file cards ─────────────────────────────────────────────────────
    output_paths = st.session_state.get("output_paths", {})
    if output_paths:
        st.divider()
        st.subheader("📂 Output Files")

        file_labels = {
            "all":          ("All samples — counts + GDC clinical",
                             "⚠️ Contains both Tumor and Normal rows"),
            "tumor":        ("Tumor samples only — counts + GDC clinical", ""),
            "normal":       ("Normal samples only — counts + GDC clinical", ""),
            "metadata":     ("Metadata + GDC clinical only (no counts)", ""),
            "full_with_cdr":("All samples — counts + GDC clinical + CDR",
                             "🏆 Recommended for survival and subtype analyses"),
        }

        for key, (description, note) in file_labels.items():
            path = output_paths.get(key)
            if path and Path(path).exists():
                size_mb = Path(path).stat().st_size / (1024**2)
                col1, col2 = st.columns([3, 1])
                note_html = f"<br><small style='color:#888'>{note}</small>" if note else ""
                col1.markdown(
                    f"**{description}**{note_html}  \n"
                    f"`{Path(path).name}` ({size_mb:.1f} MB)",
                    unsafe_allow_html=True,
                )
                with open(path, "rb") as f:
                    col2.download_button(
                        label     = "⬇ Download",
                        data      = f,
                        file_name = Path(path).name,
                        mime      = "text/tab-separated-values",
                        key       = f"dl_{key}",
                    )

        # Also show CDR-specific files if they exist
        cdr_extra = {
            "cdr_annotations": "CDR annotations only (no counts)",
            "coverage_report": "CDR field coverage report",
            "complete_cases":  "CDR-complete cases only",
            "unmatched_cases": "Unmatched cases (not in CDR)",
        }
        cdr_paths = st.session_state.get("cdr_result", {}).get("output_paths", {})
        if cdr_paths:
            st.markdown("**CDR supporting files:**")
            for key, label in cdr_extra.items():
                path = cdr_paths.get(key)
                if path and Path(path).exists():
                    size_mb = Path(path).stat().st_size / (1024**2)
                    col1, col2 = st.columns([3, 1])
                    col1.markdown(f"{label}  \n`{Path(path).name}` ({size_mb:.1f} MB)")
                    with open(path, "rb") as f:
                        col2.download_button(
                            label     = "⬇ Download",
                            data      = f,
                            file_name = Path(path).name,
                            mime      = "text/tab-separated-values",
                            key       = f"dl_cdr_{key}",
                        )

        st.divider()
        st.markdown("""
**How to open your files:**

| Tool | Command |
|------|---------|
| Python | `df = pd.read_csv("file.tsv", sep="\\t", index_col=0)` |
| R | `df <- read.table("file.tsv", sep="\\t", header=TRUE, row.names=1)` |
| Excel | File → Open → select file → choose Tab delimited |

> **Recommended for survival analysis:** use the `_FULL_merged_with_CDR.tsv` file
> and the `cdr_OS`, `cdr_OS.time`, `cdr_PFI`, `cdr_PFI.time` columns.
> Filter to `cdr_survival_complete == True` for a clean analytical cohort.
>
> **Always** verify the `sample_category` column (Tumor / Normal) before
> running differential expression analyses.
""")

        if st.button("🔁 Download another project"):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.rerun()
