"""
cli.py
======
Command-line interface for tcga_downloader.

Pipeline step order
-------------------
1. auth          — load token (optional)
2. project       — validate project ID, detect program name
3. files_found   — discover STAR-Counts files via GDC /files API
4. downloaded    — bulk download via GDC /data endpoint (UUID list, NOT manifest)
5. clinical      — fetch clinical data via GDC /cases API
6. cdr           — download & merge PanCanAtlas CDR annotations (TCGA only)
7. matrix        — parse TSV files, build gene × sample count matrix
8. merged        — join counts + GDC clinical + CDR, write output files

Usage examples
--------------
tcga-download --gui
tcga-download --project TCGA-BRCA --output ~/my_data
tcga-download --project TCGA-BRCA --token ~/gdc-user-token.txt
tcga-download --project TCGA-BRCA --dry-run
tcga-download --project TCGA-BRCA --fresh
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

from tcga_downloader.exceptions import GDCError


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tcga-download",
        description=(
            "Download TCGA STAR RNA-seq counts from the GDC portal.\n"
            "Run without arguments to use the interactive wizard.\n"
            "Use --gui for the graphical interface."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--project", "-p", metavar="PROJECT_ID",
                        help="TCGA project ID, e.g. TCGA-BRCA")
    parser.add_argument("--output",  "-o", metavar="DIR",
                        help="Output directory (default: ~/TCGA_downloads/PROJECT_ID)")
    parser.add_argument("--token",   "-t", metavar="TOKEN_FILE",
                        help="Path to your GDC authentication token file")
    parser.add_argument("--gui", action="store_true",
                        help="Launch the Streamlit graphical interface")
    parser.add_argument("--dry-run", action="store_true",
                        help="Discover files and report counts without downloading")
    parser.add_argument("--fresh", action="store_true",
                        help="Ignore any existing checkpoint and start from scratch")
    parser.add_argument("--no-cdr", action="store_true",
                        help="Skip the PanCanAtlas CDR annotation step")
    parser.add_argument("--version", "-v", action="version", version="%(prog)s 2.1.2")
    return parser


def run_gui():
    """Launch the Streamlit GUI."""
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print("\n❌  The GUI requires Streamlit.")
        print("\n    Install it with:\n\n        pip install tcga-gdc-downloader[gui]\n")
        sys.exit(1)
    import subprocess
    from tcga_downloader import _gui_path
    subprocess.run([sys.executable, "-m", "streamlit", "run", str(_gui_path())], check=False)


def run_cli(args: argparse.Namespace):
    """Run non-interactively with command-line arguments."""
    import pandas as pd
    import zipfile

    from tcga_downloader.client     import GDCClient
    from tcga_downloader.merge      import (flatten_file_metadata, flatten_clinical,
                                            merge_outputs, save_full_merged_with_cdr)
    from tcga_downloader.matrix     import build_count_matrix
    from tcga_downloader.checkpoint import Checkpoint
    from tcga_downloader.cdr        import run_cdr_pipeline

    # ── Setup ─────────────────────────────────────────────────────────────────
    project_id = args.project.upper()
    output_dir = (Path(args.output) if args.output
                  else Path.home() / "TCGA_downloads" / project_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    cp = Checkpoint(output_dir)
    if args.fresh:
        cp.reset_from("auth")
        print("  🔄  Starting fresh (checkpoint cleared).")

    # ── Step 1: Client ────────────────────────────────────────────────────────
    client = GDCClient.from_file(args.token) if args.token else GDCClient()

    # ── Step 2: Validate project ──────────────────────────────────────────────
    print(f"\n  🔍  Validating project: {project_id}")
    info = client.validate_project(project_id)
    program_name = info.get("program", {}).get("name", "TCGA")
    print(f"  ✅  {info.get('name', project_id)}  (program: {program_name})")

    # ── Step 3: Discover files ────────────────────────────────────────────────
    if not cp.is_done("files_found"):
        print(f"\n  🔍  Discovering STAR-Counts files...")
        hits = client.discover_star_files(project_id)
        metadata_df = flatten_file_metadata(hits)
        from tcga_downloader.sample import summarise_sample_types
        summary = summarise_sample_types(metadata_df)
        print(f"  ✅  {len(hits)} files found.")
        print(f"      Tumor: {summary.get('Tumor', 0)}   "
              f"Normal: {summary.get('Normal', 0)}")
        cp.save("files_found", {"n_files": len(hits), "hits": hits})
    else:
        print("  ✅  File discovery: already done.")
        saved = cp.get("files_found")
        hits = saved.get("hits", [])
        metadata_df = flatten_file_metadata(hits)

    if args.dry_run:
        print("\n  ℹ️  Dry run complete. No files downloaded.")
        return

    # ── Step 4: Download ──────────────────────────────────────────────────────
    counts_dir = output_dir / "raw_counts"
    if not cp.is_done("downloaded"):
        counts_dir.mkdir(exist_ok=True)
        zip_path   = output_dir / "gdc_download.zip"
        file_ids   = [f["file_id"] for f in hits]
        print(f"\n  📦  Downloading {len(file_ids)} files...")
        client.stream_download(file_ids, str(zip_path))
        print("  📂  Extracting archive...")
        # GDC returns ZIP for multi-file downloads but occasionally tar.gz
        # Detect actual format by reading magic bytes — never assume extension
        with open(zip_path, "rb") as fh:
            magic = fh.read(4)

        if magic[:2] == b'\x1f\x8b':
            # gzip / tar.gz
            import tarfile
            archive_path = zip_path.with_suffix(".tar.gz")
            zip_path.rename(archive_path)
            with tarfile.open(archive_path, "r:gz") as t:
                t.extractall(counts_dir)
            archive_path.unlink(missing_ok=True)
        elif magic[:4] == b'PK\x03\x04':
            # standard ZIP
            with zipfile.ZipFile(zip_path, "r") as z:
                z.extractall(counts_dir)
            zip_path.unlink(missing_ok=True)
        else:
            file_size = zip_path.stat().st_size
            raise GDCError(
                f"Downloaded file is not a ZIP or tar.gz "
                f"(size: {file_size/1024:.1f} KB, magic: {magic.hex()}).",
                fix=(
                    "The GDC may have returned an error page instead of data.\n"
                    "Delete the partial file and retry:\n"
                    f"    rm {zip_path}\n"
                    "    tcga-download --project TCGA-CHOL --output ~/TCGA_test\n"
                    "(Completed steps are checkpointed and will not re-run.)"
                ),
                step="extraction",
            )
        cp.save("downloaded", {"counts_dir": str(counts_dir)})
    else:
        print("  ✅  Download: already done.")

    # ── Step 5: Clinical data ─────────────────────────────────────────────────
    clinical_cache = output_dir / "clinical_data.tsv"
    if not cp.is_done("clinical"):
        print("\n  🏥  Fetching clinical data...")
        case_hits    = client.fetch_clinical_data(project_id)
        clinical_df  = flatten_clinical(case_hits)
        clinical_df.to_csv(clinical_cache, sep="\t", index=False)
        cp.save("clinical", {"n_cases": len(clinical_df)})
        print(f"  ✅  {len(clinical_df)} clinical records fetched.")
    else:
        print("  ✅  Clinical data: already done.")
        clinical_df = pd.read_csv(clinical_cache, sep="\t", dtype=str)

    # ── Step 6: CDR annotations ───────────────────────────────────────────────
    metadata_with_cdr = metadata_df.copy()   # will be enriched if CDR runs

    cdr_cache = output_dir / "metadata_with_cdr.tsv"
    if not args.no_cdr and not cp.is_done("cdr"):
        print(f"\n  📚  Fetching PanCanAtlas CDR annotations...")
        try:
            cdr_result = run_cdr_pipeline(
                client       = client,
                metadata_df  = metadata_df,
                project_id   = project_id,
                program_name = program_name,
                output_dir   = output_dir,
                verbose      = True,
            )
            metadata_with_cdr = cdr_result["merged_df"]
            metadata_with_cdr.to_csv(cdr_cache, sep="\t", index=False)
            cp.save("cdr", {
                "n_matched":   int(metadata_with_cdr.get("cdr_matched", pd.Series()).sum()),
                "n_total":     len(metadata_with_cdr),
                "cache":       str(cdr_cache),
                "output_paths": {k: str(v) for k, v in cdr_result["output_paths"].items()},
            })
        except GDCError as e:
            # Non-TCGA project or other CDR error — warn but continue without CDR
            print(f"\n  ⚠️  CDR step skipped: {e.message}")
            print(f"      Pipeline will continue without CDR annotations.")
            cp.save("cdr", {"skipped": True, "reason": e.message})

    elif cp.is_done("cdr"):
        saved = cp.get("cdr")
        if saved.get("skipped"):
            print(f"  ⚠️  CDR: previously skipped ({saved.get('reason', 'unknown reason')})")
        elif cdr_cache.exists():
            print(f"  ✅  CDR: already done "
                  f"({saved.get('n_matched', '?')}/{saved.get('n_total', '?')} matched).")
            metadata_with_cdr = pd.read_csv(cdr_cache, sep="\t", dtype=str)
        else:
            print("  ✅  CDR: already done.")
    else:
        print("  ⏭   CDR step skipped (--no-cdr flag).")

    # ── Step 7: Count matrix ──────────────────────────────────────────────────
    matrix_cache = output_dir / "count_matrix_cache.tsv"
    if not cp.is_done("matrix"):
        print("\n  🔧  Building count matrix...")
        matrix = build_count_matrix(counts_dir, metadata_df, verbose=True)
        matrix.to_csv(matrix_cache, sep="\t")
        cp.save("matrix", {"cache": str(matrix_cache),
                            "shape": list(matrix.shape)})
    else:
        saved = cp.get("matrix")
        shape = saved.get("shape", ["?", "?"])
        print(f"  ✅  Count matrix: already built ({shape[0]} genes × {shape[1]} samples).")
        matrix = pd.read_csv(matrix_cache, sep="\t", index_col=0)

    # ── Step 8: Merge and save ────────────────────────────────────────────────
    if not cp.is_done("merged"):
        print("\n  🗂️   Merging and saving output files...")

        # Standard output files (counts + GDC clinical, no CDR)
        paths = merge_outputs(matrix, metadata_df, clinical_df,
                              project_id, output_dir, verbose=True)

        # CDR-enriched output files (if CDR ran successfully)
        has_cdr = any(c.startswith("cdr_") for c in metadata_with_cdr.columns)
        if has_cdr:
            cdr_paths = save_full_merged_with_cdr(
                count_matrix     = matrix,
                metadata_with_cdr = metadata_with_cdr,
                project_id       = project_id,
                output_dir       = output_dir,
                verbose          = True,
            )
            paths.update(cdr_paths)

        cp.save("merged", {k: str(v) for k, v in paths.items()})
        print(f"\n  🎉  Done! All files saved to: {output_dir}")
        print(f"\n  Output files:")
        for key, path in paths.items():
            p = Path(path)
            if p.exists():
                size_mb = p.stat().st_size / (1024**2)
                print(f"    • {p.name}  ({size_mb:.1f} MB)")
    else:
        print("  ✅  Already merged. All done.")
        print(f"  📂  Output directory: {output_dir}")


def main():
    """Entry point registered by pyproject.toml as 'tcga-download'."""
    parser = _build_parser()
    args   = parser.parse_args()

    try:
        if args.gui:
            run_gui()
        elif args.project:
            run_cli(args)
        else:
            # No arguments — print help
            parser.print_help()

    except GDCError as e:
        print("\n" + e.formatted())
        sys.exit(1)

    except KeyboardInterrupt:
        print("\n\n  ⛔  Cancelled. Progress is saved — re-run to resume.")
        sys.exit(0)

    except Exception as e:
        print(f"\n  ❌  Unexpected error: {type(e).__name__}: {e}")
        print("\n  Technical detail:")
        traceback.print_exc()
        print("\n  Please report this at: "
              "https://github.com/yourusername/tcga-gdc-downloader/issues")
        sys.exit(1)
