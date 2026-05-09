"""
Download OSD-569 / OSD-570 tables from URLs in data/raw/Torchlight_Hackathon_2026.ipynb,
then run the full pipeline. Writes a timestamped log under outputs/logs/.

Usage (repo root, venv activated):
    python scripts/run_real_data_pipeline.py
"""

from __future__ import annotations

import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

LOG_DIR = ROOT / "outputs" / "logs"


def main() -> int:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"real_data_run_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%SZ')}.log"

    lines: list[str] = []

    def log(msg: str) -> None:
        lines.append(msg)
        print(msg)

    log(f"=== Real-data pipeline run {datetime.now(timezone.utc).isoformat()} ===")

    try:
        from osdr_notebook_fetch import fetch_gene_expression_tables

        log("\n-- Step 1: fetch from notebook URLs --")
        fetch_results = fetch_gene_expression_tables(ROOT)
        for key in sorted(fetch_results.keys()):
            log(f"  {key}: {fetch_results[key]}")

        raw_dir = ROOT / "data" / "raw"
        for name in (
            "GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx",
            "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx",
        ):
            p = raw_dir / name
            if p.is_file():
                log(f"  [OK] on disk: {name} ({p.stat().st_size} bytes)")
            else:
                log(f"  [MISSING] {name}")

        if "_error" in fetch_results:
            log("\nFetch reported an error; fix it before expecting real DE outputs.")
            log_path.write_text("\n".join(lines), encoding="utf-8")
            return 1

        log("\n-- Step 2: main pipeline --")
        import main as pipeline_main

        code = pipeline_main.main()
        log(f"\nmain.py exited with {code}")
        log_path.write_text("\n".join(lines), encoding="utf-8")
        return int(code)

    except Exception as exc:
        log(f"\nFATAL: {exc}")
        log(traceback.format_exc())
        log_path.write_text("\n".join(lines), encoding="utf-8")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
