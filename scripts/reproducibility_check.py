"""
Compare key pipeline outputs between two runs (e.g. outputs_run1 vs outputs_run2).

Usage:
    python scripts/reproducibility_check.py outputs_run1 outputs_run2

Does not compare PNG files by default (binary metadata may differ).
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path


REL_FILES = [
    Path("tables") / "pathway_scores_evidence_filtered.csv",
    Path("tables") / "main_candidate_signals.csv",
    Path("tables") / "top_genes_from_pathway_level_signals.csv",
    Path("tables") / "top_damage_response_genes.csv",
    Path("logs") / "initial_analysis_summary.json",
    Path("logs") / "analysis_narrative.md",
]


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare pipeline outputs between two output directories."
    )
    parser.add_argument("dir_a", type=Path, help="First outputs folder (e.g. outputs_run1)")
    parser.add_argument("dir_b", type=Path, help="Second outputs folder (e.g. outputs_run2)")
    args = parser.parse_args()
    a = args.dir_a.resolve()
    b = args.dir_b.resolve()

    print(f"Comparing:\n  A: {a}\n  B: {b}\n")
    for rel in REL_FILES:
        pa, pb = a / rel, b / rel
        if not pa.is_file() and not pb.is_file():
            status = "MISSING (both)"
        elif not pa.is_file():
            status = "MISSING (A)"
        elif not pb.is_file():
            status = "MISSING (B)"
        else:
            sa, sb = file_sha256(pa), file_sha256(pb)
            status = "MATCH" if sa == sb else "DIFFERENT"
        print(f"  {rel}: {status}")

    print(
        "\nNote: CSV/JSON/Markdown should match for deterministic re-runs with the same inputs. "
        "PNG hashes are not compared here (metadata may differ)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
