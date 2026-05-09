"""
Download OSD-569 / OSD-570 gene-expression XLSX files using URLs embedded in the Torchlight notebook.

Usage:
  python scripts/fetch_osdr_from_notebook.py
  python scripts/fetch_osdr_from_notebook.py "C:\\path\\to\\Torchlight_Hackathon_2026.ipynb"

Put **`Torchlight_Hackathon_2026.ipynb`** in **`data/raw/`**, or pass an explicit path as the first argument.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from osdr_notebook_fetch import fetch_gene_expression_tables  # noqa: E402


def main() -> int:
    nb = Path(sys.argv[1]).expanduser() if len(sys.argv) > 1 else None
    results = fetch_gene_expression_tables(ROOT, notebook_path=nb)
    print("OSDR fetch results:")
    for key in sorted(results.keys()):
        print(f"  {key}: {results[key]}")
    return 1 if "_error" in results else 0


if __name__ == "__main__":
    raise SystemExit(main())
