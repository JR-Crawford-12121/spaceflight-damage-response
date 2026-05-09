"""
Create a starter Ensembl → HGNC mapping file for OSD-569 whole-blood RNA.

The main pipeline does not download mappings from the network. Fill in rows
offline (e.g. from Ensembl BioMart, HGNC, or study metadata) and save as:

  data/processed/ensembl_to_symbol.tsv

Then re-run: python main.py
"""

from __future__ import annotations

from pathlib import Path

_TEMPLATE_LINES = """# OSD-569 uses Ensembl gene IDs (e.g. ENSG00000000003.15).
# Curated pathway gene sets in this project use HGNC symbols (e.g. TP53).
# Provide a two-column table mapping Ensembl stable IDs to HGNC symbols.
# Version suffixes on Ensembl IDs (.15, etc.) are stripped automatically when matching.
#
# Save a completed copy as: data/processed/ensembl_to_symbol.tsv
# See docs/ENSEMBL_MAPPING.md for full instructions.
#
ensembl_gene_id\thgnc_symbol
ENSG00000141510\tTP53
ENSG00000139618\tBRCA2
"""


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "data" / "processed"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "ensembl_to_symbol_TEMPLATE.tsv"
    path.write_text(_TEMPLATE_LINES, encoding="utf-8")
    print(f"Wrote {path}")
    print("Copy/rename to ensembl_to_symbol.tsv and add your mappings, then run python main.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
