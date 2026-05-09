"""
Validate data/processed/ensembl_to_symbol.tsv before running python main.py.

Usage:
    python scripts/validate_ensembl_mapping.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from gene_id_mapping import load_ensembl_symbol_mapping, strip_ensembl_version  # noqa: E402

PATHWAY_SYMBOL_CHECKS = [
    "TP53",
    "ATM",
    "ATR",
    "BRCA1",
    "BRCA2",
    "RAD51",
    "XPC",
    "DDB2",
    "GADD45A",
    "CDKN1A",
    "BAX",
    "SOD1",
    "SOD2",
    "GPX4",
    "RB1",
]


def main() -> int:
    root = ROOT
    mapping_path = root / "data" / "processed" / "ensembl_to_symbol.tsv"
    out_path = root / "outputs" / "logs" / "ensembl_mapping_validation.txt"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("Ensembl→HGNC mapping validation")
    lines.append(f"file: {mapping_path}")
    lines.append("")

    if not mapping_path.is_file():
        msg = "ERROR: mapping file not found."
        lines.append(msg)
        lines.append("Run: python scripts/build_ensembl_to_symbol_mapping.py")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        print(msg)
        return 1

    if mapping_path.stat().st_size == 0:
        msg = (
            "Mapping file exists but is empty.\n"
            "Rerun: python scripts/build_ensembl_to_symbol_mapping.py "
            "(try --curated-only or --insecure-ssl if SSL failed)."
        )
        lines.append(msg.replace("\n", " "))
        out_path.write_text("\n".join(lines), encoding="utf-8")
        print(msg)
        return 1

    raw_preview = mapping_path.read_text(encoding="utf-8", errors="replace").strip()
    if not raw_preview:
        msg = (
            "Mapping file exists but is empty.\n"
            "Rerun: python scripts/build_ensembl_to_symbol_mapping.py"
        )
        lines.append("Mapping file has no readable content.")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        print(msg)
        return 1

    try:
        df = pd.read_csv(mapping_path, sep="\t", dtype=str, comment="#")
    except pd.errors.EmptyDataError:
        msg = (
            "Mapping file could not be parsed (no columns).\n"
            "Rerun: python scripts/build_ensembl_to_symbol_mapping.py"
        )
        lines.append(msg.replace("\n", " "))
        out_path.write_text("\n".join(lines), encoding="utf-8")
        print(msg)
        return 1
    except Exception as exc:
        lines.append(f"ERROR reading TSV: {exc}")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        print(lines[-1])
        return 1

    if df.empty or len(df.columns) == 0:
        msg = (
            "Mapping file has no data rows.\n"
            "Rerun: python scripts/build_ensembl_to_symbol_mapping.py"
        )
        lines.append("ERROR: zero rows in mapping table.")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        print(msg)
        return 1

    lines.append(f"rows read from file: {len(df)}")
    ens_col = sym_col = None
    for cand in ["ensembl_gene_id", "ensembl_id", "gene_id"]:
        if cand in df.columns:
            ens_col = cand
            break
    for cand in ["hgnc_symbol", "gene_symbol", "symbol"]:
        if cand in df.columns:
            sym_col = cand
            break
    if ens_col is None or sym_col is None:
        lines.append(
            f"ERROR: could not detect Ensembl + symbol columns. Found columns: {list(df.columns)}"
        )
        out_path.write_text("\n".join(lines), encoding="utf-8")
        print(lines[-1])
        return 1

    lines.append(f"using columns: {ens_col} , {sym_col}")

    bad_sym = df[sym_col].isna() | (df[sym_col].astype(str).str.strip() == "")
    n_bad_sym = int(bad_sym.sum())
    lines.append(f"rows with missing symbol: {n_bad_sym}")

    keys = df[ens_col].map(lambda x: strip_ensembl_version(str(x).strip()))
    dup_mask = keys.duplicated(keep=False)
    n_dup = int(dup_mask.sum())
    lines.append(f"rows involved in duplicate Ensembl keys: {n_dup}")

    mapping, meta = load_ensembl_symbol_mapping(mapping_path)
    lines.append(f"dictionary pairs loaded (gene_id_mapping loader): {meta.get('n_pairs_loaded', 0)}")
    lines.append("")

    sym_upper = {str(v).strip().upper() for v in mapping.values()}
    lines.append("pathway reference symbols present in mapping values:")
    for g in PATHWAY_SYMBOL_CHECKS:
        ok = g.upper() in sym_upper
        lines.append(f"  {g}: {'yes' if ok else 'no'}")

    lines.append("")
    lines.append("OK — validation checks completed (review duplicates/missing above).")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
