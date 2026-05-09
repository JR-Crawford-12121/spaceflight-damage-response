"""
Build data/processed/ensembl_to_symbol.tsv for OSD-569 Ensembl gene IDs.

Priority (full mode):
  A) Workbook symbol column (offline).
  B) Ensembl REST batch (requests + certifi).
  C) MyGene.info ensemblgene → symbol batch.

Curated-only mode (--curated-only):
  Map only genes in src/gene_sets.py pathway lists → Ensembl IDs present in OSD-569 (small TSV).

Usage:
    python scripts/build_ensembl_to_symbol_mapping.py
    python scripts/build_ensembl_to_symbol_mapping.py --curated-only
    python scripts/build_ensembl_to_symbol_mapping.py --skip-network
    python scripts/build_ensembl_to_symbol_mapping.py --curated-only --insecure-ssl

main.py does not call this script.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from ensembl_mapping_http import (  # noqa: E402
    lookup_ensembl_rest_requests,
    lookup_mygene_ensembl_to_symbol,
    lookup_mygene_symbols_to_ensembl_in_osd569,
    ssl_verify_bundle,
)
from ensembl_mapping_sources import extract_mapping_from_wide_table  # noqa: E402
from gene_id_mapping import strip_ensembl_version  # noqa: E402
from gene_sets import union_pathway_gene_symbols  # noqa: E402
from load_data import OSD569_FILENAME, locate_file  # noqa: E402
from osd569_io import read_osd569_processed_wide, unique_ensembl_ids_from_wide  # noqa: E402


def _normalize_mapping(raw: dict[str, str]) -> dict[str, str]:
    cleaned: dict[str, str] = {}
    for k, v in raw.items():
        kk = strip_ensembl_version(str(k).strip())
        if kk.startswith("ENSG") and str(v).strip():
            cleaned[kk] = str(v).strip()
    return cleaned


def main() -> int:
    parser = argparse.ArgumentParser(description="Build Ensembl→HGNC TSV for OSD-569.")
    parser.add_argument(
        "--skip-network",
        action="store_true",
        help="Only workbook-derived mappings (no REST/MyGene).",
    )
    parser.add_argument(
        "--xlsx",
        type=Path,
        default=None,
        help="Override path to GLDS-561 XLSX.",
    )
    parser.add_argument(
        "--insecure-ssl",
        action="store_true",
        help="Disable SSL verification for mapping HTTP requests (local troubleshooting only).",
    )
    parser.add_argument(
        "--curated-only",
        action="store_true",
        help=(
            "Only map curated pathway gene symbols (gene_sets.py) to Ensembl IDs present "
            "in OSD-569 — small mapping, not genome-wide."
        ),
    )
    args = parser.parse_args()

    report_path = ROOT / "outputs" / "logs" / "ensembl_mapping_build_report.txt"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    map_out = ROOT / "data" / "processed" / "ensembl_to_symbol.tsv"
    unmapped_out = ROOT / "outputs" / "tables" / "unmapped_osd569_ensembl_ids.csv"

    verify = ssl_verify_bundle(args.insecure_ssl)
    methods_used: list[str] = []

    lines: list[str] = []
    lines.append("OSD-569 Ensembl → HGNC mapping build")
    lines.append(f"repo root: {ROOT}")
    lines.append(f"SSL verify: {verify} (certifi bundle)" if verify else "SSL verify: False (--insecure-ssl)")
    lines.append(f"--insecure-ssl: {args.insecure_ssl}")
    lines.append(f"--curated-only: {args.curated_only}")
    lines.append(f"--skip-network: {args.skip_network}")
    lines.append("")

    if args.insecure_ssl:
        warn = (
            "WARNING: SSL verification disabled for mapping requests. "
            "Use only for local troubleshooting."
        )
        lines.append(warn)
        print(warn)
        try:
            import urllib3

            urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        except Exception:
            pass
        lines.append("")

    path = args.xlsx if args.xlsx else locate_file(OSD569_FILENAME, ROOT)
    if path is not None:
        path = Path(path).resolve()

    if path is None or not path.is_file():
        msg = f"ERROR: OSD-569 XLSX not found ({OSD569_FILENAME}). Place it under data/raw/."
        lines.append(msg)
        report_path.write_text("\n".join(lines), encoding="utf-8")
        print(msg)
        return 1

    lines.append(f"Input XLSX: {path}")
    try:
        wide = read_osd569_processed_wide(path)
    except Exception as exc:
        lines.append(f"ERROR reading processed wide table: {exc}")
        lines.append(traceback.format_exc())
        report_path.write_text("\n".join(lines), encoding="utf-8")
        print(lines[-2])
        return 1

    all_ids_list = unique_ensembl_ids_from_wide(wide)
    all_ids_set = set(all_ids_list)
    lines.append(f"total unique Ensembl gene IDs (stripped index): {len(all_ids_list)}")
    lines.append("")

    workbook_map = extract_mapping_from_wide_table(wide)
    merged: dict[str, str] = dict(workbook_map)
    if workbook_map:
        methods_used.append("workbook_columns")

    lines.append(f"mappings from workbook symbol column: {len(workbook_map)}")
    if workbook_map:
        for k, v in list(sorted(workbook_map.items()))[:20]:
            lines.append(f"  example: {k} -> {v}")

    network_errors: list[str] = []

    if args.curated_only:
        lines.append("")
        lines.append(
            "Curated-only mapping mode was used. This mapping is not genome-wide; "
            "it only supports the current pathway gene sets in gene_sets.py."
        )
        curated_syms = union_pathway_gene_symbols()
        lines.append(f"unique curated pathway symbols: {len(curated_syms)}")
        if not args.skip_network:
            cm, err = lookup_mygene_symbols_to_ensembl_in_osd569(
                curated_syms,
                all_ids_set,
                verify=verify,
            )
            merged.update(cm)
            if cm:
                methods_used.append("mygene_symbol_to_ensembl_curated")
            if err:
                network_errors.append(f"MyGene curated symbol lookup: {err}")
                lines.append(f"MyGene (curated) note: {err}")
        else:
            lines.append("(skip-network: no MyGene symbol lookup attempted)")
    elif not args.skip_network:
        need_rest = [e for e in all_ids_list if e not in merged]
        lines.append("")
        lines.append(f"IDs needing symbols after workbook (full mode): {len(need_rest)}")

        if need_rest:
            rest_map, rest_err = lookup_ensembl_rest_requests(
                need_rest,
                verify=verify,
            )
            if rest_map:
                merged.update(rest_map)
                methods_used.append("ensembl_rest")
            if rest_err:
                network_errors.append(f"Ensembl REST: {rest_err}")
                lines.append(f"Ensembl REST error (may be partial): {rest_err}")

        need_mg = [e for e in all_ids_list if e not in merged]
        lines.append(f"IDs still needing symbols before MyGene: {len(need_mg)}")
        if need_mg:
            mg_map, mg_err = lookup_mygene_ensembl_to_symbol(
                need_mg,
                verify=verify,
            )
            if mg_map:
                merged.update(mg_map)
                methods_used.append("mygene_ensembl_to_symbol")
            if mg_err:
                network_errors.append(f"MyGene ensembl→symbol: {mg_err}")
                lines.append(f"MyGene error (may be partial): {mg_err}")
    else:
        lines.append("")
        lines.append("(skip-network: no Ensembl REST / MyGene attempted)")

    cleaned = _normalize_mapping(merged)
    unmapped = sorted(all_ids_set - set(cleaned.keys()))

    lines.append("")
    lines.append("--- Summary ---")
    lines.append(f"source methods: {', '.join(methods_used) if methods_used else '(none)'}")
    lines.append(f"mapped unique IDs: {len(cleaned)}")
    lines.append(f"unmapped unique IDs: {len(unmapped)}")
    pct = (100.0 * len(cleaned) / len(all_ids_list)) if all_ids_list else 0.0
    lines.append(f"percent of OSD-569 IDs mapped: {pct:.4f}%")
    lines.append("")
    lines.append("First 20 mappings (sorted by Ensembl ID):")
    for k, v in list(sorted(cleaned.items()))[:20]:
        lines.append(f"  {k} -> {v}")

    if network_errors:
        lines.append("")
        lines.append("Network / API messages:")
        for ne in network_errors:
            lines.append(f"  - {ne}")

    lines.append("")
    lines.append(f"output path (if written): {map_out}")

    unmapped_out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"ensembl_gene_id": unmapped}).to_csv(unmapped_out, index=False)
    lines.append(f"wrote unmapped IDs list: {unmapped_out}")

    pre_existing = map_out.is_file()
    pre_size = map_out.stat().st_size if pre_existing else 0

    if len(cleaned) == 0:
        lines.append("")
        lines.append("FAILURE: No mappings were produced.")
        lines.append("Existing ensembl_to_symbol.tsv was NOT overwritten.")
        if pre_existing:
            lines.append(f"Previous file kept ({pre_size} bytes).")
        report_path.write_text("\n".join(lines), encoding="utf-8")

        err_msg = (
            "No mappings were produced. Existing ensembl_to_symbol.tsv was not overwritten."
        )
        print(err_msg)
        print(f"\nFull report: {report_path}")
        return 1

    map_out.parent.mkdir(parents=True, exist_ok=True)
    df_out = pd.DataFrame(
        [{"ensembl_gene_id": k, "hgnc_symbol": v} for k, v in sorted(cleaned.items())]
    )
    df_out.to_csv(map_out, sep="\t", index=False)
    lines.append("")
    lines.append(f"Wrote mapping file: {map_out} ({len(cleaned)} rows)")

    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Mapped {len(cleaned)} Ensembl IDs — wrote {map_out}")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
