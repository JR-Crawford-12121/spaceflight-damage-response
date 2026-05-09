# Ensembl → HGNC mapping for OSD-569

## Why this file exists

**OSD-569** whole-blood RNA differential-expression exports typically use **Ensembl gene IDs** (for example `ENSG00000000003.15`).

The curated DNA damage–response pathway lists in this project use **HGNC gene symbols** (for example `TP53`, `BRCA2`).

Until Ensembl IDs are mapped to symbols, **overlap between OSD-569 genes and those curated lists is usually zero**, even when the biology is present in the table.

## Expected mapping file

**Path:**

`data/processed/ensembl_to_symbol.tsv`

**Recommended columns:**

| Column | Purpose |
|--------|---------|
| `ensembl_gene_id` | Ensembl stable gene ID (version suffix optional) |
| `hgnc_symbol` | HGNC gene symbol |

The pipeline accepts several alternate header names; see `src/gene_id_mapping.py`.

## Version suffixes

Ensembl IDs may include a version suffix (for example `.15`). The pipeline **strips these automatically** when matching both the mapping file and the expression table (`ENSG00000000003.15` → `ENSG00000000003`).

## Main pipeline stays offline

**`python main.py` never downloads mapping data.** If `ensembl_to_symbol.tsv` is already present, scoring runs fully offline.

Optional helper scripts may use the network **only when you run them** — not when you run `main.py`.

---

## Option 1 — Inspect the official workbook first (recommended)

The GLDS-561 processed XLSX may already contain gene symbols alongside Ensembl IDs.

```bash
python scripts/inspect_osd569_workbook.py
```

This writes **`outputs/logs/osd569_workbook_inspection.txt`** (sheet names, columns, samples).

If the processed DE sheet contains a detectable symbol column, the script may **automatically generate** `data/processed/ensembl_to_symbol.tsv` (when enough pairs are found).

---

## Option 2 — Build mapping from the workbook + optional APIs

If the workbook does **not** expose symbols (or coverage is incomplete), run:

```bash
python scripts/build_ensembl_to_symbol_mapping.py
```

The builder uses **`requests`** + **`certifi`** for TLS verification by default (fixes many `CERTIFICATE_VERIFY_FAILED` issues on macOS/Python).

Behavior:

1. Loads the same OSD-569 wide table as **`src/load_data.py`**.
2. Extracts unique **`ENSG`** IDs from the gene index.
3. Tries **workbook symbol columns** first (offline).
4. Then **Ensembl REST** batch (`https://rest.ensembl.org/lookup/id`).
5. Then **MyGene.info** batch (`ensemblgene` → `symbol`).
6. Writes **`data/processed/ensembl_to_symbol.tsv`** only if **at least one** mapping row is produced (**never overwrites with an empty file**).
7. Writes **`outputs/logs/ensembl_mapping_build_report.txt`** and **`outputs/tables/unmapped_osd569_ensembl_ids.csv`**.

### Curated genes only (recommended when full genome mapping is slow or SSL is flaky)

Maps only symbols from **`src/gene_sets.py`** that appear in the OSD-569 gene list (small TSV, enough for pathway overlap):

```bash
python scripts/build_ensembl_to_symbol_mapping.py --curated-only
```

### SSL troubleshooting (local hackathons only)

```bash
python scripts/build_ensembl_to_symbol_mapping.py --curated-only --insecure-ssl
```

Offline-only (workbook extraction only):

```bash
python scripts/build_ensembl_to_symbol_mapping.py --skip-network
```

If every API fails, the script **keeps any existing** `ensembl_to_symbol.tsv` and exits with an error instead of writing an empty file.

---

## Validate the mapping file

```bash
python scripts/validate_ensembl_mapping.py
```

Writes **`outputs/logs/ensembl_mapping_validation.txt`** (columns, duplicates, missing symbols, presence of common DDR-related symbols in the mapping values).

---

## Run the science pipeline

```bash
python main.py
```

---

## Manual template (no automation)

Generate a starter TSV with comments:

```bash
python scripts/create_ensembl_mapping_template.py
```

Copy or rename to **`ensembl_to_symbol.tsv`**, fill rows offline (BioMart, HGNC, etc.), then run **`python main.py`**.

## Interpretation

A missing or incomplete mapping file does **not** prove that OSD-569 lacks DNA-damage-response–related genes; it only means this pipeline **cannot align** those rows to HGNC-based pathway lists until identifiers are mapped.
