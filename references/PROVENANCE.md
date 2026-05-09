# Reference material

Keep the official starter notebook export as **`data/raw/Torchlight_Hackathon_2026.ipynb`** (File → Download → .ipynb from Google Colab, or copy from your machine into **`data/raw/`**).

The notebook defines **many** OSDR Excel downloads (RNA-seq, m6A, ATAC, VDJ, spatial, other studies). **This repository’s pathway RNA pipeline only ingests two gene-expression tables** from that notebook; see the full inventory and contrasts in **`docs/NOTEBOOK_DATA_SOURCES.md`** and kwargs in **`src/notebook_manifest.py`**.

The pipeline can read OSDR URLs from the notebook (`python main.py --download-from-notebook`). Canonical **`pd.read_excel`** patterns for the **two loaded files** are implemented in **`src/load_data.py`** and match the notebook JSON for those cells.

OSDR base URL:

`https://osdr.nasa.gov/geode-py/ws/studies/OSD-{study}/download?source=datamanager&file={filename}`

**Loaded by `main.py` (gene-level pathway RNA scoring):**

- **OSD-569:** `GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx` — contrast **(R+82) vs (L-92, L-44, L-3)**
- **OSD-570:** `GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx` — contrast **(R+45) vs (R+1)**

**Also in the notebook but not loaded by this pipeline (different modality / row type):** e.g. OSD-569 m6A (`GLDS-561_directm6Aseq_*`), OSD-570 snATAC (`GLDS-562_snATAC-Seq_*`), VDJ, spatial — see **`docs/NOTEBOOK_DATA_SOURCES.md`**.
