"""
Extract NASA OSDR download URLs from a Torchlight Hackathon .ipynb and fetch files.

URLs are parsed from notebook cell sources — the same strings as in the official Colab.
Only gene-expression targets used by this pipeline are downloaded by default.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, unquote, urlparse
from urllib.request import Request, urlopen

# Must match load_data.py (official starter notebook filenames)
OSD569_FILENAME = "GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx"
OSD570_FILENAME = "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx"

DEFAULT_NOTEBOOK_NAMES = (
    "Torchlight_Hackathon_2026.ipynb",
    "Torchlight Hackathon 2026.ipynb",
)

# OSDR file downloads from the starter notebook / geode-py
OSDR_URL_RE = re.compile(
    r"https://osdr\.nasa\.gov/geode-py/ws/studies/OSD-\d+/download\?[^\s'\"\)\]\}]+"
)


def iter_notebook_text(notebook_path: Path) -> Iterable[str]:
    """Yield concatenated source text from code and markdown cells."""
    raw = notebook_path.read_text(encoding="utf-8")
    nb = json.loads(raw)
    for cell in nb.get("cells", []):
        src = cell.get("source", [])
        if isinstance(src, list):
            yield "".join(src)
        else:
            yield str(src)


def extract_osdr_urls(notebook_path: Path) -> list[str]:
    """Return unique OSDR download URLs found in the notebook."""
    text = "\n".join(iter_notebook_text(notebook_path))
    found = []
    for m in OSDR_URL_RE.finditer(text):
        url = m.group(0).rstrip("\\")
        # Trim trailing punctuation sometimes captured from markdown
        while url.endswith((")", ",", ".", ";")):
            url = url[:-1]
        if url not in found:
            found.append(url)
    return found


def filename_from_osdr_url(url: str) -> str | None:
    """Return the `file=` query value if present."""
    q = urlparse(url).query
    params = parse_qs(q)
    vals = params.get("file")
    if not vals:
        return None
    return unquote(vals[0])


def urls_for_target_files(urls: list[str], targets: set[str]) -> dict[str, str]:
    """
    Map official filename -> full URL for each target that appears in the URL list.
    If multiple URLs map to the same file, the first wins.
    """
    out: dict[str, str] = {}
    for url in urls:
        fn = filename_from_osdr_url(url)
        if fn and fn in targets and fn not in out:
            out[fn] = url
    return out


def find_torchlight_notebook(root: Path) -> Path | None:
    """
    Locate the starter notebook under data/raw/ only.

    Expected: data/raw/Torchlight_Hackathon_2026.ipynb (or a matching *torchlight*hackathon*.ipynb).
    """
    base = root / "data" / "raw"
    if not base.is_dir():
        return None
    for name in DEFAULT_NOTEBOOK_NAMES:
        p = base / name
        if p.is_file():
            return p
    for p in base.glob("*.ipynb"):
        n = p.name.lower()
        if "torchlight" in n and "hackathon" in n:
            return p
    return None


def download_url_to(url: str, dest: Path, timeout: int = 120) -> None:
    """GET url and write to dest (streaming)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = Request(url, headers={"User-Agent": "spaceflight-damage-response/1.0"})
    with urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    dest.write_bytes(data)


def fetch_gene_expression_tables(
    root: Path | None = None,
    notebook_path: Path | None = None,
    dest_dir: Path | None = None,
) -> dict[str, str]:
    """
    Download OSD-569 / OSD-570 gene expression XLSX files referenced in the notebook.

    Returns a dict filename -> status message ("saved", or error text).
    """
    root = root or Path(__file__).resolve().parent.parent
    dest_dir = dest_dir or (root / "data" / "raw")

    if notebook_path is not None:
        nb = Path(notebook_path).expanduser().resolve()
        if not nb.is_file():
            return {
                "_error": f"Notebook not found: {nb}",
            }
    else:
        nb = find_torchlight_notebook(root)
        if nb is None:
            return {
                "_error": (
                    "No Torchlight notebook found in data/raw/. "
                    "Place Torchlight_Hackathon_2026.ipynb there (Colab export)."
                ),
            }

    targets = {OSD569_FILENAME, OSD570_FILENAME}
    urls = extract_osdr_urls(nb)
    mapping = urls_for_target_files(urls, targets)

    results: dict[str, str] = {}
    for fn in targets:
        url = mapping.get(fn)
        if not url:
            results[fn] = (
                f"not found in {nb.name} — URL for this file may be missing from the notebook export."
            )
            continue
        out_path = dest_dir / fn
        try:
            download_url_to(url, out_path)
            results[fn] = f"saved -> {out_path}"
        except HTTPError as exc:
            results[fn] = f"HTTP error {exc.code}: {exc.reason}"
        except URLError as exc:
            results[fn] = f"network error: {exc.reason}"
        except OSError as exc:
            results[fn] = f"write error: {exc}"

    results["_notebook"] = str(nb)
    results["_urls_found"] = str(len(urls))
    return results
