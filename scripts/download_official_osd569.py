"""
Download the official NASA OSD-569 / GLDS-561 processed gene-expression workbook.

Optional: run manually when you need a fresh copy. main.py does not require network.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent.parent
DEST_REL = Path("data/raw/GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx")
OFFICIAL_URL = (
    "https://osdr.nasa.gov/geode-py/ws/studies/OSD-569/download?"
    "source=datamanager&file=GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx"
)


def download(url: str, dest: Path, timeout: int = 180) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = Request(url, headers={"User-Agent": "spaceflight-damage-response/download_official_osd569"})
    with urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    dest.write_bytes(data)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download official OSD-569 GLDS-561 Excel from OSDR.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite data/raw file if it already exists.",
    )
    args = parser.parse_args()

    dest = ROOT / DEST_REL
    log_dir = ROOT / "outputs" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "osd569_download_report.txt"

    lines = [
        f"osd569_download_report generated at {datetime.now(timezone.utc).isoformat()}",
        f"url: {OFFICIAL_URL}",
        f"destination: {dest}",
        "",
    ]

    if dest.is_file() and not args.force:
        sz = dest.stat().st_size
        msg = f"File already exists ({sz} bytes). Use --force to overwrite."
        lines.append(msg)
        log_path.write_text("\n".join(lines), encoding="utf-8")
        print(msg)
        print(f"Log: {log_path}")
        return 0

    try:
        download(OFFICIAL_URL, dest)
    except HTTPError as exc:
        err = f"HTTP {exc.code}: {exc.reason}"
        lines.append(f"ERROR: {err}")
        log_path.write_text("\n".join(lines), encoding="utf-8")
        print(err, file=sys.stderr)
        return 1
    except URLError as exc:
        err = f"Network error: {exc.reason}"
        lines.append(f"ERROR: {err}")
        log_path.write_text("\n".join(lines), encoding="utf-8")
        print(err, file=sys.stderr)
        return 1

    sz = dest.stat().st_size
    lines.append(f"status: saved")
    lines.append(f"file_size_bytes: {sz}")
    log_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Downloaded {sz} bytes -> {dest}")
    print(f"Log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
