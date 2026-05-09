"""
Internet-based Ensembl ↔ HGNC lookups for optional mapping scripts only.

Uses requests + certifi for SSL. main.py does not import this module.
"""

from __future__ import annotations

import json
from typing import Any

import certifi
import requests

from gene_id_mapping import strip_ensembl_version

USER_AGENT = (
    "spaceflight-damage-response/1.0 "
    "(mapping helper; https://github.com/)"
)

ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id"
MYGENE_QUERY_URL = "https://mygene.info/v3/query"


def ssl_verify_bundle(insecure_ssl: bool) -> bool | str:
    """Return certifi bundle path, or False if explicitly insecure."""
    if insecure_ssl:
        return False
    return certifi.where()


def _chunks(xs: list[str], n: int) -> list[list[str]]:
    return [xs[i : i + n] for i in range(0, len(xs), n)]


def lookup_ensembl_rest_requests(
    ensembl_ids: list[str],
    *,
    verify: bool | str,
    chunk_size: int = 400,
    timeout: int = 180,
) -> tuple[dict[str, str], str | None]:
    """
    Batch POST https://rest.ensembl.org/lookup/id
    Returns (mapping ensembl -> display_name/symbol, error_message or None).
    """
    mapping: dict[str, str] = {}
    last_err: str | None = None
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": USER_AGENT,
    }

    for chunk in _chunks(ensembl_ids, chunk_size):
        try:
            r = requests.post(
                ENSEMBL_LOOKUP_URL,
                json={"ids": chunk},
                headers=headers,
                verify=verify,
                timeout=timeout,
            )
            if r.status_code >= 400:
                last_err = f"HTTP {r.status_code}: {r.text[:400]}"
                continue
            payload: Any = r.json()
        except requests.RequestException as exc:
            last_err = f"{type(exc).__name__}: {exc}"
            continue

        if isinstance(payload, list):
            for j, item in enumerate(payload):
                if not isinstance(item, dict) or item.get("error"):
                    continue
                rid = item.get("id") or (chunk[j] if j < len(chunk) else None)
                if rid is None:
                    continue
                key = strip_ensembl_version(str(rid))
                sym = item.get("display_name") or item.get("display_label")
                if sym:
                    mapping[key] = str(sym).strip()
        elif isinstance(payload, dict):
            for k, item in payload.items():
                if isinstance(item, dict):
                    key = strip_ensembl_version(str(item.get("id", k)))
                    sym = item.get("display_name") or item.get("display_label")
                    if sym:
                        mapping[key] = str(sym).strip()

    return mapping, last_err


def _parse_mygene_hits(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [x for x in payload if isinstance(x, dict)]
    if isinstance(payload, dict):
        if "hits" in payload and isinstance(payload["hits"], list):
            return [x for x in payload["hits"] if isinstance(x, dict)]
        if "results" in payload and isinstance(payload["results"], list):
            return [x for x in payload["results"] if isinstance(x, dict)]
    return []


def _ensembl_from_hit(hit: dict[str, Any]) -> str | None:
    eg = hit.get("ensemblgene")
    if isinstance(eg, dict):
        eg = eg.get("gene") or eg.get("id")
    emb = hit.get("ensembl")
    if eg is None and isinstance(emb, dict):
        eg = emb.get("gene") or emb.get("id")
    if eg is None and isinstance(emb, str):
        eg = emb
    if eg is None:
        return None
    s = strip_ensembl_version(str(eg).strip())
    return s if s.startswith("ENSG") else None


def lookup_mygene_ensembl_to_symbol(
    ensembl_ids: list[str],
    *,
    verify: bool | str,
    chunk_size: int = 1000,
    timeout: int = 180,
) -> tuple[dict[str, str], str | None]:
    """
    Map Ensembl gene IDs -> HGNC symbol via MyGene.info (batch JSON POST).
    """
    mapping: dict[str, str] = {}
    last_err: str | None = None
    headers = {"Content-Type": "application/json", "User-Agent": USER_AGENT}

    for chunk in _chunks(ensembl_ids, chunk_size):
        body = {
            "q": chunk,
            "scopes": "ensemblgene",
            "fields": "symbol",
            "species": "human",
        }
        try:
            r = requests.post(
                MYGENE_QUERY_URL,
                json=body,
                headers=headers,
                verify=verify,
                timeout=timeout,
            )
            if r.status_code >= 400:
                last_err = f"MyGene HTTP {r.status_code}: {r.text[:400]}"
                continue
            data = r.json()
        except requests.RequestException as exc:
            last_err = f"MyGene {type(exc).__name__}: {exc}"
            continue

        rows: list[Any] = []
        if isinstance(data, list):
            rows = data
        elif isinstance(data, dict):
            rows = _parse_mygene_hits(data)

        got = 0
        for j, item in enumerate(rows):
            if not isinstance(item, dict):
                continue
            if item.get("notfound"):
                continue
            sym = item.get("symbol")
            if not sym:
                continue
            rid = item.get("query") or (chunk[j] if j < len(chunk) else None)
            key = strip_ensembl_version(str(rid)) if rid else None
            if key and key.startswith("ENSG"):
                mapping[key] = str(sym).strip()
                got += 1
                continue
            eg = _ensembl_from_hit(item)
            if eg and eg.startswith("ENSG"):
                mapping[eg] = str(sym).strip()
                got += 1

        if got == 0 and isinstance(data, dict) and not rows:
            try:
                r2 = requests.post(
                    MYGENE_QUERY_URL,
                    data={
                        "q": ",".join(chunk),
                        "scopes": "ensemblgene",
                        "fields": "symbol",
                        "species": "human",
                    },
                    headers={"User-Agent": USER_AGENT},
                    verify=verify,
                    timeout=timeout,
                )
                if r2.ok:
                    data2 = r2.json()
                    for hit in _parse_mygene_hits(data2):
                        sym = hit.get("symbol")
                        eg = _ensembl_from_hit(hit)
                        if eg and sym:
                            mapping[eg] = str(sym).strip()
            except requests.RequestException as exc:
                last_err = f"MyGene form fallback {type(exc).__name__}: {exc}"

    return mapping, last_err


def lookup_mygene_symbols_to_ensembl_in_osd569(
    symbols_upper: set[str],
    osd569_ensembl_ids: set[str],
    *,
    verify: bool | str,
    chunk_size: int = 500,
    timeout: int = 180,
) -> tuple[dict[str, str], str | None]:
    """
    Resolve HGNC symbols -> Ensembl gene IDs, keep pairs where Ensembl ID is in OSD-569 DE table.
    Returns mapping ensembl_id -> symbol (as in pipeline TSV).
    """
    mapping: dict[str, str] = {}
    last_err: str | None = None
    sym_list = sorted(symbols_upper)
    headers = {"Content-Type": "application/json", "User-Agent": USER_AGENT}

    for chunk in _chunks(sym_list, chunk_size):
        body = {
            "q": chunk,
            "scopes": "symbol",
            "fields": "ensemblgene,symbol",
            "species": "human",
        }
        try:
            r = requests.post(
                MYGENE_QUERY_URL,
                json=body,
                headers=headers,
                verify=verify,
                timeout=timeout,
            )
            if r.status_code >= 400:
                last_err = f"MyGene symbol lookup HTTP {r.status_code}: {r.text[:400]}"
                continue
            data = r.json()
        except requests.RequestException as exc:
            last_err = f"MyGene symbol lookup {type(exc).__name__}: {exc}"
            continue

        rows: list[dict[str, Any]] = []
        if isinstance(data, list):
            rows = [x for x in data if isinstance(x, dict)]
        else:
            rows = _parse_mygene_hits(data)

        for item in rows:
            if item.get("notfound"):
                continue
            sym = item.get("symbol") or item.get("query")
            if sym:
                sym_u = str(sym).strip().upper()
            else:
                continue
            eg = _ensembl_from_hit(item)
            if eg is None:
                continue
            if eg in osd569_ensembl_ids:
                mapping[eg] = sym_u if len(sym_u) <= 40 else str(sym).strip()

    return mapping, last_err
