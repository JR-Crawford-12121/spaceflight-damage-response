"""
Curated starter gene sets for DNA damage-response–associated pathways.

These are exploratory signatures (gene expression / regulation), not direct
measurements of DNA lesions.
"""

from __future__ import annotations

PATHWAY_GENE_SETS: dict[str, list[str]] = {
    "DNA_repair": [
        "TP53",
        "ATM",
        "ATR",
        "BRCA1",
        "BRCA2",
        "RAD51",
        "XPC",
        "XPA",
        "ERCC1",
        "DDB2",
        "GADD45A",
        "GADD45B",
        "MSH2",
        "MSH6",
        "MLH1",
        "PARP1",
        "XRCC1",
        "XRCC5",
        "XRCC6",
        "CHEK1",
        "CHEK2",
    ],
    "p53_pathway": [
        "TP53",
        "CDKN1A",
        "MDM2",
        "BAX",
        "BBC3",
        "PMAIP1",
        "GADD45A",
        "DDB2",
        "SESN1",
        "SESN2",
        "RRM2B",
        "FAS",
        "ZMAT3",
        "BTG2",
    ],
    "oxidative_stress": [
        "NFE2L2",
        "KEAP1",
        "HMOX1",
        "NQO1",
        "SOD1",
        "SOD2",
        "CAT",
        "GPX1",
        "GPX4",
        "PRDX1",
        "PRDX2",
        "TXN",
        "TXNRD1",
        "GSR",
        "SRXN1",
    ],
    "apoptosis": [
        "BAX",
        "BCL2",
        "CASP3",
        "CASP7",
        "CASP8",
        "CASP9",
        "APAF1",
        "FAS",
        "FASLG",
        "TNFRSF10B",
        "BBC3",
        "PMAIP1",
        "BID",
        "BAD",
    ],
    "cell_cycle_checkpoint": [
        "CDKN1A",
        "CHEK1",
        "CHEK2",
        "ATM",
        "ATR",
        "CCNG1",
        "GADD45A",
        "GADD45B",
        "WEE1",
        "CDC25A",
        "CDC25C",
        "RB1",
        "E2F1",
    ],
}


def union_pathway_gene_symbols() -> set[str]:
    """All HGNC symbols appearing in any curated pathway list (uppercase)."""
    out: set[str] = set()
    for genes in PATHWAY_GENE_SETS.values():
        out.update(str(g).strip().upper() for g in genes)
    return out
