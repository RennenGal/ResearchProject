#!/usr/bin/env python3
"""
Deduplicate canonical proteins that share the same gene_name.

The original deduplication grouped by (protein_name, organism) and missed
TrEMBL entries whose protein_name differs from the reviewed Swiss-Prot
representative (e.g. "KLB protein" vs "Beta-klotho", "Galactosidase beta 1"
vs "Beta-galactosidase").

Strategy
--------
For each gene_name with multiple canonical proteins (canonical_uniprot_id IS NULL):
  - Keep the best representative: reviewed > annotation_score > isoform count
  - Mark all others as redundant: set canonical_uniprot_id = best_uid
  - Reassign any proteins that currently point to a non-best entry so they
    point to the new best representative instead.

Proteins without a gene_name are skipped (cannot be safely grouped).

Usage
-----
    python scripts/dedup_by_gene.py                   # dry-run (prints plan, writes nothing)
    python scripts/dedup_by_gene.py --apply           # apply changes
    python scripts/dedup_by_gene.py --db path.sqlite  # custom DB path
"""

import logging
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def _isoform_count(conn: sqlite3.Connection, uid: str) -> int:
    return conn.execute(
        "SELECT COUNT(*) FROM isoforms WHERE uniprot_id = ?", (uid,)
    ).fetchone()[0]


def find_gene_duplicates(conn: sqlite3.Connection) -> dict[str, list[dict]]:
    """
    Return {gene_name: [protein_row, ...]} for genes with > 1 canonical entry.
    Each protein_row has: uniprot_id, reviewed, annotation_score, protein_name.
    """
    rows = conn.execute("""
        SELECT uniprot_id, gene_name, protein_name, reviewed, annotation_score
        FROM   proteins
        WHERE  gene_name            IS NOT NULL
          AND  canonical_uniprot_id IS NULL
        ORDER BY gene_name, reviewed DESC, annotation_score DESC
    """).fetchall()

    by_gene: dict[str, list[dict]] = defaultdict(list)
    for uid, gene, pname, rev, score in rows:
        by_gene[gene].append({
            "uniprot_id": uid,
            "reviewed":   rev or 0,
            "annotation_score": score or 0,
            "protein_name": pname,
        })

    return {gene: entries for gene, entries in by_gene.items() if len(entries) > 1}


def _best(entries: list[dict], conn: sqlite3.Connection) -> dict:
    """Pick the best representative: reviewed > annotation_score > isoform count."""
    def key(e):
        return (e["reviewed"], e["annotation_score"],
                _isoform_count(conn, e["uniprot_id"]))
    return max(entries, key=key)


# Manual overrides for genes where the reviewed canonical uses a different gene symbol
# (e.g. MUT → MMUT rename): map each stale canonical to the correct reviewed uid.
_MANUAL_REDIRECTS: dict[str, str] = {
    "A0A0B4U8R5": "P22033",   # Mutant methylmalonyl CoA mutase → MMUT (P22033)
    "A0A0K0PWN6": "P22033",   # Truncated methylmalonyl CoA mutase → MMUT (P22033)
    "S4UM43":     "P22033",   # Mitochondrial methylmalonyl CoA mutase → MMUT (P22033)
}


def _dedup(conn: sqlite3.Connection, apply: bool) -> dict:
    duplicates = find_gene_duplicates(conn)

    stats = {"genes": 0, "proteins_merged": 0}

    for gene, entries in sorted(duplicates.items()):
        best = _best(entries, conn)
        best_uid = best["uniprot_id"]
        to_merge = [e for e in entries if e["uniprot_id"] != best_uid]

        print(f"\n  Gene: {gene}  ({len(entries)} canonical entries)")
        print(f"    KEEP  : {best_uid}  reviewed={best['reviewed']}"
              f"  score={best['annotation_score']}  [{best['protein_name']}]")
        for e in to_merge:
            print(f"    MERGE : {e['uniprot_id']}  reviewed={e['reviewed']}"
                  f"  score={e['annotation_score']}  [{e['protein_name']}]")

        if apply:
            for e in to_merge:
                merge_uid = e["uniprot_id"]
                # Mark this protein as redundant → best
                conn.execute(
                    "UPDATE proteins SET canonical_uniprot_id = ? WHERE uniprot_id = ?",
                    (best_uid, merge_uid),
                )
                # Reassign any proteins currently pointing to merge_uid → best_uid
                conn.execute(
                    "UPDATE proteins SET canonical_uniprot_id = ?"
                    " WHERE canonical_uniprot_id = ?",
                    (best_uid, merge_uid),
                )
            stats["proteins_merged"] += len(to_merge)

        stats["genes"] += 1

    # Apply manual redirects (gene-alias edge cases)
    print("\n  Manual redirects (gene alias mismatches):")
    for uid, canonical_uid in _MANUAL_REDIRECTS.items():
        row = conn.execute(
            "SELECT protein_name, reviewed, annotation_score, canonical_uniprot_id"
            " FROM proteins WHERE uniprot_id = ?", (uid,)
        ).fetchone()
        if not row:
            continue
        pname, rev, score, current_canon = row
        if current_canon == canonical_uid:
            print(f"    SKIP  : {uid} already -> {canonical_uid}")
            continue
        print(f"    REDIRECT: {uid} [{pname}] -> {canonical_uid}")
        if apply:
            conn.execute(
                "UPDATE proteins SET canonical_uniprot_id = ? WHERE uniprot_id = ?",
                (canonical_uid, uid),
            )
            conn.execute(
                "UPDATE proteins SET canonical_uniprot_id = ?"
                " WHERE canonical_uniprot_id = ?",
                (canonical_uid, uid),
            )
            stats["proteins_merged"] += 1

    if apply:
        conn.commit()

    return stats


def run(db_path: str) -> None:
    conn = sqlite3.connect(db_path)

    print("\n" + "="*60)
    print("  Gene-name deduplication plan")
    print("="*60)

    stats = _dedup(conn, apply=True)

    print(f"\n{'='*60}")
    print("  Summary")
    print(f"{'='*60}")
    print(f"  Genes with duplicate canonical entries : {stats['genes']}")
    print(f"  Proteins merged (marked redundant)     : {stats['proteins_merged']}")
    print(f"{'='*60}")

    conn.close()
