#!/usr/bin/env python3
"""
Cross-validate DSSP-based TIM barrel motif annotations using per-family profile HMMs.

TIM barrels are so sequence-divergent that no single HMM covers all families.
This script builds one HMM per family (grouped by tim_barrel_accession), then
uses each family's HMM to validate proteins in that family.

Strategy per family
-------------------
  - Families with >= 5 complete (8-motif) proteins:
      jackhmmer seeded from the best-scoring complete protein; max 5 iterations.
  - Families with 1-4 complete proteins:
      phmmer (single-sequence HMM) seeded from the single best-scoring protein.
  - Families with 0 complete proteins:
      No HMM built; proteins are marked 'no_profile'.

For each protein the result stored in hmmer_annotations (JSON) contains:
    hit        : bool — whether the family HMM found a significant hit
    score      : float — bit score (even for non-significant hits if any)
    evalue     : float — E-value
    env_from   : int — envelope start in full-protein coords (1-based)
    env_to     : int — envelope end in full-protein coords
    family     : str — the tim_barrel_accession for this protein
    profile    : str — seed protein UID used to build the HMM

hmmer_source is set to 'family_jackhmmer' or 'family_phmmer' or 'no_profile'.

Usage
-----
    python scripts/cross_validate_hmmer.py
    python scripts/cross_validate_hmmer.py --db db/protein_data.db
    python scripts/cross_validate_hmmer.py --iterations 5
    python scripts/cross_validate_hmmer.py --evalue 1e-3
"""

import argparse
import json
import logging
import sqlite3
import sys
from pathlib import Path

import pyhmmer

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.config import get_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

_ALPHABET        = pyhmmer.easel.Alphabet.amino()
_JACKHMMER_MIN   = 5     # min complete proteins in family to use jackhmmer vs phmmer
_DEFAULT_EVALUE  = 1e-3


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(tb_canonical_analysis)")}
    if "hmmer_annotations" not in cols:
        conn.execute("ALTER TABLE tb_canonical_analysis ADD COLUMN hmmer_annotations TEXT")
        logger.info("Added hmmer_annotations to tb_canonical_analysis")
    conn.commit()


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def _digitize(uid: str, seq: str) -> pyhmmer.easel.DigitalSequence:
    return pyhmmer.easel.TextSequence(
        name=uid.encode(), sequence=seq
    ).digitize(_ALPHABET)


# ---------------------------------------------------------------------------
# Per-family data loading
# ---------------------------------------------------------------------------

def load_families(conn: sqlite3.Connection) -> dict[str, dict]:
    """
    Return {accession: {
        'name': str,
        'complete': [(uid, ds, de, seq), ...],   # 8-motif proteins
        'all':      [(uid, ds, de, seq), ...],   # all proteins in family
    }}
    """
    rows = conn.execute("""
        SELECT ca.uniprot_id, ca.domain_start, ca.domain_end, ca.domain_sequence,
               json_array_length(ca.motif_annotations) as n_motifs,
               p.tim_barrel_accession, e.name as entry_name
        FROM tb_canonical_analysis ca
        JOIN tb_proteins p ON p.uniprot_id = ca.uniprot_id
        JOIN tb_entries e ON e.accession = p.tim_barrel_accession
        WHERE ca.domain_sequence IS NOT NULL
        ORDER BY p.tim_barrel_accession, ca.uniprot_id
    """).fetchall()

    families: dict[str, dict] = {}
    for uid, ds, de, seq, n_motifs, accession, entry_name in rows:
        fam = families.setdefault(accession, {
            "name":     entry_name,
            "complete": [],
            "all":      [],
        })
        try:
            dseq = _digitize(uid, seq)
        except Exception:
            continue
        fam["all"].append((uid, ds, de, dseq))
        if n_motifs == 8:
            fam["complete"].append((uid, ds, de, dseq))

    return families


# ---------------------------------------------------------------------------
# Profile building
# ---------------------------------------------------------------------------

def _pick_seed(complete: list) -> tuple:
    """Return the first entry from the complete list (ordered by uid)."""
    return complete[0]


def build_jackhmmer(complete: list, max_iterations: int) -> pyhmmer.plan7.HMM | None:
    seed_uid, _, _, seed_seq = _pick_seed(complete)
    db_seqs = [s for _, _, _, s in complete]
    final_hmm = None
    try:
        for result in pyhmmer.hmmer.jackhmmer(
            seed_seq, db_seqs, max_iterations=max_iterations
        ):
            if result.hmm is not None:
                final_hmm = result.hmm
            if result.converged:
                break
    except Exception as e:
        logger.warning("jackhmmer failed (seed=%s): %s", seed_uid, e)
        return None
    return final_hmm


def build_phmmer_profile(complete: list) -> pyhmmer.plan7.HMM | None:
    """Build a single-sequence HMM (phmmer-style) from the best complete protein."""
    seed_uid, _, _, seed_seq = _pick_seed(complete)
    builder = pyhmmer.plan7.Builder(_ALPHABET)
    bg = pyhmmer.plan7.Background(_ALPHABET)
    try:
        hmm, _, _ = builder.build(seed_seq, bg)
        return hmm
    except Exception as e:
        logger.warning("phmmer build failed (seed=%s): %s", seed_uid, e)
        return None


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def search_family(
    hmm: pyhmmer.plan7.HMM,
    targets: list,
    evalue_thr: float,
) -> dict[str, dict]:
    """
    Run hmmsearch with *hmm* against *targets* = [(uid, domain_start, domain_end, dseq), ...].
    Returns {uid: annotation_dict}.
    """
    coord_map  = {uid: (ds, de) for uid, ds, de, _ in targets}
    target_seqs = [s for _, _, _, s in targets]
    results: dict[str, dict] = {}

    for hits in pyhmmer.hmmer.hmmsearch([hmm], target_seqs):
        for hit in hits:
            uid = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            ds, de = coord_map.get(uid, (None, None))

            if not hit.included or hit.evalue > evalue_thr:
                results[uid] = {
                    "hit":   False,
                    "score": round(hit.score, 2),
                    "evalue": hit.evalue,
                }
                continue

            best = hit.best_domain
            env_from = (ds + best.env_from - 1) if ds else best.env_from
            env_to   = (ds + best.env_to   - 1) if ds else best.env_to

            results[uid] = {
                "hit":      True,
                "score":    round(hit.score, 2),
                "evalue":   hit.evalue,
                "env_from": env_from,
                "env_to":   env_to,
                "n_domains": len(hit.domains),
            }

    return results


# ---------------------------------------------------------------------------
# Storage
# ---------------------------------------------------------------------------

def store_family_results(
    conn: sqlite3.Connection,
    accession: str,
    search_res: dict[str, dict],
    source: str,
    seed_uid: str,
) -> None:
    for uid, ann in search_res.items():
        ann["family"]  = accession
        ann["profile"] = seed_uid
        conn.execute("""
            UPDATE tb_canonical_analysis
            SET hmmer_annotations=?, hmmer_source=?
            WHERE uniprot_id=?
        """, (json.dumps(ann), source, uid))

    # Proteins in this family with no hit at all (not in hmmsearch output)
    all_uids = {uid for uid in search_res}
    family_uids = conn.execute("""
        SELECT ca.uniprot_id FROM tb_canonical_analysis ca
        JOIN tb_proteins p ON p.uniprot_id = ca.uniprot_id
        WHERE p.tim_barrel_accession=? AND ca.domain_sequence IS NOT NULL
          AND ca.hmmer_annotations IS NULL
    """, (accession,)).fetchall()
    for (uid,) in family_uids:
        if uid not in all_uids:
            conn.execute("""
                UPDATE tb_canonical_analysis
                SET hmmer_annotations=?, hmmer_source=?
                WHERE uniprot_id=?
            """, (json.dumps({"hit": False, "family": accession, "profile": seed_uid}), source, uid))

    conn.commit()


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def print_comparison(conn: sqlite3.Connection) -> None:
    rows = conn.execute("""
        SELECT
            json_array_length(motif_annotations)          AS n_motifs,
            json_extract(hmmer_annotations, '$.hit')      AS hmmer_hit,
            hmmer_source,
            json_extract(hmmer_annotations, '$.score')    AS score
        FROM tb_canonical_analysis
        WHERE hmmer_annotations IS NOT NULL
          AND hmmer_source != 'no_profile'
    """).fetchall()

    buckets: dict = {}
    for n_motifs, hmmer_hit, source, score in rows:
        key = n_motifs if n_motifs is not None else -1
        b = buckets.setdefault(key, {"total": 0, "hit": 0, "scores": []})
        b["total"] += 1
        if hmmer_hit:
            b["hit"] += 1
        if score is not None:
            b["scores"].append(score)

    print(f"\n{'='*72}")
    print("  DSSP vs. family HMM cross-validation")
    print(f"{'='*72}")
    print(f"  {'DSSP motifs':>12}  {'proteins':>9}  {'HMM hit':>9}  {'HMM hit%':>9}  {'med.score':>9}")
    print(f"  {'-'*12}  {'-'*9}  {'-'*9}  {'-'*9}  {'-'*9}")
    for key in sorted(buckets):
        b = buckets[key]
        label = str(key) if key >= 0 else "none"
        pct   = 100 * b["hit"] / b["total"] if b["total"] else 0
        sc = b["scores"]
        med = sorted(sc)[len(sc) // 2] if sc else float("nan")
        print(f"  {label:>12}  {b['total']:>9}  {b['hit']:>9}  {pct:>8.0f}%  {med:>9.1f}")
    print(f"{'='*72}")

    total   = conn.execute("SELECT COUNT(*) FROM tb_canonical_analysis WHERE hmmer_source != 'no_profile'").fetchone()[0]
    hit     = conn.execute("SELECT COUNT(*) FROM tb_canonical_analysis WHERE json_extract(hmmer_annotations,'$.hit')=1").fetchone()[0]
    no_prof = conn.execute("SELECT COUNT(*) FROM tb_canonical_analysis WHERE hmmer_source='no_profile'").fetchone()[0]
    full8   = conn.execute("SELECT COUNT(*) FROM tb_canonical_analysis WHERE json_array_length(motif_annotations)=8 AND json_extract(hmmer_annotations,'$.hit')=1").fetchone()[0]
    full8_t = conn.execute("SELECT COUNT(*) FROM tb_canonical_analysis WHERE json_array_length(motif_annotations)=8 AND hmmer_source != 'no_profile'").fetchone()[0]
    print(f"\n  Proteins with a family HMM       : {total}")
    print(f"  Proteins without a profile        : {no_prof}")
    print(f"  Significant HMM hits              : {hit} / {total}")
    print(f"  8-motif confirmed by family HMM   : {full8} / {full8_t}")
    print(f"{'='*72}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Per-family HMM cross-validation of TIM barrel motif annotations"
    )
    parser.add_argument("--db",         default=None)
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--evalue",     type=float, default=_DEFAULT_EVALUE)
    parser.add_argument("--log-level",  default="INFO")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    ensure_columns(conn)

    # Reset previous results
    conn.execute("UPDATE tb_canonical_analysis SET hmmer_annotations=NULL, hmmer_source=NULL")
    conn.commit()

    families = load_families(conn)
    logger.info("Families to process: %d", len(families))

    for accession, fam in sorted(families.items()):
        n_all      = len(fam["all"])
        n_complete = len(fam["complete"])
        name_short = fam["name"][:40] if fam["name"] else accession

        if n_complete == 0:
            logger.info("  [no_profile] %s — %s (%d proteins, 0 complete)",
                        accession, name_short, n_all)
            for uid, ds, de, _ in fam["all"]:
                conn.execute("""
                    UPDATE tb_canonical_analysis
                    SET hmmer_annotations=?, hmmer_source='no_profile'
                    WHERE uniprot_id=?
                """, (json.dumps({"hit": False, "family": accession}), uid))
            conn.commit()
            continue

        seed_uid = fam["complete"][0][0]

        if n_complete >= _JACKHMMER_MIN:
            logger.info("  [jackhmmer ] %s — %s (%d total, %d complete, seed=%s)",
                        accession, name_short, n_all, n_complete, seed_uid)
            hmm = build_jackhmmer(fam["complete"], args.iterations)
            source = "family_jackhmmer"
        else:
            logger.info("  [phmmer    ] %s — %s (%d total, %d complete, seed=%s)",
                        accession, name_short, n_all, n_complete, seed_uid)
            hmm = build_phmmer_profile(fam["complete"])
            source = "family_phmmer"

        if hmm is None:
            logger.warning("    HMM build failed — marking as no_profile")
            for uid, _, _, _ in fam["all"]:
                conn.execute("""
                    UPDATE tb_canonical_analysis
                    SET hmmer_annotations=?, hmmer_source='no_profile'
                    WHERE uniprot_id=?
                """, (json.dumps({"hit": False, "family": accession}), uid))
            conn.commit()
            continue

        search_res = search_family(hmm, fam["all"], args.evalue)
        hit_count  = sum(1 for v in search_res.values() if v.get("hit"))
        logger.info("    hits: %d / %d", hit_count, n_all)

        store_family_results(conn, accession, search_res, source, seed_uid)

    print_comparison(conn)
    conn.close()


if __name__ == "__main__":
    main()
