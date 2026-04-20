#!/usr/bin/env python3
"""
HMMER3-based domain boundary detection for canonical isoforms.

Uses pyhmmer (Python bindings for HMMER3) to scan canonical protein sequences
against the HMM profiles for all entries in the domain's entries table.
Updates tim_barrel_location and tim_barrel_sequence in the isoforms table.

HMMs are fetched from the InterPro API and cached in data/hmm/.

Usage:
    python scripts/run_hmmer.py
    python scripts/run_hmmer.py --domain tim_barrel
    python scripts/run_hmmer.py --domain tim_barrel --evalue 1e-5
    python scripts/run_hmmer.py --rebuild-hmms   # re-fetch HMMs even if cached
"""

import argparse
import gzip
import io
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests
import sqlite3

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.config import DOMAINS, ORGANISMS, get_config
from protein_data_collector.database.connection import get_connection

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# HMM fetching
# ---------------------------------------------------------------------------

_INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api"


def _fetch_hmm_bytes(db: str, accession: str) -> Optional[bytes]:
    """Fetch the gzipped HMM for a single entry from InterPro API."""
    url = f"{_INTERPRO_BASE}/entry/{db}/{accession}/?annotation=hmm"
    try:
        r = requests.get(url, timeout=60)
        if r.status_code == 200 and r.content:
            return r.content
        logger.debug("No HMM for %s (%s): HTTP %d", accession, db, r.status_code)
        return None
    except requests.RequestException as e:
        logger.warning("Failed to fetch HMM for %s: %s", accession, e)
        return None


def _get_pfam_members_for_interpro(accession: str) -> List[str]:
    """Return Pfam member accessions of an InterPro entry."""
    url = f"{_INTERPRO_BASE}/entry/interpro/{accession}/"
    try:
        r = requests.get(url, timeout=30)
        if r.ok:
            members = r.json().get("metadata", {}).get("member_databases") or {}
            return list(members.get("pfam", {}).keys())
    except requests.RequestException:
        pass
    return []


def fetch_all_hmms(entries_table: str, hmm_dir: Path, db_path: str, rebuild: bool = False) -> Path:
    """
    Fetch HMMs for all entries in entries_table and write to a single combined
    HMM file in hmm_dir.  Returns the path to the combined file.

    Skips fetch if the file exists and rebuild=False.
    """
    hmm_dir.mkdir(parents=True, exist_ok=True)
    combined_path = hmm_dir / f"{entries_table}.hmm"

    if combined_path.exists() and not rebuild:
        logger.info("HMM database already exists at %s — skipping fetch (use --rebuild-hmms to re-fetch)", combined_path)
        return combined_path

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    entries = conn.execute(f"SELECT accession, entry_type FROM {entries_table}").fetchall()
    conn.close()

    # Collect all (db, accession) pairs to fetch
    to_fetch: List[Tuple[str, str]] = []
    pfam_seen: set = set()

    for entry in entries:
        acc, et = entry["accession"], entry["entry_type"]
        if et == "pfam":
            if acc not in pfam_seen:
                to_fetch.append(("pfam", acc))
                pfam_seen.add(acc)
        elif et == "cathgene3d":
            to_fetch.append(("cathgene3d", acc))
        elif et == "interpro":
            # Use Pfam member HMMs (InterPro entries are groupings, not HMM families)
            for pfam_acc in _get_pfam_members_for_interpro(acc):
                if pfam_acc not in pfam_seen:
                    to_fetch.append(("pfam", pfam_acc))
                    pfam_seen.add(pfam_acc)
            time.sleep(0.1)

    logger.info("Fetching %d HMMs for %s ...", len(to_fetch), entries_table)

    hmm_bytes_list: List[bytes] = []
    failed = 0
    for i, (db, acc) in enumerate(to_fetch, 1):
        raw = _fetch_hmm_bytes(db, acc)
        if raw:
            # Decompress gzip and keep the plain HMM text
            try:
                hmm_bytes_list.append(gzip.decompress(raw))
            except (gzip.BadGzipFile, OSError):
                hmm_bytes_list.append(raw)  # already plain text
        else:
            failed += 1
        if i % 10 == 0 or i == len(to_fetch):
            logger.info("  %d/%d fetched (%d failed so far)", i, len(to_fetch), failed)
        time.sleep(0.15)

    combined = b"\n".join(hmm_bytes_list)
    combined_path.write_bytes(combined)
    logger.info("Wrote %d HMMs to %s (%d failed)", len(hmm_bytes_list), combined_path, failed)
    return combined_path


# ---------------------------------------------------------------------------
# Sequence extraction
# ---------------------------------------------------------------------------

def load_canonical_sequences(isoform_table: str, db_path: str) -> List[Tuple[str, str, str]]:
    """Return list of (isoform_id, uniprot_id, sequence) for all canonical, non-fragment isoforms."""
    conn = sqlite3.connect(db_path)
    rows = conn.execute(
        f"SELECT isoform_id, uniprot_id, sequence FROM {isoform_table} "
        f"WHERE is_canonical=1 AND is_fragment=0 AND sequence IS NOT NULL"
    ).fetchall()
    conn.close()
    return [(r[0], r[1], r[2]) for r in rows]


# ---------------------------------------------------------------------------
# HMMER scan
# ---------------------------------------------------------------------------

def run_hmmscan(
    sequences: List[Tuple[str, str, str]],
    hmm_path: Path,
    evalue_threshold: float = 1e-5,
) -> Dict[str, dict]:
    """
    Scan all sequences against the HMM database.

    Returns a dict keyed by uniprot_id with the best-scoring domain hit:
        { uniprot_id: {start, end, length, hmm_name, score, evalue} }
    """
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextSequence

    alphabet = pyhmmer.easel.Alphabet.amino()

    # Load HMMs
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        profiles = list(hf)
    logger.info("Loaded %d HMM profiles from %s", len(profiles), hmm_path)

    # Build digital sequences
    digital_seqs = []
    for iso_id, uid, seq in sequences:
        ts = TextSequence(name=iso_id.encode(), description=uid.encode(), sequence=seq)
        digital_seqs.append(ts.digitize(alphabet))

    logger.info("Scanning %d sequences ...", len(digital_seqs))

    best_hits: Dict[str, dict] = {}

    # hmmscan: query=sequences, target=profiles
    for hits in pyhmmer.hmmer.hmmscan(digital_seqs, profiles, E=evalue_threshold, cpus=0):
        raw_uid    = hits.query.description
        raw_iso_id = hits.query.name
        uid    = raw_uid.decode()    if isinstance(raw_uid,    bytes) else raw_uid
        iso_id = raw_iso_id.decode() if isinstance(raw_iso_id, bytes) else raw_iso_id

        best_domain = None
        best_score = -1.0

        for hit in hits:
            if hit.evalue > evalue_threshold:
                continue
            for domain in hit.domains.included:
                if domain.score > best_score:
                    best_score = domain.score
                    best_domain = {
                        "uniprot_id":  uid,
                        "isoform_id":  iso_id,
                        "hmm_name":    hit.name.decode() if isinstance(hit.name, bytes) else hit.name,
                        "start":       domain.env_from,   # 1-based
                        "end":         domain.env_to,     # 1-based
                        "length":      domain.env_to - domain.env_from + 1,
                        "score":       round(domain.score, 2),
                        "evalue":      hit.evalue,
                        "source":      "hmmer3",
                    }

        if best_domain:
            # Keep only the best hit per protein
            if uid not in best_hits or best_domain["score"] > best_hits[uid]["score"]:
                best_hits[uid] = best_domain

    logger.info("HMMER found domain hits for %d / %d proteins", len(best_hits), len(sequences))
    return best_hits


# ---------------------------------------------------------------------------
# DB update
# ---------------------------------------------------------------------------

def update_db(
    hits: Dict[str, dict],
    isoform_table: str,
    db_path: str,
    overwrite: bool = False,
) -> Tuple[int, int]:
    """
    Write HMMER-derived domain locations to the isoforms table.

    If overwrite=False (default), only updates rows where tim_barrel_location is NULL
    or the existing source is not 'interpro_api' (i.e. don't clobber good API data).
    If overwrite=True, updates all rows regardless of existing source.

    Returns (updated, skipped).
    """
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    updated = skipped = 0

    for uid, hit in hits.items():
        row = conn.execute(
            f"SELECT isoform_id, sequence, tim_barrel_location FROM {isoform_table} "
            f"WHERE uniprot_id=? AND is_canonical=1",
            (uid,)
        ).fetchone()
        if not row:
            continue

        existing_loc = row["tim_barrel_location"]
        if existing_loc and not overwrite:
            try:
                existing_source = json.loads(existing_loc).get("source", "")
            except (json.JSONDecodeError, AttributeError):
                existing_source = ""
            # Keep InterPro API data; only fill in where it's missing or non-API
            if existing_source == "interpro_api":
                skipped += 1
                continue

        start = hit["start"]
        end   = hit["end"]
        seq   = row["sequence"]
        domain_seq = seq[start - 1:end] if seq and len(seq) >= end else None

        loc_json = json.dumps({
            "domain_id": hit["hmm_name"],
            "start":     start,
            "end":       end,
            "length":    end - start + 1,
            "score":     hit["score"],
            "evalue":    hit["evalue"],
            "source":    "hmmer3",
        })

        conn.execute(
            f"UPDATE {isoform_table} "
            f"SET tim_barrel_location=?, tim_barrel_sequence=? "
            f"WHERE uniprot_id=? AND is_canonical=1",
            (loc_json, domain_seq, uid)
        )
        updated += 1

    conn.commit()
    conn.close()
    return updated, skipped


# ---------------------------------------------------------------------------
# Report helpers
# ---------------------------------------------------------------------------

def _compare_locations(hits: Dict[str, dict], isoform_table: str, db_path: str) -> None:
    """Log a comparison between HMMER boundaries and existing InterPro API boundaries."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    diffs = []
    for uid, hit in hits.items():
        row = conn.execute(
            f"SELECT tim_barrel_location FROM {isoform_table} WHERE uniprot_id=? AND is_canonical=1",
            (uid,)
        ).fetchone()
        if not row or not row["tim_barrel_location"]:
            continue
        try:
            existing = json.loads(row["tim_barrel_location"])
        except json.JSONDecodeError:
            continue
        if existing.get("source") != "interpro_api":
            continue
        d_start = abs(hit["start"] - existing["start"])
        d_end   = abs(hit["end"]   - existing["end"])
        if d_start > 10 or d_end > 10:
            diffs.append((uid, existing["start"], existing["end"], hit["start"], hit["end"]))

    conn.close()

    if diffs:
        logger.info("%d proteins have >10 aa boundary difference between InterPro API and HMMER:", len(diffs))
        for uid, es, ee, hs, he in diffs[:10]:
            logger.info("  %s: API=[%d,%d] HMMER=[%d,%d]", uid, es, ee, hs, he)
    else:
        logger.info("All HMMER boundaries agree with InterPro API within 10 aa.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Run HMMER3 domain scan on canonical isoforms")
    parser.add_argument("--domain",       default="tim_barrel", choices=list(DOMAINS))
    parser.add_argument("--organism",     default="homo_sapiens", choices=list(ORGANISMS))
    parser.add_argument("--db",           default=None)
    parser.add_argument("--evalue",       type=float, default=1e-5)
    parser.add_argument("--overwrite",    action="store_true",
                        help="Overwrite existing InterPro API boundaries with HMMER results")
    parser.add_argument("--rebuild-hmms", action="store_true",
                        help="Re-fetch HMMs even if cached file exists")
    parser.add_argument("--hmm-dir",      default="data/hmm")
    parser.add_argument("--log-level",    default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)-8s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    domain_cfg   = DOMAINS[args.domain]
    organism_cfg = ORGANISMS[args.organism]
    db_path      = args.db or get_config().db_path
    hmm_dir      = Path(args.hmm_dir)

    isoform_table = organism_cfg.isoform_table(domain_cfg)

    # Step 1 — fetch / load HMMs
    hmm_path = fetch_all_hmms(
        entries_table=domain_cfg.entries_table,
        hmm_dir=hmm_dir,
        db_path=db_path,
        rebuild=args.rebuild_hmms,
    )

    # Step 2 — load sequences
    sequences = load_canonical_sequences(isoform_table, db_path)
    logger.info("Loaded %d canonical non-fragment sequences from %s", len(sequences), isoform_table)

    # Step 3 — scan
    hits = run_hmmscan(sequences, hmm_path, evalue_threshold=args.evalue)

    # Step 4 — compare with existing boundaries (info only)
    _compare_locations(hits, isoform_table, db_path)

    # Step 5 — update DB
    updated, skipped = update_db(hits, isoform_table, db_path, overwrite=args.overwrite)

    print(f"\n{'='*55}")
    print(f"  Domain   : {domain_cfg.display_name}")
    print(f"  Organism : {organism_cfg.display_name}")
    print(f"  E-value  : {args.evalue}")
    print(f"{'='*55}")
    print(f"  Sequences scanned     : {len(sequences)}")
    print(f"  Proteins with hit     : {len(hits)}")
    print(f"  Proteins with no hit  : {len(sequences) - len(hits)}")
    print(f"  DB rows updated       : {updated}")
    print(f"  DB rows skipped       : {skipped}  (existing interpro_api boundary kept)")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
