#!/usr/bin/env python3
"""
Annotate TIM barrel (beta-alpha)_8 motifs for all proteins in tb_canonical_analysis.

For each protein:
  1. Download the AlphaFold F1 PDB structure (latest version, cached locally).
  2. Run pydssp to assign per-residue secondary structure (H / E / -).
  3. Identify the 8 beta-alpha motifs inside the stored domain region.
  4. Write the result to tb_canonical_analysis.motif_annotations and dssp_source.

PDB files are cached under data/alphafold_pdb/ to avoid re-downloading.

Usage
-----
    python scripts/annotate_motifs.py
    python scripts/annotate_motifs.py --db db/protein_data.db
    python scripts/annotate_motifs.py --limit 50        # first N proteins (testing)
    python scripts/annotate_motifs.py --rerun           # overwrite existing annotations
    python scripts/annotate_motifs.py --log-level DEBUG
"""

import argparse
import json
import logging
import sqlite3
import sys
import time
from pathlib import Path

import numpy as np
import pydssp
import requests

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.analysis.motif_annotator import identify_ba_motifs
from protein_data_collector.config import get_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

_AF_VERSIONS   = [6, 5, 4, 3, 2]   # try newest first
_AF_URL        = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{ver}.pdb"
_REQUEST_DELAY = 0.15               # seconds between downloads (be polite)


# ---------------------------------------------------------------------------
# AlphaFold PDB fetching with local cache
# ---------------------------------------------------------------------------

def _pdb_cache_path(uid: str, cache_dir: Path) -> Path:
    return cache_dir / f"{uid}.pdb"


def fetch_alphafold_pdb(uid: str, cache_dir: Path, session: requests.Session) -> str | None:
    """
    Return PDB text for *uid*, downloading and caching if not already present.
    Tries AlphaFold model versions newest-to-oldest.  Returns None on failure.
    """
    path = _pdb_cache_path(uid, cache_dir)
    if path.exists():
        return path.read_text()

    for ver in _AF_VERSIONS:
        url = _AF_URL.format(uid=uid, ver=ver)
        try:
            r = session.get(url, timeout=30)
            time.sleep(_REQUEST_DELAY)
        except requests.exceptions.RequestException as e:
            logger.debug("Network error for %s v%d: %s", uid, ver, e)
            continue

        if r.status_code == 200:
            path.write_text(r.text)
            logger.debug("Downloaded %s (v%d, %d bytes)", uid, ver, len(r.text))
            return r.text
        elif r.status_code == 404:
            continue
        else:
            logger.warning("Unexpected status %d for %s v%d", r.status_code, uid, ver)

    logger.warning("No AlphaFold structure found for %s", uid)
    return None


# ---------------------------------------------------------------------------
# Secondary structure + motif identification
# ---------------------------------------------------------------------------

def run_dssp(pdb_text: str) -> np.ndarray | None:
    """
    Run pydssp on *pdb_text* and return the c3 SS array (H / E / -).
    Returns None if pydssp fails (e.g. missing backbone atoms).
    """
    try:
        coord = pydssp.read_pdbtext(pdb_text)
        return pydssp.assign(coord, out_type='c3')
    except Exception as e:
        logger.debug("pydssp failed: %s", e)
        return None


def annotate_protein(
    uid: str,
    domain_start: int,
    domain_end: int,
    pdb_text: str,
) -> tuple[list[dict] | None, str]:
    """
    Run DSSP and motif identification for one protein.

    Returns (motifs, dssp_source) where dssp_source is 'alphafold_v{N}' or
    the empty string on failure.
    """
    ss = run_dssp(pdb_text)
    if ss is None:
        return None, ""

    protein_len = len(ss)
    if domain_end > protein_len:
        logger.warning(
            "%s domain_end=%d > PDB length=%d — clamping", uid, domain_end, protein_len
        )
        domain_end = protein_len

    motifs = identify_ba_motifs(ss, domain_start, domain_end)

    # Infer version from cached filename — stored as just the pdb; version tag
    # is embedded in the URL we used; for the DB record we store a generic tag.
    dssp_source = "alphafold"
    return motifs, dssp_source


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def run(
    conn: sqlite3.Connection,
    cache_dir: Path,
    limit: int | None,
    rerun: bool,
) -> dict:
    where = "" if rerun else "AND motif_annotations IS NULL"
    query = f"""
        SELECT ca.uniprot_id, iso.alphafold_id, ca.domain_start, ca.domain_end
        FROM tb_canonical_analysis ca
        JOIN tb_isoforms iso ON iso.uniprot_id = ca.uniprot_id
          AND iso.is_canonical = 1
        WHERE ca.domain_start IS NOT NULL
          AND ca.domain_end   IS NOT NULL
          {where}
        ORDER BY ca.uniprot_id
    """
    if limit:
        query += f" LIMIT {limit}"

    rows = conn.execute(query).fetchall()
    total = len(rows)
    logger.info("Proteins to annotate: %d", total)

    session = requests.Session()
    session.headers["User-Agent"] = "TIMBarrelResearch/1.0"

    stats = {"annotated": 0, "no_structure": 0, "no_motifs": 0, "partial": 0}

    for i, (uniprot_id, alphafold_id, domain_start, domain_end) in enumerate(rows, 1):
        af_id = alphafold_id or uniprot_id

        pdb_text = fetch_alphafold_pdb(af_id, cache_dir, session)
        if pdb_text is None:
            stats["no_structure"] += 1
            conn.execute(
                "UPDATE tb_canonical_analysis SET dssp_source='not_found' WHERE uniprot_id=?",
                (uniprot_id,),
            )
            if i % 100 == 0:
                conn.commit()
            continue

        motifs, dssp_source = annotate_protein(uniprot_id, domain_start, domain_end, pdb_text)

        if motifs is None:
            stats["no_motifs"] += 1
            conn.execute(
                "UPDATE tb_canonical_analysis SET dssp_source='dssp_failed' WHERE uniprot_id=?",
                (uniprot_id,),
            )
        else:
            if len(motifs) < 8:
                stats["partial"] += 1
                logger.debug("%s: only %d motifs found", uniprot_id, len(motifs))
            else:
                stats["annotated"] += 1

            conn.execute(
                """UPDATE tb_canonical_analysis
                   SET motif_annotations=?, dssp_source=?
                   WHERE uniprot_id=?""",
                (json.dumps(motifs), dssp_source, uniprot_id),
            )

        if i % 50 == 0:
            conn.commit()
            logger.info(
                "  %d / %d | annotated=%d partial=%d no_structure=%d",
                i, total, stats["annotated"], stats["partial"], stats["no_structure"],
            )

    conn.commit()
    return stats


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(conn: sqlite3.Connection, stats: dict) -> None:
    total = conn.execute(
        "SELECT COUNT(*) FROM tb_canonical_analysis WHERE domain_start IS NOT NULL"
    ).fetchone()[0]
    with_8 = conn.execute(
        "SELECT COUNT(*) FROM tb_canonical_analysis WHERE motif_annotations IS NOT NULL "
        "AND json_array_length(motif_annotations) = 8"
    ).fetchone()[0]
    with_any = conn.execute(
        "SELECT COUNT(*) FROM tb_canonical_analysis WHERE motif_annotations IS NOT NULL "
        "AND json_array_length(motif_annotations) > 0"
    ).fetchone()[0]

    print(f"\n{'='*60}")
    print("  Motif annotation summary")
    print(f"{'='*60}")
    print(f"  Proteins with domain location  : {total}")
    print(f"  Full 8 motifs annotated        : {with_8}")
    print(f"  Partial (1-7 motifs)           : {with_any - with_8}")
    print(f"  No AlphaFold structure         : {stats['no_structure']}")
    print(f"  DSSP failed                    : {stats['no_motifs']}")
    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate TIM barrel beta-alpha motifs via AlphaFold + pydssp"
    )
    parser.add_argument("--db",        default=None)
    parser.add_argument("--cache-dir", default="data/alphafold_pdb",
                        help="Directory for cached AlphaFold PDB files")
    parser.add_argument("--limit",     type=int, default=None,
                        help="Process at most N proteins (useful for testing)")
    parser.add_argument("--rerun",     action="store_true",
                        help="Re-annotate proteins that already have motif_annotations")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    stats = run(conn, cache_dir, args.limit, args.rerun)
    print_summary(conn, stats)

    conn.close()


if __name__ == "__main__":
    main()
