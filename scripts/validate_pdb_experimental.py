#!/usr/bin/env python3
"""
Validate TIM barrel motif annotations using experimental PDB structures.

For each protein in tb_canonical_analysis, queries the PDBe best_structures
API to find the highest-quality experimental structure (X-ray or EM, resolution
<= 3.0 Å) that covers the TIM barrel domain.  The selected structure is
downloaded from RCSB, filtered to the relevant chain, run through pydssp, and
motifs are identified using the same identify_ba_motifs() function used for
AlphaFold validation.

UniProt position mapping
------------------------
PDBe reports unp_start/unp_end (UniProt positions) and start/end (PDB auth
sequence numbers).  The offset

    offset = unp_start - start

converts a PDB auth_seq_num to a UniProt position:

    uniprot_pos = auth_seq_num + offset

Results stored in tb_canonical_analysis
----------------------------------------
    pdb_motif_annotations : JSON — same format as motif_annotations
    pdb_source            : "<PDB_ID>_<chain>_<resolution>Å"
                            or 'no_structure' / 'dssp_failed'

Usage
-----
    python scripts/validate_pdb_experimental.py
    python scripts/validate_pdb_experimental.py --db db/protein_data.db
    python scripts/validate_pdb_experimental.py --resolution 2.5
    python scripts/validate_pdb_experimental.py --rerun
    python scripts/validate_pdb_experimental.py --limit 50
"""

import argparse
import json
import logging
import sqlite3
import sys
import time
from collections import OrderedDict
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

_PDBE_API     = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/{uid}"
_RCSB_URL     = "https://files.rcsb.org/download/{pdb_id}.pdb"
_REQUEST_DELAY = 0.15
_DEFAULT_RES   = 3.0
_BACKBONE      = {"N", "CA", "C", "O"}
_ACCEPTED_METHODS = {"X-ray diffraction", "Electron Microscopy"}


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(tb_canonical_analysis)")}
    for col, typedef in [
        ("pdb_motif_annotations", "TEXT"),
        ("pdb_source",            "TEXT"),
    ]:
        if col not in cols:
            conn.execute(f"ALTER TABLE tb_canonical_analysis ADD COLUMN {col} {typedef}")
            logger.info("Added %s to tb_canonical_analysis", col)
    conn.commit()


# ---------------------------------------------------------------------------
# PDBe — pick best structure
# ---------------------------------------------------------------------------

def _covers_domain(entry: dict, domain_start: int, domain_end: int) -> bool:
    unp_s = entry.get("unp_start", 0)
    unp_e = entry.get("unp_end", 0)
    return unp_s <= domain_start and unp_e >= domain_end


def best_structure(
    uid: str,
    domain_start: int,
    domain_end: int,
    max_resolution: float,
    session: requests.Session,
) -> dict | None:
    """
    Return the best PDBe structure entry for *uid* that:
      - is X-ray or EM
      - has resolution <= max_resolution
      - covers [domain_start, domain_end] in UniProt coordinates

    Returns the entry dict (keys: pdb_id, chain_id, resolution, unp_start,
    unp_end, start, end) or None.
    """
    url = _PDBE_API.format(uid=uid)
    try:
        resp = session.get(url, timeout=30)
        time.sleep(_REQUEST_DELAY)
    except requests.exceptions.RequestException as e:
        logger.debug("PDBe request failed for %s: %s", uid, e)
        return None

    if resp.status_code == 404:
        return None
    if not resp.ok:
        logger.warning("PDBe returned %d for %s", resp.status_code, uid)
        return None

    data = resp.json()
    entries = data.get(uid, [])

    candidates = []
    for e in entries:
        method = e.get("experimental_method", "")
        res    = e.get("resolution") or 99.0
        if method not in _ACCEPTED_METHODS:
            continue
        if res > max_resolution:
            continue
        if not _covers_domain(e, domain_start, domain_end):
            continue
        candidates.append(e)

    if not candidates:
        return None

    # Primary sort: resolution; secondary: coverage (larger = better)
    candidates.sort(key=lambda e: (e.get("resolution") or 99.0,
                                   -(e.get("unp_end", 0) - e.get("unp_start", 0))))
    return candidates[0]


# ---------------------------------------------------------------------------
# PDB download + cache
# ---------------------------------------------------------------------------

def fetch_pdb(pdb_id: str, cache_dir: Path, session: requests.Session) -> str | None:
    path = cache_dir / f"{pdb_id}.pdb"
    if path.exists():
        return path.read_text(encoding="utf-8", errors="replace")

    url = _RCSB_URL.format(pdb_id=pdb_id.lower())
    try:
        resp = session.get(url, timeout=60)
        time.sleep(_REQUEST_DELAY)
    except requests.exceptions.RequestException as e:
        logger.debug("RCSB download failed %s: %s", pdb_id, e)
        return None

    if not resp.ok:
        logger.warning("RCSB returned %d for %s", resp.status_code, pdb_id)
        return None

    text = resp.text
    path.write_text(text, encoding="utf-8")
    return text


# ---------------------------------------------------------------------------
# Chain filtering
# ---------------------------------------------------------------------------

def filter_chain_pdb(pdb_text: str, chain: str) -> str:
    """Keep only ATOM records for *chain* (exclude HETATM, alt conformations != ' '/'A')."""
    kept = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            if len(line) <= 21 or line[21] != chain:
                continue
            altloc = line[16] if len(line) > 16 else " "
            if altloc not in (" ", "A"):
                continue
            kept.append(line)
        elif not line.startswith("HETATM"):
            kept.append(line)
    return "\n".join(kept)


def parse_complete_residues(pdb_text: str, chain: str) -> list[tuple[int, str]]:
    """
    Return list of (auth_seq_num, icode) for residues with all 4 backbone atoms
    (N, CA, C, O) present in *chain*.  Ordered by sequence.
    """
    residues: dict[tuple[int, str], set] = OrderedDict()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 27:
            continue
        if line[21] != chain:
            continue
        atom = line[12:16].strip()
        if atom not in _BACKBONE:
            continue
        try:
            seq_num = int(line[22:26])
        except ValueError:
            continue
        icode = line[26]
        residues.setdefault((seq_num, icode), set()).add(atom)
    return [(s, ic) for (s, ic), atoms in residues.items() if len(atoms) == 4]


# ---------------------------------------------------------------------------
# DSSP + motif annotation
# ---------------------------------------------------------------------------

def run_dssp_chain(pdb_text: str, chain: str) -> np.ndarray | None:
    filtered = filter_chain_pdb(pdb_text, chain)
    try:
        coord = pydssp.read_pdbtext(filtered)
        ss    = pydssp.assign(coord, out_type="c3")
        return ss
    except Exception as e:
        logger.debug("pydssp failed: %s", e)
        return None


def annotate_with_pdb(
    uid: str,
    pdb_text: str,
    entry: dict,
    domain_start: int,
    domain_end: int,
) -> list[dict] | None:
    """
    Run DSSP on the PDB chain, map positions to UniProt coordinates, and
    identify β-α motifs within the TIM barrel domain.

    Returns motif list or None on failure.
    """
    chain     = entry["chain_id"]
    unp_start = entry["unp_start"]
    pdb_start = entry["start"]      # auth_seq_num of first residue in mapping
    offset    = unp_start - pdb_start  # auth_seq_num + offset = uniprot_pos

    complete = parse_complete_residues(pdb_text, chain)
    if not complete:
        logger.debug("%s: no complete residues in chain %s", uid, chain)
        return None

    ss = run_dssp_chain(pdb_text, chain)
    if ss is None:
        return None

    if len(ss) != len(complete):
        logger.debug(
            "%s: pydssp output length %d != complete residues %d",
            uid, len(ss), len(complete),
        )
        # Align by position index anyway — lengths can differ if insertion codes exist
        if abs(len(ss) - len(complete)) > 5:
            return None

    # Build a uniprot_pos → ss_char mapping
    pdb_positions = [seq_num for seq_num, _ in complete]
    uniprot_positions = [s + offset for s in pdb_positions]

    if not uniprot_positions:
        return None

    min_uni = min(uniprot_positions)
    max_uni = max(uniprot_positions)
    arr_len = max_uni - min_uni + 1
    ss_full = np.full(arr_len, "-", dtype="U1")

    for i, (uni_pos, ss_char) in enumerate(zip(uniprot_positions, ss)):
        idx = uni_pos - min_uni
        if 0 <= idx < arr_len:
            ss_full[idx] = ss_char

    # Adjust domain boundaries to this array's coordinate frame
    adj_start = max(domain_start - min_uni + 1, 1)
    adj_end   = min(domain_end   - min_uni + 1, arr_len)

    if adj_start >= adj_end:
        logger.debug("%s: domain outside PDB coverage", uid)
        return None

    # identify_ba_motifs works in 1-based coords relative to the ss array
    # Pass adj_start/adj_end and then shift motif coords back to UniProt space
    motifs_raw = identify_ba_motifs(ss_full, adj_start, adj_end)

    # Shift from local (1-based in ss_full) to UniProt coordinates
    shift = min_uni - 1
    motifs = []
    for m in motifs_raw:
        motifs.append({
            "motif":       m["motif"],
            "start":       m["start"]       + shift,
            "end":         m["end"]         + shift,
            "beta_start":  m["beta_start"]  + shift,
            "beta_end":    m["beta_end"]    + shift,
            "alpha_start": m["alpha_start"] + shift,
            "alpha_end":   m["alpha_end"]   + shift,
        })

    return motifs


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def run(
    conn: sqlite3.Connection,
    cache_dir: Path,
    max_resolution: float,
    rerun: bool,
    limit: int | None,
) -> dict:
    where = "" if rerun else "AND pdb_motif_annotations IS NULL AND pdb_source IS NULL"
    query = f"""
        SELECT uniprot_id, domain_start, domain_end
        FROM tb_canonical_analysis
        WHERE domain_start IS NOT NULL AND domain_end IS NOT NULL
          {where}
        ORDER BY uniprot_id
    """
    if limit:
        query += f" LIMIT {limit}"

    rows = conn.execute(query).fetchall()
    logger.info("Proteins to validate: %d", len(rows))

    session = requests.Session()
    session.headers["User-Agent"] = "TIMBarrelResearch/1.0"

    stats = {
        "no_structure": 0,
        "download_failed": 0,
        "dssp_failed": 0,
        "full_8": 0,
        "partial": 0,
    }

    for i, (uid, domain_start, domain_end) in enumerate(rows, 1):
        # 1. Find best experimental structure
        entry = best_structure(uid, domain_start, domain_end, max_resolution, session)
        if entry is None:
            stats["no_structure"] += 1
            conn.execute(
                "UPDATE tb_canonical_analysis SET pdb_source='no_structure' WHERE uniprot_id=?",
                (uid,),
            )
            if i % 100 == 0:
                conn.commit()
            continue

        pdb_id     = entry["pdb_id"]
        chain_id   = entry["chain_id"]
        resolution = entry.get("resolution", 0.0)

        # 2. Download PDB
        pdb_text = fetch_pdb(pdb_id, cache_dir, session)
        if pdb_text is None:
            stats["download_failed"] += 1
            conn.execute(
                "UPDATE tb_canonical_analysis SET pdb_source='download_failed' WHERE uniprot_id=?",
                (uid,),
            )
            if i % 100 == 0:
                conn.commit()
            continue

        # 3. Run DSSP + motif annotation
        motifs = annotate_with_pdb(uid, pdb_text, entry, domain_start, domain_end)
        source = f"{pdb_id.upper()}_{chain_id}_{resolution:.2f}A"

        if motifs is None:
            stats["dssp_failed"] += 1
            conn.execute("""
                UPDATE tb_canonical_analysis
                SET pdb_source='dssp_failed'
                WHERE uniprot_id=?
            """, (uid,))
        else:
            n = len(motifs)
            if n == 8:
                stats["full_8"] += 1
            else:
                stats["partial"] += 1
            conn.execute("""
                UPDATE tb_canonical_analysis
                SET pdb_motif_annotations=?, pdb_source=?
                WHERE uniprot_id=?
            """, (json.dumps(motifs), source, uid))
            logger.debug("%s: %s chain=%s res=%.2f motifs=%d",
                         uid, pdb_id, chain_id, resolution, n)

        if i % 50 == 0:
            conn.commit()
            logger.info(
                "  %d / %d | full8=%d partial=%d no_struct=%d dssp_fail=%d",
                i, len(rows),
                stats["full_8"], stats["partial"],
                stats["no_structure"], stats["dssp_failed"],
            )

    conn.commit()
    return stats


# ---------------------------------------------------------------------------
# Comparison report
# ---------------------------------------------------------------------------

def print_comparison(conn: sqlite3.Connection) -> None:
    rows = conn.execute("""
        SELECT
            json_array_length(motif_annotations)     AS af_motifs,
            json_array_length(pdb_motif_annotations) AS pdb_motifs,
            pdb_source
        FROM tb_canonical_analysis
        WHERE pdb_motif_annotations IS NOT NULL
    """).fetchall()

    buckets: dict = {}
    for af_n, pdb_n, source in rows:
        key = af_n if af_n is not None else -1
        b = buckets.setdefault(key, {"total": 0, "agree8": 0, "pdb_counts": []})
        b["total"] += 1
        if pdb_n == 8:
            b["agree8"] += 1
        if pdb_n is not None:
            b["pdb_counts"].append(pdb_n)

    print(f"\n{'='*76}")
    print("  AlphaFold DSSP vs. experimental PDB validation")
    print(f"{'='*76}")
    print(f"  {'AF motifs':>10}  {'proteins':>9}  {'PDB=8':>7}  {'PDB=8%':>7}  {'med PDB':>8}")
    print(f"  {'-'*10}  {'-'*9}  {'-'*7}  {'-'*7}  {'-'*8}")
    for key in sorted(buckets):
        b = buckets[key]
        label = str(key) if key >= 0 else "none"
        pct   = 100 * b["agree8"] / b["total"] if b["total"] else 0
        sc    = sorted(b["pdb_counts"])
        med   = sc[len(sc) // 2] if sc else float("nan")
        print(f"  {label:>10}  {b['total']:>9}  {b['agree8']:>7}  {pct:>6.0f}%  {med:>8.1f}")

    total_pdb = conn.execute(
        "SELECT COUNT(*) FROM tb_canonical_analysis WHERE pdb_motif_annotations IS NOT NULL"
    ).fetchone()[0]
    pdb8 = conn.execute(
        "SELECT COUNT(*) FROM tb_canonical_analysis WHERE json_array_length(pdb_motif_annotations)=8"
    ).fetchone()[0]
    no_struct = conn.execute(
        "SELECT COUNT(*) FROM tb_canonical_analysis WHERE pdb_source='no_structure'"
    ).fetchone()[0]
    af8_pdb8 = conn.execute("""
        SELECT COUNT(*) FROM tb_canonical_analysis
        WHERE json_array_length(motif_annotations)=8
          AND json_array_length(pdb_motif_annotations)=8
    """).fetchone()[0]
    af8_total = conn.execute("""
        SELECT COUNT(*) FROM tb_canonical_analysis
        WHERE json_array_length(motif_annotations)=8
          AND pdb_motif_annotations IS NOT NULL
    """).fetchone()[0]

    print(f"{'='*76}")
    print(f"  Proteins with experimental PDB motifs : {total_pdb}")
    print(f"  No acceptable PDB structure           : {no_struct}")
    print(f"  PDB gives full 8 motifs               : {pdb8} / {total_pdb}")
    if af8_total:
        print(f"  AF=8 confirmed by PDB=8               : {af8_pdb8} / {af8_total} "
              f"({100*af8_pdb8/af8_total:.0f}%)")
    print(f"{'='*76}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate TIM barrel motif annotations using experimental PDB structures"
    )
    parser.add_argument("--db",          default=None)
    parser.add_argument("--cache-dir",   default="data/pdb_experimental",
                        help="Directory for cached PDB files")
    parser.add_argument("--resolution",  type=float, default=_DEFAULT_RES,
                        help="Maximum resolution in Angstrom (default: %(default)s)")
    parser.add_argument("--rerun",       action="store_true",
                        help="Re-annotate proteins that already have pdb results")
    parser.add_argument("--limit",       type=int, default=None)
    parser.add_argument("--log-level",   default="INFO")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    ensure_columns(conn)
    run(conn, cache_dir, args.resolution, args.rerun, args.limit)
    print_comparison(conn)

    conn.close()


if __name__ == "__main__":
    main()
