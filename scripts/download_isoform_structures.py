#!/usr/bin/env python3
"""
Download AlphaFold structure predictions for AS isoforms from AlphaFold DB.

For each non-canonical isoform in view_noncanonical, queries the AlphaFold
API (https://alphafold.ebi.ac.uk/api/prediction/{isoform_id}) to obtain the
versioned PDB download URL, then fetches the file.

Also records per-isoform metadata returned by the API (global pLDDT,
pLDDT-fraction breakdown, sequence length) in a TSV log for later use.

Output:
  data/alphafold_isoforms/{isoform_id}.pdb   — downloaded structures
  data/alphafold_isoforms/download_log.tsv   — per-isoform result + metadata

Usage:
    python scripts/download_isoform_structures.py [--delay 0.5] [--force]
"""

import argparse
import json
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config


AF_API  = "https://alphafold.ebi.ac.uk/api/prediction"
OUT_DIR = Path("data/alphafold_isoforms")

LOG_HEADER = "\t".join([
    "isoform_id", "uniprot_id", "status", "af_entry_id",
    "latest_version", "global_plddt",
    "frac_very_low", "frac_low", "frac_confident", "frac_very_high",
    "seq_length", "pdb_file",
])


def query_api(isoform_id):
    """Return parsed JSON list from AF API, or None on 404."""
    url = f"{AF_API}/{isoform_id}"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise


def download_pdb(pdb_url, dest):
    urllib.request.urlretrieve(pdb_url, dest)


def load_isoforms(conn):
    return conn.execute("""
        SELECT DISTINCT i.isoform_id, i.uniprot_id
        FROM   isoforms i
        JOIN   view_noncanonical nc ON nc.isoform_id = i.isoform_id
        WHERE  i.is_canonical = 0
        ORDER  BY i.isoform_id
    """).fetchall()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",    default=None)
    parser.add_argument("--delay", type=float, default=0.4,
                        help="Seconds to wait between API calls (default 0.4)")
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if file already exists")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log_path = OUT_DIR / "download_log.tsv"

    db_path  = args.db or get_config().db_path
    conn     = sqlite3.connect(db_path)
    isoforms = load_isoforms(conn)
    conn.close()

    print(f"AS isoforms to process: {len(isoforms)}")
    print(f"Output directory:       {OUT_DIR}")
    print(f"Request delay:          {args.delay}s\n")

    log_rows = []
    n_ok = n_miss = n_skip = n_err = 0

    for isoform_id, uid in isoforms:
        dest = OUT_DIR / f"{isoform_id}.pdb"

        if dest.exists() and not args.force:
            print(f"  SKIP  {isoform_id}")
            log_rows.append(f"{isoform_id}\t{uid}\tskipped" + "\t" * 8)
            n_skip += 1
            continue

        # Query API
        try:
            data = query_api(isoform_id)
        except Exception as exc:
            print(f"  ERR   {isoform_id}  {exc}")
            log_rows.append(f"{isoform_id}\t{uid}\terror" + "\t" * 8)
            n_err += 1
            time.sleep(args.delay)
            continue

        if data is None:
            print(f"  MISS  {isoform_id}  (not in AlphaFold DB)")
            log_rows.append(f"{isoform_id}\t{uid}\tnot_found" + "\t" * 8)
            n_miss += 1
            time.sleep(args.delay)
            continue

        d = data[0]
        pdb_url     = d.get("pdbUrl", "")
        entry_id    = d.get("entryId", "")
        version     = d.get("latestVersion", "")
        g_plddt     = d.get("globalMetricValue", "")
        frac_vl     = d.get("fractionPlddtVeryLow",  "")
        frac_lo     = d.get("fractionPlddtLow",       "")
        frac_co     = d.get("fractionPlddtConfident", "")
        frac_vh     = d.get("fractionPlddtVeryHigh",  "")
        seq_len     = d.get("sequenceEnd", "")

        # Download PDB
        try:
            download_pdb(pdb_url, dest)
            size_kb = dest.stat().st_size // 1024
            print(f"  OK    {isoform_id}  v{version}  pLDDT={g_plddt:.1f}  ({size_kb} KB)")
            log_rows.append("\t".join(str(x) for x in [
                isoform_id, uid, "ok", entry_id, version,
                g_plddt, frac_vl, frac_lo, frac_co, frac_vh,
                seq_len, dest.name,
            ]))
            n_ok += 1
        except Exception as exc:
            print(f"  ERR   {isoform_id}  download failed: {exc}")
            log_rows.append(f"{isoform_id}\t{uid}\tdownload_error" + "\t" * 8)
            n_err += 1

        time.sleep(args.delay)

    # Write log
    with open(log_path, "w") as f:
        f.write(LOG_HEADER + "\n")
        f.write("\n".join(log_rows) + "\n")

    print(f"\n{'='*52}")
    print(f"  Downloaded:     {n_ok}")
    print(f"  Not in AF DB:   {n_miss}")
    print(f"  Skipped:        {n_skip}")
    print(f"  Errors:         {n_err}")
    print(f"  Log:            {log_path}")


if __name__ == "__main__":
    main()
