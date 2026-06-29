#!/usr/bin/env python3
"""
Build (or rebuild) affected_isoforms.

Primary detection: VSP overlap — any alternative isoform with at least one
UniProt VSP feature overlapping the canonical domain is included.
Fallback for isoforms with no VSP data: sliding window (12.5% <= id < 95%).

After building affected_isoforms, also backfills fragment isoforms
(< 200 aa) that have Ensembl transcripts and pass the strict motif-core
filter into analysis_proteins.
"""

import json
import logging
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.analysis.tim_barrel_alignment import populate_tim_barrel_isoforms
from protein_data_collector.database.connection import get_connection
from protein_data_collector.database.schema import init_db

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Fragment isoform helpers (merged from backfill_fragment_isoforms.py)
# ---------------------------------------------------------------------------

def _motif_core(motifs):
    if not motifs:
        return None, None
    return motifs[0]["beta_start"], motifs[-1]["alpha_end"]


def _build_vsp_domain_events(splice_variants, domain_start, domain_end):
    """Convert raw splice_variants to vsp_domain_events format expected by analysis."""
    events = []
    for sv in splice_variants:
        loc  = sv.get("location", {})
        vs   = loc.get("start", {}).get("value")
        ve   = loc.get("end",   {}).get("value")
        fid  = sv.get("featureId", "")
        if vs is None or ve is None:
            continue
        # Only deletions (Missing) overlap with the domain meaningfully
        ol_s = max(vs, domain_start)
        ol_e = min(ve, domain_end)
        if ol_s > ol_e:
            continue
        events.append({
            "feature_id":       fid,
            "can_start":        vs,
            "can_end":          ve,
            "overlap_start":    ol_s,
            "overlap_end":      ol_e,
        })
    return events


def _backfill_fragment_isoforms(conn: sqlite3.Connection) -> None:
    """
    Add fragment isoforms (< 200 aa) that have Ensembl transcripts and pass
    the strict motif-core filter to analysis_proteins.
    """
    conn.row_factory = sqlite3.Row

    candidates = conn.execute("""
        SELECT i.isoform_id, i.uniprot_id, i.sequence_length,
               i.splice_variants, i.ensembl_transcript_id,
               ca.gene_name, ca.domain_start, ca.domain_end,
               ca.motif_annotations,
               p.protein_name
        FROM   isoforms i
        JOIN   canonical_analysis ca ON ca.uniprot_id = i.uniprot_id
                                     AND ca.domain_index = 1
        LEFT JOIN proteins p         ON p.uniprot_id    = i.uniprot_id
        WHERE  i.is_canonical  = 0
          AND  i.is_fragment   = 1
          AND  i.sequence      IS NOT NULL
          AND  i.ensembl_transcript_id IS NOT NULL
          AND  ca.motif_annotations   IS NOT NULL
    """).fetchall()

    # Already-inserted isoform IDs
    existing = {r[0] for r in conn.execute(
        "SELECT isoform_id FROM analysis_proteins WHERE is_canonical=0"
    ).fetchall()}

    inserted = 0
    skipped  = 0

    for row in candidates:
        iso_id    = row["isoform_id"]
        uid       = row["uniprot_id"]
        gene      = row["gene_name"] or ""
        pname     = row["protein_name"] or ""
        ds        = row["domain_start"]
        de        = row["domain_end"]
        sv_raw    = row["splice_variants"]
        ma_raw    = row["motif_annotations"]

        if iso_id in existing:
            continue
        if not sv_raw:
            continue

        motifs    = json.loads(ma_raw)
        core_s, core_e = _motif_core(motifs)
        if core_s is None:
            continue

        svs = json.loads(sv_raw)
        vsp_events = _build_vsp_domain_events(svs, ds, de)
        if not vsp_events:
            continue

        # Strict filter: at least one VSP boundary within motif core
        core_vsps = [
            v for v in vsp_events
            if (core_s <= v["can_start"] <= core_e) or
               (core_s <= v["can_end"]   <= core_e)
        ]
        if not core_vsps:
            skipped += 1
            continue

        conn.execute("""
            INSERT OR IGNORE INTO analysis_proteins (
                uniprot_id, domain_index, isoform_id,
                gene_name, protein_name,
                is_canonical, canonical_id,
                has_experimental_structure, pdb_accession,
                num_motifs, num_domain_junctions, num_as_variants,
                identity_pct, num_vsp_domain_events,
                domain_start, domain_end,
                exon_annotations, motif_annotations,
                vsp_domain_events
            ) VALUES (
                :uniprot_id, 0, :isoform_id,
                :gene_name, :protein_name,
                0, :canonical_id,
                NULL, NULL,
                NULL, NULL, :num_as_variants,
                NULL, :num_vsp_domain_events,
                NULL, NULL,
                NULL, NULL,
                :vsp_domain_events
            )
        """, {
            "uniprot_id":           uid,
            "isoform_id":           iso_id,
            "gene_name":            gene,
            "protein_name":         pname,
            "canonical_id":         uid,
            "num_as_variants":      len(svs),
            "num_vsp_domain_events": len(core_vsps),
            "vsp_domain_events":    json.dumps(core_vsps),
        })
        inserted += 1
        print(f"  Added {iso_id} ({gene}, {row['sequence_length']} aa,"
              f" {len(core_vsps)} core VSP(s))")

    conn.row_factory = None
    conn.commit()
    logger.info("Fragment isoforms backfill: inserted=%d skipped=%d", inserted, skipped)
    print(f"\nInserted: {inserted}   Skipped (no core VSP): {skipped}")


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run(db_path: str) -> None:
    with get_connection(db_path) as conn:
        init_db(conn)
        inserted, skipped_no_overlap, skipped_fallback = populate_tim_barrel_isoforms(
            conn,
            isoform_table="isoforms",
            output_table="affected_isoforms",
        )

    print(f"\n{'='*55}")
    print(f"  AS-affected domain isoforms   : {inserted}")
    print(f"  Skipped — no domain overlap   : {skipped_no_overlap}")
    print(f"  Skipped — no VSP + low id     : {skipped_fallback}")
    print(f"  Total alternative isoforms    : {inserted + skipped_no_overlap + skipped_fallback}")
    print(f"{'='*55}")

    # Merge B: backfill fragment isoforms into analysis_proteins
    logger.info("=== Backfilling fragment isoforms into analysis_proteins ===")
    frag_conn = sqlite3.connect(db_path)
    _backfill_fragment_isoforms(frag_conn)
    frag_conn.close()
