"""
Shared utility for loading canonical exon-junction positions from Ensembl.

All splice-junction analyses should import load_canonical_junctions() from
this module rather than parsing view_canonical.exon_annotations inline.

Why Ensembl is the authoritative source
----------------------------------------
view_canonical.exon_annotations (sourced from the isoforms table) is
incomplete for ~11 % of canonical proteins: it either misses domain-internal
junctions entirely (e.g. ENGASE missing junctions at positions 346, 381) or
records systematically offset positions from a different transcript version.
ensembl_transcripts.exon_annotations was collected directly from Ensembl and
is more complete and consistent.

Exclusions
----------
13 canonical proteins have no matched Ensembl transcript and are excluded from
all Ensembl-based analyses (see exclusion note in Results.md).  These include
GLB1 (P16278) and IMPDH1 (P20839), which have domain-level AS isoforms, plus
11 proteins with no isoforms at all.
"""

import json


def load_canonical_junctions(conn):
    """
    Return {uniprot_id: sorted list of domain-internal junction positions}.

    Junction position = end residue of each non-last exon (1-indexed,
    same convention as view_canonical.exon_annotations e["end"]).
    Only positions satisfying  domain_start <= j < domain_end  are included.

    Source: ensembl_transcripts, matched to canonical proteins by protein
    sequence identity.  Proteins without a sequence match are omitted.

    Parameters
    ----------
    conn : sqlite3.Connection

    Returns
    -------
    dict[str, list[int]]
    """
    # Build lookup: protein_sequence -> list of junction positions
    by_seq = {}
    for _enst_id, seq, ea in conn.execute(
        "SELECT enst_id, sequence, exon_annotations "
        "FROM   ensembl_transcripts "
        "WHERE  exon_annotations IS NOT NULL AND sequence IS NOT NULL"
    ):
        if seq not in by_seq:
            by_seq[seq] = json.loads(ea)   # already a flat list of ints

    # Load canonical domain boundaries and match to Ensembl junctions.
    # Keyed by (uniprot_id, domain_index) so multi-domain proteins each get
    # the junctions that fall within their own domain range.
    rows = conn.execute("""
        SELECT vc.uniprot_id, vc.domain_index, vc.domain_start, vc.domain_end,
               ca.sequence
        FROM   view_canonical vc
        JOIN   canonical_analysis ca
               ON ca.uniprot_id  = vc.uniprot_id
              AND ca.domain_index = vc.domain_index
        WHERE  ca.sequence       IS NOT NULL
          AND  vc.domain_start   IS NOT NULL
          AND  vc.domain_end     IS NOT NULL
    """).fetchall()

    result = {}
    for uid, didx, ds, de, seq in rows:
        if seq not in by_seq:
            continue                        # no Ensembl match — excluded
        all_jcts = by_seq[seq]
        result[(uid, didx)] = sorted(j for j in all_jcts if ds <= j < de)

    return result


def load_isoform_junctions(conn):
    """
    Return {isoform_id: sorted list of junction positions in isoform protein
    coordinates}, sourced from ensembl_transcripts matched by isoform sequence.

    Used by the transcript-comparison AS junction analysis.  Only isoforms
    that (a) are non-canonical and (b) have a matched Ensembl transcript are
    included.

    Returns
    -------
    dict[str, list[int]]
        Keys are UniProt isoform IDs (e.g. 'Q8NFI3-3').
    """
    by_seq = {}
    for _enst_id, seq, ea in conn.execute(
        "SELECT enst_id, sequence, exon_annotations "
        "FROM   ensembl_transcripts "
        "WHERE  exon_annotations IS NOT NULL AND sequence IS NOT NULL"
    ):
        if seq not in by_seq:
            by_seq[seq] = json.loads(ea)

    rows = conn.execute("""
        SELECT i.isoform_id, i.sequence, i.ensembl_transcript_id
        FROM   isoforms i
        WHERE  i.is_canonical = 0 AND i.sequence IS NOT NULL
    """).fetchall()

    by_id = {}
    for _enst_id, seq, ea in conn.execute(
        "SELECT enst_id, sequence, exon_annotations "
        "FROM   ensembl_transcripts "
        "WHERE  exon_annotations IS NOT NULL"
    ):
        by_id[_enst_id.split('.')[0]] = json.loads(ea)

    result = {}
    for iso_id, seq, enst_hint in rows:
        if seq in by_seq:
            result[iso_id] = sorted(by_seq[seq])
        elif enst_hint:
            base = enst_hint.split('.')[0]
            if base in by_id:
                result[iso_id] = sorted(by_id[base])

    return result
