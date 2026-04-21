"""Export query results to FASTA, CSV, and JSON."""

import csv
import io
import json
from typing import Any, Dict, List


def to_fasta(isoforms: List[Dict[str, Any]]) -> str:
    """
    Convert isoform rows to FASTA format.

    Header: >{isoform_id} {uniprot_id} canonical={0|1}
    """
    lines = []
    for iso in isoforms:
        header = (
            f">{iso['isoform_id']} "
            f"uniprot_id={iso['uniprot_id']} "
            f"canonical={iso.get('is_canonical', 0)}"
        )
        seq = iso.get("sequence", "")
        lines.append(header)
        # Wrap at 60 characters (standard FASTA)
        for i in range(0, len(seq), 60):
            lines.append(seq[i : i + 60])
    return "\n".join(lines)


def to_csv(isoforms: List[Dict[str, Any]]) -> str:
    """
    Convert isoform rows to CSV (excludes sequence to keep output readable).

    Includes: isoform_id, uniprot_id, is_canonical, sequence_length,
              exon_count, ensembl_transcript_id, alphafold_id,
              has_tim_barrel, has_splice_variants
    """
    fields = [
        "isoform_id", "uniprot_id", "is_canonical", "sequence_length",
        "exon_count", "ensembl_transcript_id", "alphafold_id",
        "has_tim_barrel", "has_splice_variants",
    ]
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    for iso in isoforms:
        row = {k: iso.get(k) for k in fields}
        row["has_tim_barrel"] = 1 if iso.get("tim_barrel_location") else 0
        sv = iso.get("splice_variants")
        row["has_splice_variants"] = 1 if sv else 0
        writer.writerow(row)
    return buf.getvalue()


def to_json(data: Any, indent: int = 2) -> str:
    """Serialize *data* to a JSON string."""
    return json.dumps(data, indent=indent, default=str)


def write_fasta(isoforms: List[Dict[str, Any]], path: str) -> None:
    with open(path, "w") as fh:
        fh.write(to_fasta(isoforms))


def write_csv(isoforms: List[Dict[str, Any]], path: str) -> None:
    with open(path, "w", newline="") as fh:
        fh.write(to_csv(isoforms))


def write_json(data: Any, path: str) -> None:
    with open(path, "w") as fh:
        fh.write(to_json(data))
