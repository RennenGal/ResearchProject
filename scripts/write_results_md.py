#!/usr/bin/env python3
"""Generate Results.md from cached raw script outputs (results_raw.json)."""

import json
import re
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config

if "--update-summary" not in sys.argv:
    with open("results_raw.json", encoding="utf-8") as f:
        o = json.load(f)
else:
    o = {}  # not needed for --update-summary mode


# ---------------------------------------------------------------------------
# Live dataset summary — queried from DB so excluded proteins are never counted
# ---------------------------------------------------------------------------

def _is_core_affected(vsps, motifs):
    if not motifs:
        return False
    core_s = motifs[0]["beta_start"]
    core_e = motifs[-1]["alpha_end"]
    return any(
        (core_s <= v["can_start"] <= core_e) or (core_s <= v["can_end"] <= core_e)
        for v in vsps
    )


def compute_dataset_summary():
    conn = sqlite3.connect(get_config().db_path)

    # Use view_canonical (domain_index=1) — same scope as the analysis scripts
    canon_rows = conn.execute("""
        SELECT vc.uniprot_id, ca.motif_annotations
        FROM   view_canonical vc
        JOIN   canonical_analysis ca ON ca.uniprot_id = vc.uniprot_id
                                    AND ca.domain_index = vc.domain_index
        WHERE  vc.domain_index = 1
    """).fetchall()
    n_canonical = len(canon_rows)
    motif_map = {uid: json.loads(mj) for uid, mj in canon_rows if mj}

    all_isos = conn.execute("""
        SELECT nc.isoform_id, nc.uniprot_id, nc.vsp_domain_events
        FROM   view_noncanonical nc
        JOIN   isoforms i ON i.isoform_id = nc.isoform_id
        WHERE  nc.vsp_domain_events != '[]'
          AND  i.sequence IS NOT NULL
    """).fetchall()
    conn.close()

    # Restrict to proteins present in view_canonical
    all_isos = [(id_, uid, vj) for id_, uid, vj in all_isos if uid in motif_map]
    n_iso_primary = len(all_isos)
    uids_primary  = {uid for _, uid, _ in all_isos}

    strict = [
        (id_, uid)
        for id_, uid, vj in all_isos
        if _is_core_affected(json.loads(vj), motif_map[uid])
    ]
    n_iso_strict = len(strict)
    uids_strict  = {uid for _, uid in strict}

    return dict(
        n_canonical    = n_canonical,
        n_iso_primary  = n_iso_primary,
        n_can_primary  = len(uids_primary),
        n_iso_strict   = n_iso_strict,
        n_can_strict   = len(uids_strict),
        n_can_no_strict= n_canonical - len(uids_strict),
        n_removed      = n_iso_primary - n_iso_strict,
    )


def _build_summary_block():
    ds = compute_dataset_summary()
    n_total = ds["n_canonical"] + ds["n_iso_strict"]
    rows = [
        "## Dataset summary",
        "",
        "| | Count |",
        "|---|---|",
        f"| Total proteins in analysis (canonical + isoforms) | **{n_total}** |",
        f"| &ensp;— canonical proteins | {ds['n_canonical']} |",
        f"| &ensp;— AS isoforms (strict) | {ds['n_iso_strict']} |",
        f"| Canonical proteins with ≥ 1 domain-level AS isoform (primary) | {ds['n_can_primary']} |",
        f"| Canonical proteins with ≥ 1 motif-core AS isoform (strict) | **{ds['n_can_strict']}** |",
        f"| Canonical proteins with no motif-core AS isoform (strict) | {ds['n_can_no_strict']} |",
        f"| Total AS isoforms (primary) | {ds['n_iso_primary']} |",
        f"| Total AS isoforms (strict) | **{ds['n_iso_strict']}** |",
        f"| Removed by strict filter | {ds['n_removed']} |",
        "",
    ]
    return "\n".join(rows)


def update_summary_in_place(md_path="Results.md"):
    """Replace the Dataset summary section in an existing Results.md."""
    text = Path(md_path).read_text(encoding="utf-8")
    new_block = _build_summary_block()
    text = re.sub(
        r"## Dataset summary\n.*?(?=\n## |\Z)",
        new_block,
        text,
        flags=re.DOTALL,
    )
    Path(md_path).write_text(text, encoding="utf-8")
    print(f"Dataset summary updated in {md_path}")


# ---------------------------------------------------------------------------
# Early exit for --update-summary mode (no results_raw.json needed)
# ---------------------------------------------------------------------------

if "--update-summary" in sys.argv:
    update_summary_in_place()
    sys.exit(0)


def extract(text, pattern, group=1, default="?"):
    m = re.search(pattern, text)
    return m.group(group).strip() if m else default


def sig(p):
    p = float(p)
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


ej  = o["analyze_exon_junctions"]
je  = o["analyze_junction_enrichment"]
me  = o["analyze_motif_enrichment"]
asj = o["analyze_as_splice_junctions"]
vb  = o["analyze_vsp_boundaries"]
ads = o["analyze_as_domain_structure"]
dl  = o["analyze_domain_length_subgroups"]
vc  = o["analyze_vsp_exon_count"]
jc  = o["analyze_junction_consistency"]

# ── Dataset numbers ──────────────────────────────────────────────────────────
n_canonical  = extract(ej, r"Canonical proteins in canonical_analysis\s*:\s*(\d+)")
n_genes      = extract(ej, r"Distinct gene names.*:\s*(\d+)")
n_with_as    = extract(ej, r"Proteins with AS disrupting the domain\s*:\s*(\d+)")
n_iso        = extract(ej, r"Total AS-affected isoforms.*:\s*(\d+)")
n_junctions  = extract(ej, r"Junctions inside TIM barrel domain\s*:\s*(\d+)")
n_enst       = extract(je, r"Loaded (\d+) proteins")
n_full8      = extract(me, r"(\d+) full K_p=8")
n_partial    = extract(me, r"(\d+) partial")
n_can_ads    = extract(ads, r"Canonical proteins:\s+(\d+)")
n_iso_ads    = extract(ads, r"Isoforms with VSP events:\s*(\d+)")
mean_intact  = extract(ads, r"Mean intact:\s*([\d.]+)")
med_intact   = extract(ads, r"Median:\s*([\d.]+)")
n_as_pairs   = extract(asj, r"AS boundary pairs found:\s*(\d+)")
n_vsp_spans  = extract(vb,  r"VSP spans with domain overlap:\s*(\d+)")
n_vsp_prots  = extract(vb,  r"\((\d+) distinct")
n_kp8_prot   = extract(me, r"Restricting to K_p=8 proteins:\s*(\d+)")
n_kp8_jct    = extract(me, r"Total junctions N = (\d+)")
mean_beta    = extract(ads, r"Mean intact beta-strands:\s*([\d.]+)")
mean_alpha   = extract(ads, r"alpha-helices:\s*([\d.]+)")
plddT_pre    = extract(ads, r"Pre-VSP\s+mean pLDDT:\s*([\d.]+)", default="")
plddT_vsp    = extract(ads, r"VSP reg\.\s+mean pLDDT:\s*([\d.]+)", default="")
plddT_post   = extract(ads, r"Post-VSP\s+mean pLDDT:\s*([\d.]+)", default="")
n_pdb        = extract(ads, r"Single-VSP isoforms with PDB:\s*(\d+)", default="")
vex_mean     = extract(vc,  r"Mean:\s*([\d.]+)")
vex_median   = extract(vc,  r"Median:\s*([\d.]+)")
n_vex        = extract(vc,  r"Loaded (\d+) VSP events")

# ── Junction enrichment (5-cat) ───────────────────────────────────────────
je_rows = []
for line in je.split("\n"):
    m = re.match(
        r"\s+(beta-strand|alpha-helix|α→β loop|β→α loop|Flanking)"
        r"\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\S+)",
        line,
    )
    if m:
        je_rows.append(m.groups())
je_chi2 = extract(je, r"chi2\(4\) = ([\d.]+)")
je_p    = extract(je, r"chi2\(4\) = [\d.]+, p = (\S+)")

# ── AS splice junctions ───────────────────────────────────────────────────
asj_start, asj_end = [], []
in_start = in_end = False
asj_n_start = asj_n_end = "?"
for line in asj.split("\n"):
    if "Transcript start positions" in line and "N =" in line:
        in_start, in_end = True, False
        asj_n_start = extract(line, r"N = (\d+)")
    elif "Transcript end positions" in line and "N =" in line:
        in_start, in_end = False, True
        asj_n_end = extract(line, r"N = (\d+)")
    elif "Pooled" in line or "========" in line and (asj_start or asj_end):
        pass
    m = re.match(
        r"\s+(beta-strand|alpha-helix|α→β loop|β→α loop|Flanking)"
        r"\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\S+)",
        line,
    )
    if m:
        if in_start:
            asj_start.append(m.groups())
        elif in_end:
            asj_end.append(m.groups())

# chi2/p values appear as "Global chi2(4) = X,  p = Y" in order: first=start, second=end
chi2_all = re.findall(r"Global chi2\(4\) = ([\d.]+)", asj)
p_all    = re.findall(r"Global chi2\(4\) = [\d.]+,\s+p = (\S+)", asj)
asj_start_chi2 = chi2_all[0] if len(chi2_all) > 0 else "?"
asj_start_p    = p_all[0]    if len(p_all)    > 0 else "?"
asj_end_chi2   = chi2_all[1] if len(chi2_all) > 1 else "?"
asj_end_p      = p_all[1]    if len(p_all)    > 1 else "?"

# ── VSP boundaries ────────────────────────────────────────────────────────
vb_start, vb_end = [], []
in_start = in_end = False
for line in vb.split("\n"):
    if "10A." in line:
        in_start, in_end = True, False
    elif "10B." in line:
        in_start, in_end = False, True
    elif "Pooled" in line or "VSP residue" in line:
        in_start = in_end = False
    m = re.match(
        r"\s+(beta-strand|alpha-helix|α→β loop|β→α loop|Flanking):"
        r"\s+N=(\d+),\s+f=([\d.]+),\s+pi=([\d.]+),\s+rho=([\d.]+),"
        r"\s+raw p=([\d.]+),\s+BH p=([\d.]+)",
        line,
    )
    if m:
        g = m.groups()
        row = (g[0], g[1], g[2], g[3], g[4], g[5], g[6], sig(g[6]))
        if in_start:
            vb_start.append(row)
        elif in_end:
            vb_end.append(row)

# ── Domain disruption ─────────────────────────────────────────────────────
disruption = []
for line in ads.split("\n"):
    m = re.match(r"\s+Position (\d+):\s+(\d+)/(\d+)\s+\(([\d.]+)%\)", line)
    if m:
        disruption.append(m.groups())

# ── Junction consistency ──────────────────────────────────────────────────
jc_ks  = extract(jc, r"KS:\s*D = ([\d.]+)")
jc_p   = extract(jc, r"KS:\s*D = [\d.]+,\s*p = ([\d.]+)")
jc_pp  = extract(jc, r"Permutation p = ([\d.]+)")
jc_per_elem = []
for line in jc.split("\n"):
    m = re.match(
        r"\s+(alpha|beta|loop):\s+n = \d+,\s+KS D = ([\d.]+),"
        r"\s+p = ([\S]+),\s+perm p = ([\d.]+)",
        line,
    )
    if m:
        elem, d, p, pp = m.groups()
        bh_m = re.search(rf"{elem}:\s+perm p = [\d.]+,\s+BH p = ([\d.]+)", jc)
        bh = bh_m.group(1) if bh_m else "?"
        jc_per_elem.append((elem, d, p, pp, bh))

# ── Motif distribution (from create_analysis_table run earlier) ───────────
motif_dist = {1: 2, 2: 7, 3: 1, 4: 5, 5: 2, 6: 17, 7: 27, 8: 140}

# ── Build Results.md ──────────────────────────────────────────────────────
lines = []
w = lines.append

w("# Results")
w("")
w("*Strict motif-core filter: isoforms are included only if at least one VSP overlaps*")
w("*\\[first\\_beta\\_start, last\\_alpha\\_end\\] of the annotated TIM barrel core.*")
w("")
w("---")
w("")
ds = compute_dataset_summary()
n_total_in_analysis = ds["n_canonical"] + ds["n_iso_strict"]

w("## Dataset summary")
w("")
w("| | Count |")
w("|---|---|")
w(f"| Total proteins in analysis (canonical + isoforms) | **{n_total_in_analysis}** |")
w(f"| &ensp;— canonical proteins | {ds['n_canonical']} |")
w(f"| &ensp;— AS isoforms (strict) | {ds['n_iso_strict']} |")
w(f"| Canonical proteins with ≥ 1 domain-level AS isoform (primary) | {ds['n_can_primary']} |")
w(f"| Canonical proteins with ≥ 1 motif-core AS isoform (strict) | **{ds['n_can_strict']}** |")
w(f"| Canonical proteins with no motif-core AS isoform (strict) | {ds['n_can_no_strict']} |")
w(f"| Total AS isoforms (primary) | {ds['n_iso_primary']} |")
w(f"| Total AS isoforms (strict) | **{ds['n_iso_strict']}** |")
w(f"| Removed by strict filter | {ds['n_removed']} |")
w("")
w("### Motif-count distribution")
w("")
w("| $K_p$ | Domain instances |")
w("|---|---|")
for k, v in motif_dist.items():
    w(f"| {k} | {v} |")
w("")
w("### Ensembl transcript coverage and null")
w("")
w("| | Count |")
w("|---|---|")
w(f"| Ensembl-matched canonical proteins (used in Analyses 1–4) | {n_enst} |")
w(f"| Domain-internal exon junctions (all canonical) | {n_junctions} |")
w("")
w("---")
w("")
w("## Analysis 1 — Canonical exon-junction enrichment")
w("")
w(f"**Dataset:** {n_enst} Ensembl-matched canonical proteins, {n_junctions} domain-internal junctions.")
w(f"Junction-count-weighted null; 5-category τ₅.")
w("")
w(f"**Global:** $\\chi^2(4) = {je_chi2}$, $p = {je_p}$")
w("")
w("| Element | $N_t$ | $f_t$ | $\\pi_t^0$ | $\\rho_t$ | $p$ (raw) | $p$ (BH) | Sig |")
w("|---|---|---|---|---|---|---|---|")
for row in je_rows:
    cat, nt, ft, pi, rho, praw, pbh, s = row
    w(f"| {cat} | {nt} | {ft} | {pi} | **{rho}** | {praw} | {pbh} | {s} |")
w("")
w("![Exon junction enrichment in TIM-barrel structural elements](figures/enrichment_bars.png)")
w("")
w("---")
w("")
w("## Analysis 2 — Motif-specific junction enrichment ($K_p = 8$ proteins)")
w("")
w(f"**Dataset:** {n_kp8_prot} proteins with full 8-motif annotation, {n_kp8_jct} domain-internal junctions.")
w("")
w("No individual (element, position) category survives BH correction at 5%.")
w("The global test is significant but the signal is distributed rather than concentrated in any single position.")
w("")
w("![Motif-specific junction enrichment heatmap](figures/motif_enrichment_heatmap.png)")
w("")
w("---")
w("")
w("## Analysis 3 — Transcript-derived AS boundary enrichment")
w("")
w(f"**Dataset:** {n_enst} Ensembl-matched canonical proteins; {n_iso} strict isoforms; **{n_as_pairs} AS boundary pairs**.")
w("")
if asj_start:
    w(
        f"**Transcript start positions ($D_{{\\text{{seq}}}}$),"
        f" $N = {asj_n_start}$,"
        f" global $\\chi^2(4) = {asj_start_chi2}$, $p = {asj_start_p}$:**"
    )
    w("")
    w("| Element | $N_t$ | $f_t$ | $\\pi_t^0$ | $\\rho_t$ | $p$ (raw) | $p$ (BH) | Sig |")
    w("|---|---|---|---|---|---|---|---|")
    for row in asj_start:
        cat, nt, ft, pi, rho, praw, pbh, s = row
        w(f"| {cat} | {nt} | {ft} | {pi} | **{rho}** | {praw} | {pbh} | {s} |")
    w("")
if asj_end:
    w(
        f"**Transcript end positions ($R_{{\\text{{can}}}}-1$),"
        f" $N = {asj_n_end}$,"
        f" global $\\chi^2(4) = {asj_end_chi2}$, $p = {asj_end_p}$:**"
    )
    w("")
    w("| Element | $N_t$ | $f_t$ | $\\pi_t^0$ | $\\rho_t$ | $p$ (raw) | $p$ (BH) | Sig |")
    w("|---|---|---|---|---|---|---|---|")
    for row in asj_end:
        cat, nt, ft, pi, rho, praw, pbh, s = row
        w(f"| {cat} | {nt} | {ft} | {pi} | **{rho}** | {praw} | {pbh} | {s} |")
    w("")
w("![AS boundary enrichment — start positions](figures/as_splice_junctions.png)")
w("")
w("![AS boundary enrichment — pooled](figures/as_splice_junctions_pooled.png)")
w("")
w("---")
w("")
w("## Analysis 4 — VSP boundary placement in structural elements")
w("")
w(f"**Dataset:** {n_vsp_spans} VSP spans across {n_vsp_prots} canonical proteins.")
w("")
if vb_start:
    w("**VSP start positions:**")
    w("")
    w("| Element | $N_t$ | $f_t$ | $\\pi_t^0$ | $\\rho_t$ | $p$ (raw) | $p$ (BH) | Sig |")
    w("|---|---|---|---|---|---|---|---|")
    for row in vb_start:
        cat, nt, ft, pi, rho, praw, pbh, s = row
        w(f"| {cat} | {nt} | {ft} | {pi} | **{rho}** | {praw} | {pbh} | {s} |")
    w("")
if vb_end:
    w("**VSP end positions:**")
    w("")
    w("| Element | $N_t$ | $f_t$ | $\\pi_t^0$ | $\\rho_t$ | $p$ (raw) | $p$ (BH) | Sig |")
    w("|---|---|---|---|---|---|---|---|")
    for row in vb_end:
        cat, nt, ft, pi, rho, praw, pbh, s = row
        w(f"| {cat} | {nt} | {ft} | {pi} | **{rho}** | {praw} | {pbh} | {s} |")
    w("")
w("![VSP start enrichment](figures/vsp_start_enrichment.png)")
w("")
w("![VSP end enrichment](figures/vsp_end_enrichment.png)")
w("")
w("![VSP boundary enrichment combined](figures/vsp_boundary_enrichment.png)")
w("")
w("![VSP boundary enrichment pooled](figures/vsp_boundary_pooled.png)")
w("")
w("![VSP residue coverage](figures/vsp_residue_coverage.png)")
w("")
w("---")
w("")
w("## Analysis 5 — Structural impact on barrel architecture")
w("")
w(f"**Dataset:** {n_can_ads} canonical proteins, {n_iso_ads} AS isoforms (strict filter).")
w("")
w(
    f"**Mean intact motifs:** {mean_intact} (median {med_intact})."
    " No isoform retains all 8 motifs."
)
w("")
w("### Per-position disruption rate (combined β+α)")
w("")
w("| Position | Disrupted | Total | Rate |")
w("|---|---|---|---|")
for pos, dis, tot, rate in disruption:
    w(f"| {pos} | {dis} | {tot} | {rate}% |")
w("")
w(
    f"Separate β/α analysis: mean intact β-strands = {mean_beta},"
    f" mean intact α-helices = {mean_alpha}."
)
w("")
w("![Domain disruption combined](figures/as_domain_disruption.png)")
w("")
w("![Domain disruption separate β/α](figures/as_domain_disruption_separate.png)")
w("")
w("![Motif disruption heatmap](figures/motif_disruption_heatmap.png)")
w("")
w("![Isoform disruption heatmap](figures/isoform_disruption_heatmap.png)")
w("")
w("![Canonical K_p distribution](figures/canonical_kp_distribution.png)")
w("")
if plddT_pre:
    w("### AlphaFold pLDDT by region")
    w("")
    w(f"Single-VSP isoforms with AlphaFold structure: {n_pdb}.")
    w("")
    w("| Region | Mean pLDDT |")
    w("|---|---|")
    w(f"| Pre-VSP | {plddT_pre} |")
    w(f"| VSP region | {plddT_vsp} |")
    w(f"| Post-VSP | {plddT_post} |")
    w("")
w("---")
w("")
w("## Analysis 6 — Domain-length subgroup junction enrichment")
w("")
dl_sub = extract(dl, r"Proteins with domain.*:\s*(\d+)")
w(
    f"**{dl_sub} proteins** with domain 200–400 aa and full exon+motif annotations,"
    " analysed across four 50 aa bins."
)
w("")
w(
    "Loop depletion is consistent across all four subgroups;"
    " alpha-helix enrichment strengthens with domain length."
)
w("")
w("---")
w("")
w("## Analysis 7 — Exons per VSP event")
w("")
w(
    f"**{n_vex} VSP events** analysed."
    f" Mean domain exons per VSP: {vex_mean} (median {vex_median})."
)
w("")
w("![Domain exons per VSP event](figures/vsp_exon_count.png)")
w("")
w("---")
w("")
w("## Analysis 8 — Junction positional consistency")
w("")
w(
    f"**Global KS test (all elements combined):** $D = {jc_ks}$, $p = {jc_p}$"
    f" (permutation $p = {jc_pp}$) — no global positional clustering."
)
w("")
w("| Element | KS $D$ | KS $p$ | Permutation $p$ | BH $p$ |")
w("|---|---|---|---|---|")
for elem, d, p, pp, bh in jc_per_elem:
    w(f"| {elem} | {d} | {p} | {pp} | {bh} |")
w("")
w("![Junction positional consistency — global](figures/consistency_global.png)")
w("")
w("![Junction phase within elements](figures/consistency_phase.png)")
w("")

Path("Results.md").write_text("\n".join(lines), encoding="utf-8")
print(f"Results.md written ({len(lines)} lines)")
