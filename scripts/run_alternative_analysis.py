#!/usr/bin/env python3
"""
Alternative analysis — strict motif-core isoform filter.

Restricts the study group to isoforms where at least one VSP overlaps
[first_motif_beta_start, last_motif_alpha_end] — the annotated motif core
of the TIM barrel — rather than just the Gene3D/CATH domain boundaries.

Figures written to figures/alt/.
All statistics printed to stdout for capture into alternative-results.md.

Usage:
    python scripts/run_alternative_analysis.py
"""

import json
import sqlite3
import sys

sys.stdout.reconfigure(encoding="utf-8")

from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as scipy_stats
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config
from junction_utils import load_canonical_junctions

OUTDIR = Path("figures/alt")
CATS5  = ["beta", "alpha", "inter", "loop", "flanking"]
LABELS5 = {"beta": "β-strand", "alpha": "α-helix", "inter": "α→β loop",
            "loop": "β→α loop", "flanking": "Flanking"}
COLS5  = {"beta": "#4C72B0", "alpha": "#DD8452", "inter": "#55A868",
          "loop": "#8c8c8c", "flanking": "#C44E52"}

# Motif-core 4-category version (no flanking)
CATS4   = ["beta", "loop", "alpha", "inter"]
LABELS4 = {"beta": "β-strand", "alpha": "α-helix",
           "inter": "α→β loop", "loop": "β→α loop"}

# Analysis 2 motif-specific category encoding (mirrors analyze_motif_enrichment.py)
_ELEM_TYPES  = ["beta", "loop", "alpha"]
_N_MOTIFS    = 8
_INTER_BASE  = 24
_FLANKING    = 31
_BETA_LOOP_ALPHA = [(t, k) for k in range(1, _N_MOTIFS + 1) for t in _ELEM_TYPES]
_INTER_CATS  = [("inter", k) for k in range(1, _N_MOTIFS)]
PRIMARY_CATS = _BETA_LOOP_ALPHA + _INTER_CATS   # 31 categories
_CAT_IDX = {(t, k): 3 * (k - 1) + i
            for k in range(1, _N_MOTIFS + 1) for i, t in enumerate(_ELEM_TYPES)}
_CAT_IDX.update({("inter", k): _INTER_BASE + (k - 1) for k in range(1, _N_MOTIFS)})
_IDX_CAT = {v: k for k, v in _CAT_IDX.items()}

MIN_MATCH = 15
MAX_SLIDE = 5


# ---------------------------------------------------------------------------
# Strict filter
# ---------------------------------------------------------------------------

def is_core_affected(vsps, motifs):
    """True if at least one VSP has its start or end WITHIN [first beta_start, last alpha_end].

    A VSP that merely spans the entire core (starts before AND ends after) is excluded:
    it deletes the whole barrel rather than splicing at a specific structural position.
    """
    if not motifs:
        return False
    core_s = motifs[0]["beta_start"]
    core_e = motifs[-1]["alpha_end"]
    return any(
        (core_s <= v["can_start"] <= core_e) or (core_s <= v["can_end"] <= core_e)
        for v in vsps
    )


# ---------------------------------------------------------------------------
# Structural element classifier (tau5)
# ---------------------------------------------------------------------------

def _tau5(pos, motifs):
    if not motifs or pos < motifs[0]["beta_start"]:
        return "flanking"
    for i, m in enumerate(motifs):
        if m["beta_start"] <= pos <= m["beta_end"]:
            return "beta"
        if m["beta_end"] < pos < m["alpha_start"]:
            return "loop"
        if m["alpha_start"] <= pos <= m["alpha_end"]:
            return "alpha"
        if i + 1 < len(motifs) and m["alpha_end"] < pos < motifs[i + 1]["beta_start"]:
            return "inter"
    return "flanking"


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def bh_correct(pvals_dict):
    keys  = list(pvals_dict)
    pvals = np.array([pvals_dict[k] for k in keys])
    order = np.argsort(pvals)
    m     = len(pvals)
    adj   = np.zeros(m)
    for rank, idx in enumerate(order, 1):
        adj[idx] = pvals[idx] * m / rank
    for j in range(m - 2, -1, -1):
        adj[order[j]] = min(adj[order[j]], adj[order[j + 1]])
    adj = np.minimum(adj, 1.0)
    return {k: float(adj[i]) for i, k in enumerate(keys)}


def chi_square_enrichment(counts, N, pi_null):
    """Return (f_t, rho_t, pvals_raw, pvals_bh, z_scores) all keyed by CATS5."""
    f_t, rho_t, pvals_raw, z_scores = {}, {}, {}, {}
    for t in CATS5:
        f_t[t]   = counts[t] / N if N > 0 else 0.0
        rho_t[t] = f_t[t] / pi_null[t] if pi_null[t] > 0 else float("nan")
        E_t = N * pi_null[t]
        if E_t > 0:
            z = (counts[t] - E_t) / np.sqrt(E_t)
            z_scores[t]  = float(z)
            pvals_raw[t] = float(2 * scipy_stats.norm.sf(abs(z)))
        else:
            z_scores[t]  = float("nan")
            pvals_raw[t] = 1.0
    pvals_bh = bh_correct(pvals_raw)
    return f_t, rho_t, pvals_raw, pvals_bh, z_scores


def bh_ci_halfwidths(pvals_raw, pi_null, N, alpha=0.05):
    m     = len(CATS5)
    order = sorted(CATS5, key=lambda t: pvals_raw[t])
    ranks = {t: r for r, t in enumerate(order, 1)}
    errs  = {}
    for t in CATS5:
        alpha_eff = alpha * ranks[t] / m
        z_t = float(scipy_stats.norm.ppf(1 - alpha_eff / 2))
        E_t = N * pi_null[t]
        errs[t] = z_t / np.sqrt(E_t) if E_t > 0 else 0.0
    return errs


def sig_stars(p):
    if p < 0.05: return "*"
    return "ns"


def print_table(label, counts, N, pi_null, f_t, rho_t, pvals_raw, pvals_bh, chi2, p_global):
    print(f"\n{'='*72}")
    print(f"  {label}  (N = {N})")
    print(f"  Global chi2(4) = {chi2:.2f},  p = {p_global:.4g}")
    print(f"{'='*72}")
    print(f"  {'Element':<16}  {'N_t':>5}  {'f_t':>7}  {'pi_t':>7}  {'rho':>6}  {'p_raw':>8}  {'p_BH':>8}  Sig")
    print("  " + "-"*80)
    for t in CATS5:
        print(f"  {LABELS5[t]:<16}  {counts[t]:>5}  {f_t[t]:>7.3f}  "
              f"{pi_null[t]:>7.3f}  {rho_t[t]:>6.3f}  "
              f"{pvals_raw[t]:>8.4f}  {pvals_bh[t]:>8.4f}  {sig_stars(pvals_bh[t])}")


def bh_ci_halfwidths_4cat(pvals_raw, pi_null, N, alpha=0.05):
    m     = len(CATS4)
    order = sorted(CATS4, key=lambda t: pvals_raw[t])
    ranks = {t: r for r, t in enumerate(order, 1)}
    errs  = {}
    for t in CATS4:
        alpha_eff = alpha * ranks[t] / m
        z_t = float(scipy_stats.norm.ppf(1 - alpha_eff / 2))
        E_t = N * pi_null[t]
        errs[t] = z_t / np.sqrt(E_t) if E_t > 0 else 0.0
    return errs


def chi_square_enrichment_4cat(counts, N, pi_null):
    """Return (f_t, rho_t, pvals_raw, pvals_bh, z_scores) all keyed by CATS4."""
    f_t, rho_t, pvals_raw, z_scores = {}, {}, {}, {}
    for t in CATS4:
        f_t[t]   = counts[t] / N if N > 0 else 0.0
        rho_t[t] = f_t[t] / pi_null[t] if pi_null[t] > 0 else float("nan")
        E_t = N * pi_null[t]
        if E_t > 0:
            z = (counts[t] - E_t) / np.sqrt(E_t)
            z_scores[t]  = float(z)
            pvals_raw[t] = float(2 * scipy_stats.norm.sf(abs(z)))
        else:
            z_scores[t]  = float("nan")
            pvals_raw[t] = 1.0
    pvals_bh = bh_correct(pvals_raw)
    return f_t, rho_t, pvals_raw, pvals_bh, z_scores


def print_table_4cat(label, counts, N, pi_null, f_t, rho_t, pvals_raw, pvals_bh, chi2, p_global):
    print(f"\n{'='*72}")
    print(f"  {label}  (N = {N})")
    print(f"  Global chi2(3) = {chi2:.2f},  p = {p_global:.4g}")
    print(f"{'='*72}")
    print(f"  {'Element':<16}  {'N_t':>5}  {'f_t':>7}  {'pi_t':>7}  {'rho':>6}  {'p_raw':>8}  {'p_BH':>8}  Sig")
    print("  " + "-"*80)
    for t in CATS4:
        print(f"  {LABELS4[t]:<16}  {counts[t]:>5}  {f_t[t]:>7.3f}  "
              f"{pi_null[t]:>7.3f}  {rho_t[t]:>6.3f}  "
              f"{pvals_raw[t]:>8.4f}  {pvals_bh[t]:>8.4f}  {sig_stars(pvals_bh[t])}")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_canonical(conn):
    # domain_index = 1 picks the primary (first/largest) domain for multi-domain proteins.
    rows = conn.execute("""
        SELECT vc.uniprot_id, vc.gene_name, vc.domain_start, vc.domain_end,
               vc.motif_annotations, ca.sequence
        FROM   view_canonical vc
        JOIN   canonical_analysis ca ON ca.uniprot_id  = vc.uniprot_id
                                    AND ca.domain_index = vc.domain_index
        WHERE  vc.domain_index = 1
    """).fetchall()
    out = {}
    for uid, gene, ds, de, mj, seq in rows:
        out[uid] = dict(gene=gene or uid, ds=ds, de=de,
                        motifs=json.loads(mj), seq=seq)
    return out


def load_all_isoforms(conn):
    rows = conn.execute("""
        SELECT nc.isoform_id, nc.uniprot_id, nc.vsp_domain_events,
               i.sequence_length, i.sequence
        FROM   view_noncanonical nc
        JOIN   isoforms i ON i.isoform_id = nc.isoform_id
        WHERE  nc.vsp_domain_events != '[]'
          AND  i.sequence IS NOT NULL
    """).fetchall()
    out = []
    for iso_id, uid, vsp_json, iso_len, seq in rows:
        out.append(dict(isoform_id=iso_id, uniprot_id=uid,
                        vsps=json.loads(vsp_json), iso_len=iso_len, seq=seq))
    return out


# ---------------------------------------------------------------------------
# Apply strict filter
# ---------------------------------------------------------------------------

def apply_strict_filter(isoforms, canonicals):
    """Keep only isoforms with ≥1 VSP overlapping the motif core."""
    return [iso for iso in isoforms
            if iso["uniprot_id"] in canonicals
            and is_core_affected(iso["vsps"], canonicals[iso["uniprot_id"]]["motifs"])]


# ---------------------------------------------------------------------------
# Analysis 5 helpers — motif disruption classification
# ---------------------------------------------------------------------------

def classify_element(elem_start, elem_end, vsps):
    state = "intact"
    for vsp in vsps:
        vs, ve = vsp["can_start"], vsp["can_end"]
        if ve < elem_start or vs > elem_end:
            continue
        if vs <= elem_start and ve >= elem_end:
            return "removed"
        state = "partial"
    return state


def _overlap_fraction(elem_start, elem_end, vsps):
    """Fraction of motif residues [elem_start, elem_end] covered by VSPs (merged)."""
    motif_len = elem_end - elem_start + 1
    if motif_len <= 0:
        return 0.0
    intervals = []
    for vsp in vsps:
        lo = max(vsp["can_start"], elem_start)
        hi = min(vsp["can_end"],   elem_end)
        if lo <= hi:
            intervals.append((lo, hi))
    if not intervals:
        return 0.0
    intervals.sort()
    merged = [list(intervals[0])]
    for lo, hi in intervals[1:]:
        if lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    covered = sum(hi - lo + 1 for lo, hi in merged)
    return covered / motif_len


# 4-state encoding: 0=intact, 1=partial_minor(<50%), 2=partial_major(>50%), 3=deleted
_DETAIL_ENC = {"intact": 0, "partial_minor": 1, "partial_major": 2, "deleted": 3}


def classify_element_detailed(elem_start, elem_end, vsps):
    frac = _overlap_fraction(elem_start, elem_end, vsps)
    if frac == 0.0:
        return "intact"
    if frac >= 1.0:
        return "deleted"
    if frac > 0.5:
        return "partial_major"
    return "partial_minor"


def classify_combined_detailed(iso, canonicals):
    uid = iso["uniprot_id"]
    return [classify_element_detailed(m["beta_start"], m["alpha_end"], iso["vsps"])
            for m in canonicals[uid]["motifs"]]


def classify_combined(iso, canonicals):
    uid = iso["uniprot_id"]
    return [classify_element(m["beta_start"], m["alpha_end"], iso["vsps"])
            for m in canonicals[uid]["motifs"]]


def classify_separate(iso, canonicals):
    uid = iso["uniprot_id"]
    return [
        (classify_element(m["beta_start"], m["beta_end"],   iso["vsps"]),
         classify_element(m["alpha_start"], m["alpha_end"], iso["vsps"]))
        for m in canonicals[uid]["motifs"]
    ]


# ---------------------------------------------------------------------------
# Analysis 5 figures
# ---------------------------------------------------------------------------

def plot_disruption_combined(combined_list, out):
    n = len(combined_list)
    intact_counts = [sum(1 for s in st if s == "intact") for st in combined_list]
    max_pos = max(len(st) for st in combined_list)
    pos_rate = []
    for i in range(max_pos):
        relevant = [st for st in combined_list if len(st) > i]
        pos_rate.append(100 * sum(1 for st in relevant if st[i] != "intact") / len(relevant))

    fig, ax = plt.subplots(figsize=(6, 4.5))
    max_count = max(intact_counts)
    ax.hist(intact_counts, bins=np.arange(-0.5, max_count + 1.5, 1),
            color="#4C72B0", alpha=0.85, edgecolor="white", zorder=3)
    ax.set_xlabel("Intact β/α motifs retained", fontsize=15)
    ax.set_ylabel("Number of isoforms", fontsize=15)
    ax.set_xticks(range(max_count + 1))
    ax.set_title(f"Intact β/α motifs per isoform\n{n} isoforms", fontsize=16)
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")
    return intact_counts, pos_rate


def plot_isoform_heatmap(isoforms, combined_map, canonicals, out):
    N_POS = 8
    iso_ids  = [iso["isoform_id"] for iso in isoforms if iso["isoform_id"] in combined_map]
    iso_uids = {iso["isoform_id"]: iso["uniprot_id"] for iso in isoforms}
    n = len(iso_ids)

    mat = np.full((n, N_POS), np.nan)
    for row, iso_id in enumerate(iso_ids):
        states = combined_map[iso_id]
        uid    = iso_uids[iso_id]
        kp     = len(canonicals[uid]["motifs"])
        for col in range(min(kp, N_POS)):
            mat[row, col] = 0.0 if states[col] == "intact" else 1.0

    mat_filled = np.where(np.isnan(mat), 0.0, mat)
    if n > 1:
        Z     = linkage(pdist(mat_filled, metric="euclidean"), method="ward")
        order = leaves_list(Z)
    else:
        order = [0]
    mat     = mat[order]
    iso_ids = [iso_ids[i] for i in order]

    gene_count = {}
    for iso_id in iso_ids:
        g = canonicals[iso_uids[iso_id]]["gene"]
        gene_count[g] = gene_count.get(g, 0) + 1
    gene_seen = {}
    labels = []
    for iso_id in iso_ids:
        g = canonicals[iso_uids[iso_id]]["gene"]
        if gene_count[g] > 1:
            gene_seen[g] = gene_seen.get(g, 0) + 1
            labels.append(f"{g}-{gene_seen[g]}")
        else:
            labels.append(g)

    fig_h = max(6, n * 0.18)
    fig, ax = plt.subplots(figsize=(7, fig_h))
    cmap = matplotlib.colors.ListedColormap(["#f7f7f7", "#C44E52"])
    cmap.set_bad(color="#d0d0d0")
    masked = np.ma.masked_invalid(mat)
    im = ax.imshow(masked, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")

    ax.set_xticks(range(N_POS))
    ax.set_xticklabels([f"$\\beta\\alpha${k}" for k in range(1, N_POS + 1)], fontsize=14)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel("Barrel motif position", fontsize=15)
    ax.set_title(
        f"Per-isoform motif disruption vs. canonical  ({n} isoforms, strict filter)\n"
        "red = disrupted (partial or removed);  rows clustered by Ward linkage",
        fontsize=16,
    )
    ax.set_xticks(np.arange(-0.5, N_POS, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_ticks([0.25, 0.75])
    cbar.set_ticklabels(["intact", "disrupted"], fontsize=12)
    cbar.ax.tick_params(length=0)

    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def plot_isoform_heatmap_by_gene(isoforms, combined_map, canonicals, out):
    """Per-isoform disruption heatmap sorted alphabetically by gene name."""
    N_POS   = 8
    iso_ids  = [iso["isoform_id"] for iso in isoforms if iso["isoform_id"] in combined_map]
    iso_uids = {iso["isoform_id"]: iso["uniprot_id"] for iso in isoforms}

    # Sort by (gene name, isoform_id) for stable alphabetical ordering
    iso_ids.sort(key=lambda iso_id: (
        canonicals[iso_uids[iso_id]]["gene"],
        iso_id,
    ))
    n = len(iso_ids)

    mat = np.full((n, N_POS), np.nan)
    for row, iso_id in enumerate(iso_ids):
        states = combined_map[iso_id]
        uid    = iso_uids[iso_id]
        kp     = len(canonicals[uid]["motifs"])
        for col in range(min(kp, N_POS)):
            mat[row, col] = 0.0 if states[col] == "intact" else 1.0

    # Build row labels and track gene boundaries for separator lines
    gene_count = {}
    for iso_id in iso_ids:
        g = canonicals[iso_uids[iso_id]]["gene"]
        gene_count[g] = gene_count.get(g, 0) + 1
    gene_seen   = {}
    labels      = []
    separators  = []   # row indices AFTER which to draw a line
    prev_gene   = None
    for row, iso_id in enumerate(iso_ids):
        g = canonicals[iso_uids[iso_id]]["gene"]
        if g != prev_gene and prev_gene is not None:
            separators.append(row - 0.5)
        prev_gene = g
        if gene_count[g] > 1:
            gene_seen[g] = gene_seen.get(g, 0) + 1
            labels.append(f"{g}-{gene_seen[g]}")
        else:
            labels.append(g)

    fig_h = max(6, n * 0.18)
    fig, ax = plt.subplots(figsize=(7, fig_h))
    cmap = matplotlib.colors.ListedColormap(["#f7f7f7", "#C44E52"])
    cmap.set_bad(color="#d0d0d0")
    masked = np.ma.masked_invalid(mat)
    im = ax.imshow(masked, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")

    for y in separators:
        ax.axhline(y, color="#888888", linewidth=0.6, linestyle="--")

    ax.set_xticks(range(N_POS))
    ax.set_xticklabels([f"$\\beta\\alpha${k}" for k in range(1, N_POS + 1)], fontsize=14)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel("Barrel motif position", fontsize=15)
    ax.set_title(
        f"Per-isoform motif disruption vs. canonical  ({n} isoforms, strict filter)\n"
        "red = disrupted (partial or removed);  rows sorted by gene name",
        fontsize=16,
    )
    ax.set_xticks(np.arange(-0.5, N_POS, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_ticks([0.25, 0.75])
    cbar.set_ticklabels(["intact", "disrupted"], fontsize=12)
    cbar.ax.tick_params(length=0)

    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def plot_isoform_heatmap_detailed(isoforms, detailed_map, canonicals, out):
    """Per-isoform heatmap with 4 disruption states."""
    N_POS   = 8
    iso_ids  = [iso["isoform_id"] for iso in isoforms if iso["isoform_id"] in detailed_map]
    iso_uids = {iso["isoform_id"]: iso["uniprot_id"] for iso in isoforms}
    n = len(iso_ids)

    mat = np.full((n, N_POS), np.nan)
    for row, iso_id in enumerate(iso_ids):
        states = detailed_map[iso_id]
        uid    = iso_uids[iso_id]
        kp     = len(canonicals[uid]["motifs"])
        for col in range(min(kp, N_POS)):
            mat[row, col] = float(_DETAIL_ENC[states[col]])

    mat_filled = np.where(np.isnan(mat), 0.0, mat)
    if n > 1:
        Z     = linkage(pdist(mat_filled, metric="euclidean"), method="ward")
        order = leaves_list(Z)
    else:
        order = [0]
    mat     = mat[order]
    iso_ids = [iso_ids[i] for i in order]

    gene_count = {}
    for iso_id in iso_ids:
        g = canonicals[iso_uids[iso_id]]["gene"]
        gene_count[g] = gene_count.get(g, 0) + 1
    gene_seen = {}
    labels = []
    for iso_id in iso_ids:
        g = canonicals[iso_uids[iso_id]]["gene"]
        if gene_count[g] > 1:
            gene_seen[g] = gene_seen.get(g, 0) + 1
            labels.append(f"{g}-{gene_seen[g]}")
        else:
            labels.append(g)

    cmap   = matplotlib.colors.ListedColormap(["#f7f7f7", "#FDB96D", "#C44E52", "#2D2D2D"])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm   = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    cmap.set_bad(color="#d0d0d0")

    fig_h = max(6, n * 0.18)
    fig, ax = plt.subplots(figsize=(7, fig_h))
    masked = np.ma.masked_invalid(mat)
    im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")

    ax.set_xticks(range(N_POS))
    ax.set_xticklabels([f"$\\beta\\alpha${k}" for k in range(1, N_POS + 1)], fontsize=14)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel("Barrel motif position", fontsize=15)
    ax.set_title(
        f"Per-isoform motif disruption — 4-state ({n} isoforms, strict filter)\n"
        "rows clustered by Ward linkage",
        fontsize=16,
    )
    ax.set_xticks(np.arange(-0.5, N_POS, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", length=0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, ticks=[0, 1, 2, 3])
    cbar.set_ticklabels(["intact", "<50% deleted", ">50% deleted", "deleted"], fontsize=12)
    cbar.ax.tick_params(length=0)

    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def plot_motif_heatmap(isoforms, combined_map, canonicals, out):
    uid_to_states = {}
    for iso in isoforms:
        uid    = iso["uniprot_id"]
        iso_id = iso["isoform_id"]
        if iso_id not in combined_map:
            continue
        uid_to_states.setdefault(uid, []).append(combined_map[iso_id])

    N_POS = 8
    uids  = sorted(uid_to_states)
    mat   = np.full((len(uids), N_POS), np.nan)
    for row, uid in enumerate(uids):
        kp = len(canonicals[uid]["motifs"])
        for col in range(min(kp, N_POS)):
            states  = uid_to_states[uid]
            relevant = [st for st in states if len(st) > col]
            if relevant:
                mat[row, col] = float(any(st[col] != "intact" for st in relevant))

    mat_filled = np.where(np.isnan(mat), 0.0, mat)
    if len(uids) > 1:
        dist  = pdist(mat_filled, metric="jaccard")
        dist  = np.nan_to_num(dist, nan=1.0)
        Z     = linkage(dist, method="average")
        order = leaves_list(Z)
    else:
        order = [0]
    mat  = mat[order]
    uids = [uids[i] for i in order]

    labels = [canonicals[u]["gene"] for u in uids]
    n_prot = len(uids)
    n_iso  = sum(len(uid_to_states[u]) for u in uids)

    fig_h = max(5, n_prot * 0.28)
    fig, ax = plt.subplots(figsize=(7, fig_h))
    masked = np.ma.masked_invalid(mat)
    cmap = matplotlib.colors.ListedColormap(["#f7f7f7", "#C44E52"])
    cmap.set_bad(color="#d0d0d0")
    im = ax.imshow(masked, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")

    ax.set_xticks(range(N_POS))
    ax.set_xticklabels([f"$\\beta\\alpha${k}" for k in range(1, N_POS + 1)], fontsize=14)
    ax.set_yticks(range(n_prot))
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel("Barrel motif position", fontsize=15)
    ax.set_title(
        f"Per-protein motif disruption across AS isoforms (strict filter)\n"
        f"{n_prot} canonical proteins, {n_iso} isoforms  "
        f"(red = disrupted in ≥1 isoform; white = intact in all isoforms)",
        fontsize=16,
    )
    ax.set_xticks(np.arange(-0.5, N_POS, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_prot, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_ticks([0.25, 0.75])
    cbar.set_ticklabels(["intact in all", "disrupted in ≥1"], fontsize=12)
    cbar.ax.tick_params(length=0)

    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


_DISRUPTION_COLORS = ["#f7f7f7", "#FFB3C1", "#C44E52", "#7B2FBE"]
_DISRUPTION_BOUNDS = [-0.5, 0.5, 1.5, 2.5, 3.5]
_DISRUPTION_LABELS = ["0", "1", "2", "3"]


def _gene_disruption_matrix(isoforms, combined_map, canonicals):
    """Build per-gene disruption count matrix (capped at 3) and cluster rows by Jaccard."""
    uid_to_states = {}
    for iso in isoforms:
        uid, iso_id = iso["uniprot_id"], iso["isoform_id"]
        if iso_id in combined_map:
            uid_to_states.setdefault(uid, []).append(combined_map[iso_id])

    N_POS = 8
    uids  = sorted(uid_to_states)
    mat   = np.full((len(uids), N_POS), np.nan)
    for row, uid in enumerate(uids):
        kp = len(canonicals[uid]["motifs"])
        for col in range(min(kp, N_POS)):
            relevant = [st for st in uid_to_states[uid] if len(st) > col]
            if relevant:
                count = sum(1 for st in relevant if st[col] != "intact")
                mat[row, col] = float(min(count, 3))  # cap at 3+

    mat_filled = np.where(np.isnan(mat), 0.0, mat)
    if len(uids) > 1:
        mat_bin = (mat_filled > 0).astype(float)
        dist    = pdist(mat_bin, metric="jaccard")
        Z       = linkage(np.nan_to_num(dist, nan=1.0), method="average")
        order   = leaves_list(Z)
    else:
        order = [0]

    mat    = mat[order]
    uids   = [uids[i] for i in order]
    labels = [canonicals[u]["gene"] for u in uids]
    n_iso  = sum(len(v) for v in uid_to_states.values())
    return mat, labels, n_iso


def _draw_disruption_panel(ax, mat, labels, N_POS, cmap, norm):
    masked = np.ma.masked_invalid(mat)
    im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(N_POS))
    ax.set_xticklabels([f"$\\beta\\alpha${k}" for k in range(1, N_POS + 1)], fontsize=14)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel("Barrel motif position", fontsize=15)
    ax.set_xticks(np.arange(-0.5, N_POS, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", length=0)
    return im


def plot_gene_disruption_heatmap(isoforms, combined_map, canonicals, out):
    mat, labels, n_iso = _gene_disruption_matrix(isoforms, combined_map, canonicals)
    n_prot = len(labels)
    N_POS  = 8

    cmap = matplotlib.colors.ListedColormap(_DISRUPTION_COLORS)
    norm = matplotlib.colors.BoundaryNorm(_DISRUPTION_BOUNDS, cmap.N)
    cmap.set_bad(color="#d0d0d0")

    fig_h = max(5, n_prot * 0.28)
    fig, ax = plt.subplots(figsize=(7, fig_h))
    im = _draw_disruption_panel(ax, mat, labels, N_POS, cmap, norm)
    ax.set_title(
        f"Per-gene motif disruption by isoform count  ({n_prot} genes, {n_iso} isoforms)\n"
        "colour = number of isoforms with disruption at that position",
        fontsize=16,
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.06, ticks=[0, 1, 2, 3])
    cbar.set_ticklabels(_DISRUPTION_LABELS, fontsize=13)
    cbar.ax.tick_params(length=0)
    cbar.ax.set_title("Affected\nisoforms", fontsize=12, pad=4)
    fig.tight_layout()
    renderer = fig.canvas.get_renderer()
    tb = cbar.ax.get_tightbbox(renderer)
    if tb is not None:
        cx = (tb.x0 + tb.x1) / 2 / (fig.get_figwidth() * fig.dpi)
        cy = tb.y0 / (fig.get_figheight() * fig.dpi)
        fig.text(cx, cy - 0.015, "■ gray = absent\n(< 8 motifs)",
                 ha="center", va="top", fontsize=12, color="#666666")
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def plot_gene_disruption_heatmap_split(isoforms, combined_map, canonicals, out):
    """Same data as plot_gene_disruption_heatmap but split into two side-by-side panels."""
    mat, labels, n_iso = _gene_disruption_matrix(isoforms, combined_map, canonicals)
    n_prot = len(labels)
    N_POS  = 8

    cmap = matplotlib.colors.ListedColormap(_DISRUPTION_COLORS)
    norm = matplotlib.colors.BoundaryNorm(_DISRUPTION_BOUNDS, cmap.N)
    cmap.set_bad(color="#d0d0d0")

    half         = (n_prot + 1) // 2
    mat_l, lbl_l = mat[:half],  labels[:half]
    mat_r, lbl_r = mat[half:],  labels[half:]

    fig_h = max(5, half * 0.28)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, fig_h),
                                    gridspec_kw={"wspace": 0.35})
    im = _draw_disruption_panel(ax1, mat_l, lbl_l, N_POS, cmap, norm)
    _draw_disruption_panel(ax2, mat_r, lbl_r, N_POS, cmap, norm)

    ax1.set_title("1/2", fontsize=13, pad=4)
    ax2.set_title("2/2", fontsize=13, pad=4)

    cbar = fig.colorbar(im, ax=[ax1, ax2], fraction=0.015, pad=0.08, ticks=[0, 1, 2, 3])
    cbar.set_ticklabels(_DISRUPTION_LABELS, fontsize=13)
    cbar.ax.tick_params(length=0)
    cbar.ax.set_title("Affected\nisoforms", fontsize=12, pad=4)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    # Centre title over the two panels (not the full figure including colorbar)
    renderer = fig.canvas.get_renderer()
    b1 = ax1.get_position()
    b2 = ax2.get_position()
    panels_cx = (b1.x0 + b2.x1) / 2
    fig.text(panels_cx, 0.995,
             f"Per-gene motif disruption by isoform count\n{n_prot} genes, {n_iso} isoforms",
             ha="center", va="top", fontsize=16)

    tb = cbar.ax.get_tightbbox(renderer)
    if tb is not None:
        cx = (tb.x0 + tb.x1) / 2 / (fig.get_figwidth() * fig.dpi)
        cy = tb.y0 / (fig.get_figheight() * fig.dpi)
        fig.text(cx, cy - 0.015, "■ gray = absent\n(< 8 motifs)",
                 ha="center", va="top", fontsize=12, color="#666666")
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Analysis 4 — VSP boundary enrichment
# ---------------------------------------------------------------------------

def compute_length_null(canonicals, enst_matched):
    counts = defaultdict(int)
    total  = 0
    for uid, can in canonicals.items():
        if uid not in enst_matched:
            continue
        for r in range(can["ds"], can["de"]):
            counts[_tau5(r, can["motifs"])] += 1
            total += 1
    return {t: counts[t] / total if total > 0 else 0.0 for t in CATS5}


def load_vsp_boundaries_strict(isoforms_strict, canonicals):
    """VSP boundaries for strictly-filtered isoforms, with divergence check."""
    boundaries = []
    for iso in isoforms_strict:
        uid     = iso["uniprot_id"]
        can     = canonicals[uid]
        can_seq = can["seq"]
        iso_seq = iso["seq"]
        ds, de  = can["ds"], can["de"]

        diverge = None
        for i, (ca, ia) in enumerate(zip(can_seq, iso_seq), start=1):
            if ca != ia:
                diverge = i
                break

        for vsp in iso["vsps"]:
            vs = vsp.get("can_start")
            ve = vsp.get("can_end")
            if vs is None or ve is None:
                continue
            if ve < ds or vs >= de:
                continue
            if diverge is None or diverge >= ve:
                continue
            boundaries.append(dict(uid=uid, can_start=max(vs, ds), can_end=min(ve, de - 1)))
    return boundaries


def plot_boundary_enrichment(rho_t, N, pi_null, pvals_bh, pvals_raw, title, out):
    fig, ax = plt.subplots(figsize=(8, 4))
    x        = np.arange(len(CATS5))
    rho_vals = [rho_t[t] for t in CATS5]
    errs     = bh_ci_halfwidths(pvals_raw, pi_null, N)
    err_lo   = [min(rho_vals[i], errs[t]) for i, t in enumerate(CATS5)]
    err_hi   = [errs[t] for t in CATS5]

    ax.bar(x, rho_vals, width=0.65, color=[COLS5[t] for t in CATS5], alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[err_lo, err_hi],
                fmt="none", color="black", capsize=5, lw=1.2, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, t in enumerate(CATS5):
        sl = sig_stars(pvals_bh[t])
        if sl != "ns":
            ax.text(x[i], rho_vals[i] + err_hi[i] + 0.05, sl,
                    ha="center", va="bottom", fontsize=16, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[t] for t in CATS5], fontsize=14)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=15)
    ax.set_title(title, fontsize=16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def plot_combined_boundary_enrichment(rho_s, rho_e, N, pi_null,
                                      pvals_s_bh, pvals_e_bh,
                                      pvals_s_raw, pvals_e_raw, out):
    w  = 0.35
    x  = np.arange(len(CATS5))
    xs = x - w / 2
    xe = x + w / 2

    bh_s = bh_ci_halfwidths(pvals_s_raw, pi_null, N)
    bh_e = bh_ci_halfwidths(pvals_e_raw, pi_null, N)
    rho_s_vals = [rho_s[t] for t in CATS5]
    rho_e_vals = [rho_e[t] for t in CATS5]
    err_s_lo = [min(rho_s_vals[i], bh_s[t]) for i, t in enumerate(CATS5)]
    err_s_hi = [bh_s[t] for t in CATS5]
    err_e_lo = [min(rho_e_vals[i], bh_e[t]) for i, t in enumerate(CATS5)]
    err_e_hi = [bh_e[t] for t in CATS5]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(xs, rho_s_vals, w, color=[COLS5[t] for t in CATS5], alpha=0.90, zorder=3, label="VSP start")
    ax.bar(xe, rho_e_vals, w, color=[COLS5[t] for t in CATS5], alpha=0.45, zorder=3, hatch="//", label="VSP end")
    ax.errorbar(xs, rho_s_vals, yerr=[err_s_lo, err_s_hi], fmt="none", color="black", capsize=4, lw=1.0, zorder=4)
    ax.errorbar(xe, rho_e_vals, yerr=[err_e_lo, err_e_hi], fmt="none", color="black", capsize=4, lw=1.0, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, t in enumerate(CATS5):
        sl = sig_stars(pvals_s_bh[t]); el = sig_stars(pvals_e_bh[t])
        if sl != "ns": ax.text(xs[i], rho_s_vals[i] + err_s_hi[i] + 0.04, sl, ha="center", va="bottom", fontsize=14, fontweight="bold")
        if el != "ns": ax.text(xe[i], rho_e_vals[i] + err_e_hi[i] + 0.04, el, ha="center", va="bottom", fontsize=14, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[t] for t in CATS5], fontsize=14)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=15)
    ax.set_title(f"VSP boundary enrichment (strict filter)\n$N$ = {N} VSPs  |  error bars = BH-adjusted CI", fontsize=16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=14, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def plot_pooled_enrichment(N_t_s, N_t_e, N_total, pi_null, pvals_bh, pvals_raw, out):
    N_t    = {t: N_t_s[t] + N_t_e[t] for t in CATS5}
    f_t    = {t: N_t[t] / N_total for t in CATS5}
    rho_t  = {t: f_t[t] / pi_null[t] if pi_null[t] > 0 else 0.0 for t in CATS5}

    x        = np.arange(len(CATS5))
    rho_vals = [rho_t[t] for t in CATS5]
    hw       = bh_ci_halfwidths(pvals_raw, pi_null, N_total)
    err_lo   = [min(rho_vals[i], hw[t]) for i, t in enumerate(CATS5)]
    err_hi   = [hw[t] for t in CATS5]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(x, rho_vals, width=0.6, color=[COLS5[t] for t in CATS5], alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[err_lo, err_hi], fmt="none", color="black", capsize=5, lw=1.2, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)
    for i, t in enumerate(CATS5):
        sl = sig_stars(pvals_bh[t])
        if sl != "ns":
            ax.text(x[i], rho_vals[i] + err_hi[i] + 0.04, sl, ha="center", va="bottom", fontsize=14, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[t] for t in CATS5], fontsize=14)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=15)
    ax.set_title(f"VSP boundary enrichment (strict filter, start + end pooled)\n$N$ = {N_total}  |  error bars = BH-adjusted CI", fontsize=16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")
    return N_t, f_t, rho_t


# ---------------------------------------------------------------------------
# Analysis 3 — Transcript-derived AS boundary enrichment
# ---------------------------------------------------------------------------

def find_resync(can_seq, can_end, iso_seq):
    offsets = [0] + [s for d in range(1, MAX_SLIDE + 1) for s in (d, -d)]
    for delta in offsets:
        pos = can_end + delta
        if pos < 0 or pos >= len(can_seq):
            continue
        window = min(MIN_MATCH, len(can_seq) - pos)
        if window < MIN_MATCH // 2:
            continue
        target  = can_seq[pos : pos + window]
        iso_idx = iso_seq.find(target)
        if iso_idx >= 0:
            return pos, iso_idx, delta
    return None


def find_as_boundaries_strict(isoforms_strict, canonicals):
    records = []
    for iso in isoforms_strict:
        can_id  = iso["uniprot_id"]
        can     = canonicals[can_id]
        can_seq = can["seq"]
        iso_seq = iso["seq"]
        ds, de  = can["ds"], can["de"]

        diverge = None
        for i, (ca, ia) in enumerate(zip(can_seq, iso_seq), start=1):
            if ca != ia:
                diverge = i
                break
        if diverge is None:
            continue

        for vsp in sorted(iso["vsps"], key=lambda v: v.get("can_start", 0)):
            can_end = vsp.get("can_end")
            if can_end is None:
                continue
            resync = find_resync(can_seq, can_end, iso_seq)
            if resync is not None:
                R_can = resync[0]
                end   = R_can - 1
            else:
                end = can_end

            if diverge >= end:
                continue

            s_t = max(diverge, ds)
            e_t = min(end, de - 1)
            if s_t > e_t:
                continue
            records.append((can_id, s_t, e_t))
    return records


def compute_enrichment(positions, canonicals, pi_null):
    N_t = {c: 0 for c in CATS5}
    N   = 0
    for uid, pos in positions:
        if uid not in canonicals:
            continue
        can = canonicals[uid]
        if not (can["ds"] <= pos < can["de"]):
            continue
        N_t[_tau5(pos, can["motifs"])] += 1
        N += 1
    return N_t, N


def plot_splice_junctions(res_s_f, pv_raw_s, pv_bh_s, res_e_f, pv_raw_e, pv_bh_e,
                           pi_null, N_s, N_e, out, title=None):
    w  = 0.35
    x  = np.arange(len(CATS5))
    xs = x - w / 2
    xe = x + w / 2

    rho_s = [res_s_f[t] / pi_null[t] if pi_null[t] > 0 else 0.0 for t in CATS5]
    rho_e = [res_e_f[t] / pi_null[t] if pi_null[t] > 0 else 0.0 for t in CATS5]

    hw_s    = bh_ci_halfwidths(pv_raw_s, pi_null, N_s)
    hw_e    = bh_ci_halfwidths(pv_raw_e, pi_null, N_e)
    err_s_lo = [min(rho_s[i], hw_s[t]) for i, t in enumerate(CATS5)]
    err_s_hi = [hw_s[t] for t in CATS5]
    err_e_lo = [min(rho_e[i], hw_e[t]) for i, t in enumerate(CATS5)]
    err_e_hi = [hw_e[t] for t in CATS5]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(xs, rho_s, w, color=[COLS5[c] for c in CATS5], alpha=0.90, zorder=3, label=f"AS start ($D_{{\\rm seq}}$, $N$={N_s})")
    ax.bar(xe, rho_e, w, color=[COLS5[c] for c in CATS5], alpha=0.45, zorder=3, hatch="//", label=f"AS end ($R_{{\\rm can}}-1$, $N$={N_e})")
    ax.errorbar(xs, rho_s, yerr=[err_s_lo, err_s_hi], fmt="none", color="black", capsize=4, lw=1.0, zorder=4)
    ax.errorbar(xe, rho_e, yerr=[err_e_lo, err_e_hi], fmt="none", color="black", capsize=4, lw=1.0, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, c in enumerate(CATS5):
        sl = sig_stars(pv_bh_s[c]); el = sig_stars(pv_bh_e[c])
        if sl != "ns": ax.text(xs[i], rho_s[i] + err_s_hi[i] + 0.04, sl, ha="center", va="bottom", fontsize=14, fontweight="bold")
        if el != "ns": ax.text(xe[i], rho_e[i] + err_e_hi[i] + 0.04, el, ha="center", va="bottom", fontsize=14, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[c] for c in CATS5], fontsize=14)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=15)
    if title is None:
        title = f"Transcript-derived AS boundary enrichment (strict filter, all positions)\nerror bars = BH-adjusted CI"
    ax.set_title(title, fontsize=16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=14, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def plot_pooled_splice(N_t_s, N_t_e, N_total, pi_null, pvals_bh, pvals_raw, out, title=None):
    N_t   = {c: N_t_s[c] + N_t_e[c] for c in CATS5}
    f_t   = {c: N_t[c] / N_total for c in CATS5}
    rho_t = {c: f_t[c] / pi_null[c] if pi_null[c] > 0 else 0.0 for c in CATS5}

    x        = np.arange(len(CATS5))
    rho_vals = [rho_t[c] for c in CATS5]
    hw       = bh_ci_halfwidths(pvals_raw, pi_null, N_total)
    err_lo   = [min(rho_vals[i], hw[c]) for i, c in enumerate(CATS5)]
    err_hi   = [hw[c] for c in CATS5]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(x, rho_vals, width=0.6, color=[COLS5[c] for c in CATS5], alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[err_lo, err_hi], fmt="none", color="black", capsize=5, lw=1.2, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)
    for i, c in enumerate(CATS5):
        sl = sig_stars(pvals_bh[c])
        if sl != "ns":
            ax.text(x[i], rho_vals[i] + err_hi[i] + 0.04, sl, ha="center", va="bottom", fontsize=14, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[c] for c in CATS5], fontsize=14)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=15)
    if title is None:
        title = f"Transcript-derived AS boundary enrichment (strict filter, pooled)\n$N$ = {N_total}  |  error bars = BH-adjusted CI"
    ax.set_title(title, fontsize=16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")
    return N_t, f_t, rho_t


# ---------------------------------------------------------------------------
# Motif-core null (4 categories, simple length-weighted)
# ---------------------------------------------------------------------------

def compute_motif_core_null(canonicals, enst_matched):
    """Length-weighted null from [first_beta_start, last_alpha_end] only."""
    counts = defaultdict(int)
    total  = 0
    for uid, can in canonicals.items():
        if uid not in enst_matched:
            continue
        motifs = can["motifs"]
        if not motifs:
            continue
        core_s = motifs[0]["beta_start"]
        core_e = motifs[-1]["alpha_end"]
        for r in range(core_s, core_e + 1):
            t = _tau5(r, motifs)
            if t in CATS4:
                counts[t] += 1
                total += 1
    return {t: counts[t] / total if total > 0 else 0.0 for t in CATS4}


# ---------------------------------------------------------------------------
# Analysis 1 — canonical junction enrichment (motif-core null)
# ---------------------------------------------------------------------------

def run_analysis1_core(proteins):
    """
    Junction-count-weighted enrichment restricted to [first_beta_start, last_alpha_end].
    Only junctions within the motif core contribute. 4 categories (no flanking).
    """
    N_t = {c: 0   for c in CATS4}
    wq  = {c: 0.0 for c in CATS4}
    N   = 0

    for p in proteins:
        motifs = p["motifs"]
        if not motifs:
            continue
        core_s = motifs[0]["beta_start"]
        core_e = motifs[-1]["alpha_end"]
        ep_len_core = core_e - core_s + 1
        if ep_len_core <= 0:
            continue

        core_jcts = [j for j in p["junctions"] if core_s <= j <= core_e]
        n_core = len(core_jcts)
        if n_core == 0:
            continue

        # Junction-count-weighted null from core positions
        type_counts = defaultdict(int)
        for pos in range(core_s, core_e + 1):
            t = _tau5(pos, motifs)
            if t in CATS4:
                type_counts[t] += 1
        for t in CATS4:
            wq[t] += n_core * type_counts[t] / ep_len_core

        for j in core_jcts:
            t = _tau5(j, motifs)
            if t in CATS4:
                N_t[t] += 1
                N += 1

    if N == 0:
        return None
    pi_t  = {c: wq[c] / N for c in CATS4}
    f_t   = {c: N_t[c] / N for c in CATS4}
    rho_t = {c: f_t[c] / pi_t[c] if pi_t[c] > 0 else 0.0 for c in CATS4}
    return dict(N=N, N_t=N_t, pi_t=pi_t, f_t=f_t, rho_t=rho_t)


def plot_bar_4cat(rho_t, N, pi_null, pvals_bh, pvals_raw, title, out, N_t=None):
    fig, ax = plt.subplots(figsize=(7, 4))
    x        = np.arange(len(CATS4))
    rho_vals = [rho_t[t] for t in CATS4]
    errs     = bh_ci_halfwidths_4cat(pvals_raw, pi_null, N)
    err_lo   = [min(rho_vals[i], errs[t]) for i, t in enumerate(CATS4)]
    err_hi   = [errs[t] for t in CATS4]

    ax.bar(x, rho_vals, width=0.6, color=[COLS5[t] for t in CATS4], alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[err_lo, err_hi],
                fmt="none", color="black", capsize=5, lw=1.2, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)
    for i, t in enumerate(CATS4):
        sl = sig_stars(pvals_bh[t])
        y = max(rho_vals[i] + err_hi[i] + 0.04, 1.06)
        if N_t is not None:
            ax.text(x[i], y, f"{N * pi_null[t]:.1f}/{N_t[t]}",
                    ha="center", va="bottom", fontsize=12.5, color="#444444")
            y += 0.19
        if sl != "ns":
            ax.text(x[i], y, sl, ha="center", va="bottom", fontsize=14, fontweight="bold")
    # Compute y ceiling to prevent text clipping
    _tops = []
    for i, t in enumerate(CATS4):
        y = max(rho_vals[i] + err_hi[i] + 0.04, 1.06)
        if N_t is not None: y += 0.19
        if sig_stars(pvals_bh[t]) != "ns": y += 0.18
        _tops.append(y)
    ax.set_ylim(bottom=0, top=max(_tops) + 0.12)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS4[t] for t in CATS4], fontsize=15)
    ax.set_ylabel("enrichment ratio", fontsize=15)
    ax.set_title(title, fontsize=16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Analysis 2 — motif-specific junction enrichment (motif-core null)
# ---------------------------------------------------------------------------

def _tau_motif_idx(pos, motifs):
    """Return integer category index for Analysis 2."""
    if not motifs or pos < motifs[0]["beta_start"]:
        return _FLANKING
    for i, m in enumerate(motifs):
        k = m["motif"]
        if m["beta_start"] <= pos <= m["beta_end"]:
            return _CAT_IDX.get(("beta", k), _FLANKING)
        if m["beta_end"] < pos < m["alpha_start"]:
            return _CAT_IDX.get(("loop", k), _FLANKING)
        if m["alpha_start"] <= pos <= m["alpha_end"]:
            return _CAT_IDX.get(("alpha", k), _FLANKING)
        if i + 1 < len(motifs) and m["alpha_end"] < pos < motifs[i + 1]["beta_start"]:
            return _CAT_IDX.get(("inter", k), _FLANKING)
    return _FLANKING


def run_analysis2_core(proteins_kp8):
    """
    Motif-specific enrichment for K_p=8 proteins restricted to motif core.
    31 primary categories: (beta/loop/alpha, k=1..8) + (inter, k=1..7).
    Junction-count-weighted null uses core length as denominator.
    """
    N_t = defaultdict(int)
    wq  = defaultdict(float)
    N   = 0

    for p in proteins_kp8:
        motifs = p["motifs"]
        if len(motifs) != 8:
            continue
        core_s = motifs[0]["beta_start"]
        core_e = motifs[-1]["alpha_end"]
        ep_len_core = core_e - core_s + 1
        if ep_len_core <= 0:
            continue

        core_jcts = [j for j in p["junctions"] if core_s <= j <= core_e]
        n_core = len(core_jcts)
        if n_core == 0:
            continue

        # Null: core positions only
        for pos in range(core_s, core_e + 1):
            cat_i = _tau_motif_idx(pos, motifs)
            key   = _IDX_CAT.get(cat_i)
            if key is not None:
                wq[key] += n_core / ep_len_core

        for j in core_jcts:
            cat_i = _tau_motif_idx(j, motifs)
            key   = _IDX_CAT.get(cat_i)
            if key is not None:
                N_t[key] += 1
            N += 1

    if N == 0:
        return None
    pi_t  = {cat: wq[cat] / N for cat in PRIMARY_CATS}
    f_t   = {cat: N_t[cat] / N for cat in PRIMARY_CATS}
    rho_t = {cat: f_t[cat] / pi_t[cat] if pi_t[cat] > 0 else float("nan")
             for cat in PRIMARY_CATS}
    return dict(N=N, N_t=N_t, pi_t=pi_t, f_t=f_t, rho_t=rho_t)


def plot_analysis2_heatmap(result, pvals_bh, out):
    """4×8 heatmap of rho values for beta/loop/alpha/inter × motif 1..8."""
    ELEM_ORDER  = ["beta", "loop", "alpha", "inter"]
    ELEM_LABELS = {"beta": "β-strand", "loop": "β→α loop",
                   "alpha": "α-helix", "inter": "α→β loop"}

    data = np.full((4, 8), np.nan)
    sig  = np.full((4, 8), "ns", dtype=object)
    for ri, elem in enumerate(ELEM_ORDER):
        n_k = 7 if elem == "inter" else 8   # inter only exists between motifs 1–7
        for k in range(1, n_k + 1):
            cat = (elem, k)
            rho = result["rho_t"].get(cat, np.nan)
            data[ri, k - 1] = rho
            sig[ri, k - 1]  = sig_stars(pvals_bh.get(cat, 1.0))

    valid = data[~np.isnan(data)]
    vmax = max(np.abs(valid - 1.0).max(), 0.3) + 1.0
    vmin = 2.0 - vmax

    cmap = matplotlib.cm.get_cmap("RdBu_r").copy()
    cmap.set_bad(color="#d0d0d0")
    masked = np.ma.masked_invalid(data)

    fig, ax = plt.subplots(figsize=(10, 4.2))
    im = ax.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")

    ax.set_xticks(range(8))
    ax.set_xticklabels([f"βα{k}" for k in range(1, 9)], fontsize=15)
    ax.set_yticks(range(4))
    ax.set_yticklabels([ELEM_LABELS[e] for e in ELEM_ORDER], fontsize=15)

    for ri in range(4):
        for ci in range(8):
            s = sig[ri, ci]
            if s != "ns":
                ax.text(ci, ri, s, ha="center", va="center", fontsize=14, fontweight="bold")

    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02,
                 label="enrichment ratio")
    ax.set_title(
        f"Motif-specific junction enrichment   {result['N']} junctions\n"
        "* BH-adjusted p < 0.05",
        fontsize=16)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Core-restricted enrichment for Analysis 3
# ---------------------------------------------------------------------------

def compute_enrichment_core(positions, canonicals):
    """Like compute_enrichment but restricts to motif core; uses CATS4."""
    N_t = {c: 0 for c in CATS4}
    N   = 0
    for uid, pos in positions:
        if uid not in canonicals:
            continue
        motifs = canonicals[uid]["motifs"]
        if not motifs:
            continue
        core_s = motifs[0]["beta_start"]
        core_e = motifs[-1]["alpha_end"]
        if not (core_s <= pos <= core_e):
            continue
        t = _tau5(pos, motifs)
        if t in CATS4:
            N_t[t] += 1
            N += 1
    return N_t, N


def plot_boundary_enrichment_4cat(rho_s, rho_e, N_s, N_e, pi_null,
                                   pv_s_bh, pv_e_bh, pv_s_raw, pv_e_raw, title, out,
                                   ylabel="observed / motif-core null",
                                   N_t_s=None, N_t_e=None):
    w  = 0.35
    x  = np.arange(len(CATS4))
    xs = x - w / 2
    xe = x + w / 2

    bh_s = bh_ci_halfwidths_4cat(pv_s_raw, pi_null, N_s)
    bh_e = bh_ci_halfwidths_4cat(pv_e_raw, pi_null, N_e)
    rho_s_vals = [rho_s[t] for t in CATS4]
    rho_e_vals = [rho_e[t] for t in CATS4]
    err_s_lo = [min(rho_s_vals[i], bh_s[t]) for i, t in enumerate(CATS4)]
    err_s_hi = [bh_s[t] for t in CATS4]
    err_e_lo = [min(rho_e_vals[i], bh_e[t]) for i, t in enumerate(CATS4)]
    err_e_hi = [bh_e[t] for t in CATS4]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(xs, rho_s_vals, w, color=[COLS5[t] for t in CATS4], alpha=0.90, zorder=3, label="start")
    ax.bar(xe, rho_e_vals, w, color=[COLS5[t] for t in CATS4], alpha=0.45, zorder=3,
           hatch="//", label="end")
    ax.errorbar(xs, rho_s_vals, yerr=[err_s_lo, err_s_hi], fmt="none", color="black", capsize=4, lw=1.0, zorder=4)
    ax.errorbar(xe, rho_e_vals, yerr=[err_e_lo, err_e_hi], fmt="none", color="black", capsize=4, lw=1.0, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)
    for i, t in enumerate(CATS4):
        sl = sig_stars(pv_s_bh[t]); el = sig_stars(pv_e_bh[t])
        ys = rho_s_vals[i] + err_s_hi[i] + 0.04
        ye = rho_e_vals[i] + err_e_hi[i] + 0.04
        if N_t_s is not None:
            pass  # count labels omitted from paired chart; see Results table
        if N_t_e is not None:
            pass
        if sl != "ns": ax.text(xs[i], ys, sl,
                                ha="center", va="bottom", fontsize=14, fontweight="bold")
        if el != "ns": ax.text(xe[i], ye, el,
                                ha="center", va="bottom", fontsize=14, fontweight="bold")
    # Compute y ceiling to prevent text clipping
    _tops = []
    for i, t in enumerate(CATS4):
        for rv, ev, pv in [
            (rho_s_vals[i], err_s_hi[i], pv_s_bh[t]),
            (rho_e_vals[i], err_e_hi[i], pv_e_bh[t]),
        ]:
            y = rv + ev + 0.04
            if sig_stars(pv) != "ns": y += 0.16
            _tops.append(y)
    ax.set_ylim(bottom=0, top=max(_tops) + 0.12)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS4[t] for t in CATS4], fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_title(title, fontsize=16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=14, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def plot_pooled_4cat(N_t_s, N_t_e, N_total, pi_null, pvals_bh, pvals_raw, title, out,
                     ylabel="enrichment ratio", show_counts=False):
    N_t    = {t: N_t_s[t] + N_t_e[t] for t in CATS4}
    f_t    = {t: N_t[t] / N_total for t in CATS4}
    rho_t  = {t: f_t[t] / pi_null[t] if pi_null[t] > 0 else 0.0 for t in CATS4}

    x        = np.arange(len(CATS4))
    rho_vals = [rho_t[t] for t in CATS4]
    hw       = bh_ci_halfwidths_4cat(pvals_raw, pi_null, N_total)
    err_lo   = [min(rho_vals[i], hw[t]) for i, t in enumerate(CATS4)]
    err_hi   = [hw[t] for t in CATS4]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(x, rho_vals, width=0.6, color=[COLS5[t] for t in CATS4], alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[err_lo, err_hi], fmt="none", color="black", capsize=5, lw=1.2, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)
    for i, t in enumerate(CATS4):
        sl = sig_stars(pvals_bh[t])
        y = rho_vals[i] + err_hi[i] + 0.14
        if show_counts:
            ax.text(x[i], y, f"{N_total * pi_null[t]:.1f}/{N_t[t]}",
                    ha="center", va="bottom", fontsize=12.5, color="#444444")
            y += 0.19
        if sl != "ns":
            ax.text(x[i], y, sl, ha="center", va="bottom", fontsize=14, fontweight="bold")
    # Compute y ceiling to prevent text clipping
    _tops = []
    for i, t in enumerate(CATS4):
        y = rho_vals[i] + err_hi[i] + 0.14
        if show_counts: y += 0.19
        if sig_stars(pvals_bh[t]) != "ns": y += 0.18
        _tops.append(y)
    ax.set_ylim(bottom=0, top=max(_tops) + 0.12)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS4[t] for t in CATS4], fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_title(title, fontsize=16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")
    return N_t, f_t, rho_t


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    db_path    = get_config().db_path
    conn       = sqlite3.connect(db_path)
    _raw_junctions = load_canonical_junctions(conn)   # keyed (uid, domain_index)
    # Flatten to uid -> sorted union of junctions across all domain instances
    canonical_junctions: dict[str, list[int]] = {}
    for (uid, _didx), jcts in _raw_junctions.items():
        canonical_junctions.setdefault(uid, [])
        canonical_junctions[uid] = sorted(set(canonical_junctions[uid]) | set(jcts))
    enst_matched = set(canonical_junctions.keys())
    canonicals = load_canonical(conn)
    all_isos   = load_all_isoforms(conn)
    conn.close()

    isoforms = apply_strict_filter(all_isos, canonicals)

    # Ensembl-matched canonicals (for null + analyses 3/4)
    enst_canonicals = {uid: v for uid, v in canonicals.items() if uid in enst_matched}
    # Ensembl-matched strict isoforms (for analyses 3/4)
    enst_isoforms   = [iso for iso in isoforms if iso["uniprot_id"] in enst_canonicals]

    # Motif-core null (4 categories, no flanking) — used in Analyses 1–4
    pi_null_core = compute_motif_core_null(canonicals, enst_matched)

    SEPARATOR = "\n" + "="*72

    # -----------------------------------------------------------------------
    # Dataset summary
    # -----------------------------------------------------------------------
    print(SEPARATOR)
    print("  DATASET SUMMARY (strict motif-core filter)")
    print(SEPARATOR)

    n_total_iso  = len(isoforms)
    uids_with    = {iso["uniprot_id"] for iso in isoforms}
    n_prot_with  = len(uids_with)
    n_prot_all   = len(canonicals)

    per_prot = Counter(iso["uniprot_id"] for iso in isoforms)
    iso_per_prot_dist = Counter(per_prot.values())

    print(f"Total canonical proteins:             {n_prot_all}")
    print(f"Canonical proteins with >=1 isoform:  {n_prot_with}")
    print(f"Canonical proteins with no isoform:   {n_prot_all - n_prot_with}")
    print(f"Total AS isoforms (strict):            {n_total_iso}")
    print()
    print("Isoforms per canonical protein:")
    for k in sorted(iso_per_prot_dist):
        print(f"  {k}: {iso_per_prot_dist[k]} proteins")

    print("\nMotif-core null (length-weighted, 4 categories):")
    for t in CATS4:
        print(f"  {LABELS4[t]:<16}: {pi_null_core[t]:.3f}")

    # -----------------------------------------------------------------------
    # Analysis 1 — canonical splice-junction enrichment (motif-core null)
    # -----------------------------------------------------------------------
    print(SEPARATOR)
    print("  ANALYSIS 1 — Canonical splice-junction enrichment (motif-core null)")
    print(SEPARATOR)

    proteins_a1 = []
    for uid in sorted(enst_matched):
        if uid not in canonicals:
            continue
        can  = canonicals[uid]
        if not can["motifs"]:
            continue
        jcts = canonical_junctions.get(uid, [])
        if not jcts:
            continue
        proteins_a1.append({"motifs": can["motifs"], "junctions": jcts})

    print(f"\nProteins with Ensembl junctions and motifs: {len(proteins_a1)}")

    r1 = run_analysis1_core(proteins_a1)
    if r1 is not None:
        N1, N1_t, pi1, f1, rho1 = r1["N"], r1["N_t"], r1["pi_t"], r1["f_t"], r1["rho_t"]
        pv1_raw = {}
        for t in CATS4:
            E_t = N1 * pi1[t]
            z = (N1_t[t] - E_t) / np.sqrt(E_t) if E_t > 0 else 0.0
            pv1_raw[t] = float(2 * scipy_stats.norm.sf(abs(z)))
        pv1_bh = bh_correct(pv1_raw)
        obs1 = np.array([N1_t[t] for t in CATS4], dtype=float)
        exp1 = np.array([N1 * pi1[t] for t in CATS4], dtype=float)
        chi2_1, p_gl_1 = scipy_stats.chisquare(obs1, f_exp=exp1)
        print_table_4cat("Canonical junctions (motif-core)", N1_t, N1, pi1,
                         f1, rho1, pv1_raw, pv1_bh, chi2_1, p_gl_1)
        plot_bar_4cat(rho1, N1, pi1, pv1_bh, pv1_raw,
                      f"Canonical splice junction enrichment\n"
                      f"{N1} junctions, error bars = BH-adjusted CI\n* BH-adjusted p < 0.05",
                      OUTDIR / "junction_enrichment_core.png",
                      N_t=N1_t)

    # -----------------------------------------------------------------------
    # Analysis 2 — motif-specific junction enrichment, K_p=8 (motif-core null)
    # -----------------------------------------------------------------------
    print(SEPARATOR)
    print("  ANALYSIS 2 — Motif-specific junction enrichment, K_p=8 (motif-core null)")
    print(SEPARATOR)

    proteins_a2 = [p for p in proteins_a1 if len(p["motifs"]) == 8]
    print(f"\nK_p=8 proteins with Ensembl junctions: {len(proteins_a2)}")

    r2 = run_analysis2_core(proteins_a2)
    if r2 is not None:
        N2 = r2["N"]
        pv2_raw = {}
        for cat in PRIMARY_CATS:
            pi_cat = r2["pi_t"].get(cat, 0)
            n_cat  = r2["N_t"].get(cat, 0)
            E_cat  = N2 * pi_cat
            if E_cat > 0:
                z = (n_cat - E_cat) / np.sqrt(E_cat)
                pv2_raw[cat] = float(2 * scipy_stats.norm.sf(abs(z)))
            else:
                pv2_raw[cat] = 1.0
        pv2_bh = bh_correct(pv2_raw)

        obs2 = np.array([r2["N_t"].get(cat, 0) for cat in PRIMARY_CATS], dtype=float)
        exp2 = np.array([N2 * r2["pi_t"].get(cat, 0) for cat in PRIMARY_CATS], dtype=float)
        chi2_2, p_gl_2 = scipy_stats.chisquare(obs2, f_exp=exp2)

        print(f"\nTotal core junctions used: {N2}")
        print(f"Global chi2({len(PRIMARY_CATS)-1}) = {chi2_2:.2f},  p = {p_gl_2:.4g}")

        sig_cats = [(cat, r2["rho_t"][cat], pv2_bh[cat])
                    for cat in PRIMARY_CATS if pv2_bh.get(cat, 1.0) < 0.05]
        print(f"\nSignificant categories (BH p < 0.05): {len(sig_cats)}")
        if sig_cats:
            print(f"  {'Category':<22}  {'rho':>6}  {'BH p':>8}  Sig")
            print("  " + "-"*46)
            for cat, rho, pq in sorted(sig_cats, key=lambda x: x[2]):
                elem, k = cat
                label = f"{LABELS4.get(elem, elem)} k={k}"
                print(f"  {label:<22}  {rho:>6.3f}  {pq:>8.4f}  {sig_stars(pq)}")

        print("\nAll categories:")
        print(f"  {'Category':<22}  {'N_t':>4}  {'pi_t':>7}  {'f_t':>7}  {'rho':>6}  {'BH p':>8}  Sig")
        print("  " + "-"*68)
        for cat in PRIMARY_CATS:
            elem, k = cat
            label = f"{elem} k={k}"
            print(f"  {label:<22}  {r2['N_t'].get(cat,0):>4}  "
                  f"{r2['pi_t'].get(cat,0):>7.4f}  {r2['f_t'].get(cat,0):>7.4f}  "
                  f"{r2['rho_t'].get(cat,float('nan')):>6.3f}  "
                  f"{pv2_bh.get(cat,1.0):>8.4f}  {sig_stars(pv2_bh.get(cat,1.0))}")

        plot_analysis2_heatmap(r2, pv2_bh, OUTDIR / "motif_enrichment_heatmap.png")

    # -----------------------------------------------------------------------
    # Analysis 3 — transcript-derived boundaries (motif-core null, CATS4)
    # -----------------------------------------------------------------------
    print(SEPARATOR)
    print("  ANALYSIS 3 — Transcript-derived AS boundaries (motif-core null)")
    print(SEPARATOR)

    boundaries3 = find_as_boundaries_strict(enst_isoforms, enst_canonicals)
    starts3 = [(uid, s) for uid, s, e in boundaries3]
    ends3   = [(uid, e) for uid, s, e in boundaries3]
    print(f"\nAS boundary pairs (strict, Ensembl-matched): {len(boundaries3)}")

    if boundaries3:
        N_t_s3, N_s3 = compute_enrichment_core(starts3, enst_canonicals)
        N_t_e3, N_e3 = compute_enrichment_core(ends3,   enst_canonicals)
        print(f"  Start positions within core: {N_s3},  end positions within core: {N_e3}")

        if N_s3 > 0:
            f_t_s3, rho_t_s3, pv_s3_raw, pv_s3_bh, _ = chi_square_enrichment_4cat(
                N_t_s3, N_s3, pi_null_core)
            obs_s3 = np.array([N_t_s3[t] for t in CATS4], dtype=float)
            exp_s3 = np.array([N_s3 * pi_null_core[t] for t in CATS4], dtype=float)
            chi2_s3, p_gl_s3 = scipy_stats.chisquare(obs_s3, f_exp=exp_s3)
            print_table_4cat("Transcript start positions (D_seq)", N_t_s3, N_s3, pi_null_core,
                             f_t_s3, rho_t_s3, pv_s3_raw, pv_s3_bh, chi2_s3, p_gl_s3)

        if N_e3 > 0:
            f_t_e3, rho_t_e3, pv_e3_raw, pv_e3_bh, _ = chi_square_enrichment_4cat(
                N_t_e3, N_e3, pi_null_core)
            obs_e3 = np.array([N_t_e3[t] for t in CATS4], dtype=float)
            exp_e3 = np.array([N_e3 * pi_null_core[t] for t in CATS4], dtype=float)
            chi2_e3, p_gl_e3 = scipy_stats.chisquare(obs_e3, f_exp=exp_e3)
            print_table_4cat("Transcript end positions (R_can - 1)", N_t_e3, N_e3, pi_null_core,
                             f_t_e3, rho_t_e3, pv_e3_raw, pv_e3_bh, chi2_e3, p_gl_e3)

        if N_s3 > 0 and N_e3 > 0:
            N_pool3 = N_s3 + N_e3
            pv_p3_raw = {}
            for t in CATS4:
                N_t_p = N_t_s3[t] + N_t_e3[t]
                E_t = N_pool3 * pi_null_core[t]
                z   = (N_t_p - E_t) / np.sqrt(E_t) if pi_null_core[t] > 0 else 0.0
                pv_p3_raw[t] = float(2 * scipy_stats.norm.sf(abs(z)))
            pv_p3_bh = bh_correct(pv_p3_raw)
            N_t_p3 = {t: N_t_s3[t] + N_t_e3[t] for t in CATS4}
            f_p3   = {t: N_t_p3[t] / N_pool3 for t in CATS4}
            rho_p3 = {t: f_p3[t] / pi_null_core[t] if pi_null_core[t] > 0 else 0.0 for t in CATS4}
            print("\nPooled (start + end):")
            for t in CATS4:
                print(f"  {LABELS4[t]:<16}: N={N_t_p3[t]}, f={f_p3[t]:.3f}, "
                      f"rho={rho_p3[t]:.3f}, BH p={pv_p3_bh[t]:.4f} {sig_stars(pv_p3_bh[t])}")

            plot_boundary_enrichment_4cat(
                rho_t_s3, rho_t_e3, N_s3, N_e3, pi_null_core,
                pv_s3_bh, pv_e3_bh, pv_s3_raw, pv_e3_raw,
                f"AS boundary enrichment\n"
                f"N_start = {N_s3}, N_end = {N_e3}  |  error bars = BH-adjusted CI\n* BH-adjusted p < 0.05",
                OUTDIR / "as_splice_junctions.png",
                N_t_s=N_t_s3, N_t_e=N_t_e3)
            plot_pooled_4cat(
                N_t_s3, N_t_e3, N_pool3, pi_null_core, pv_p3_bh, pv_p3_raw,
                f"AS junctions enrichment\n"
                f"{N_pool3} junctions, error bars = BH-adjusted CI\n* BH-adjusted p < 0.05",
                OUTDIR / "as_splice_junctions_pooled.png",
                show_counts=True)

        # -------------------------------------------------------------------
        # Full-domain view (CATS5, all 113 boundary positions)
        # Includes starts outside core (enters/spans) and ends outside core
        # -------------------------------------------------------------------
        print("\n-- Full-domain view (all boundary positions, including those outside core) --")

        # Classify pairs: contained / enters / exits / spans
        n_a3_both = 0; n_a3_enters = 0; n_a3_exits = 0; n_a3_spans = 0
        for uid, s, e in boundaries3:
            can = enst_canonicals.get(uid)
            if not can or not can["motifs"]:
                continue
            cs = can["motifs"][0]["beta_start"]
            ce = can["motifs"][-1]["alpha_end"]
            s_in = (cs <= s <= ce); e_in = (cs <= e <= ce)
            if   s_in and e_in:       n_a3_both   += 1
            elif not s_in and e_in:   n_a3_enters += 1
            elif s_in and not e_in:   n_a3_exits  += 1
            else:                     n_a3_spans  += 1

        print(f"\nBoundary-pair classification (relative to motif core):")
        print(f"  Contained  (both within core):           {n_a3_both}")
        print(f"  Enters     (start outside, end inside):  {n_a3_enters}")
        print(f"  Exits      (start inside, end outside):  {n_a3_exits}")
        print(f"  Spans      (both outside core):          {n_a3_spans}")

        pi_null_full = compute_length_null(canonicals, enst_matched)
        N_t_s3f, N_s3f = compute_enrichment(starts3, enst_canonicals, pi_null_full)
        N_t_e3f, N_e3f = compute_enrichment(ends3,   enst_canonicals, pi_null_full)

        f_t_s3f, rho_t_s3f, pv_s3f_raw, pv_s3f_bh, _ = chi_square_enrichment(
            N_t_s3f, N_s3f, pi_null_full)
        obs_s3f = np.array([N_t_s3f[t] for t in CATS5], dtype=float)
        exp_s3f = np.array([N_s3f * pi_null_full[t] for t in CATS5], dtype=float)
        chi2_s3f, p_gl_s3f = scipy_stats.chisquare(obs_s3f, f_exp=exp_s3f)
        print_table("Transcript start positions (all, D_seq)", N_t_s3f, N_s3f, pi_null_full,
                    f_t_s3f, rho_t_s3f, pv_s3f_raw, pv_s3f_bh, chi2_s3f, p_gl_s3f)

        f_t_e3f, rho_t_e3f, pv_e3f_raw, pv_e3f_bh, _ = chi_square_enrichment(
            N_t_e3f, N_e3f, pi_null_full)
        obs_e3f = np.array([N_t_e3f[t] for t in CATS5], dtype=float)
        exp_e3f = np.array([N_e3f * pi_null_full[t] for t in CATS5], dtype=float)
        chi2_e3f, p_gl_e3f = scipy_stats.chisquare(obs_e3f, f_exp=exp_e3f)
        print_table("Transcript end positions (all, R_can-1)", N_t_e3f, N_e3f, pi_null_full,
                    f_t_e3f, rho_t_e3f, pv_e3f_raw, pv_e3f_bh, chi2_e3f, p_gl_e3f)

        N_pool3f = N_s3f + N_e3f
        pv_p3f_raw = {}
        for t in CATS5:
            N_t_p = N_t_s3f[t] + N_t_e3f[t]
            E_t = N_pool3f * pi_null_full[t]
            z   = (N_t_p - E_t) / np.sqrt(E_t) if pi_null_full[t] > 0 else 0.0
            pv_p3f_raw[t] = float(2 * scipy_stats.norm.sf(abs(z)))
        pv_p3f_bh = bh_correct(pv_p3f_raw)
        N_t_p3f = {t: N_t_s3f[t] + N_t_e3f[t] for t in CATS5}
        f_p3f   = {t: N_t_p3f[t] / N_pool3f for t in CATS5}
        rho_p3f = {t: f_p3f[t] / pi_null_full[t] if pi_null_full[t] > 0 else 0.0 for t in CATS5}
        print("\nPooled (start + end, all positions):")
        for t in CATS5:
            print(f"  {LABELS5[t]:<16}: N={N_t_p3f[t]}, f={f_p3f[t]:.3f}, "
                  f"rho={rho_p3f[t]:.3f}, BH p={pv_p3f_bh[t]:.4f} {sig_stars(pv_p3f_bh[t])}")

        plot_splice_junctions(
            f_t_s3f, pv_s3f_raw, pv_s3f_bh,
            f_t_e3f, pv_e3f_raw, pv_e3f_bh,
            pi_null_full, N_s3f, N_e3f,
            OUTDIR / "as_splice_junctions_full.png",
            title=(f"Transcript-derived AS boundary enrichment — all positions\n"
                   f"strict isoforms, domain null  |  error bars = BH-adjusted CI\n* BH-adjusted p < 0.05"))
        plot_pooled_splice(
            N_t_s3f, N_t_e3f, N_pool3f, pi_null_full, pv_p3f_bh, pv_p3f_raw,
            OUTDIR / "as_splice_junctions_full_pooled.png",
            title=(f"Transcript-derived AS boundary enrichment — pooled, all positions\n"
                   f"$N$ = {N_pool3f}  |  error bars = BH-adjusted CI\n* BH-adjusted p < 0.05"))

        # -------------------------------------------------------------------
        # Part C — AS boundaries vs. canonical junction null (from Analysis 1)
        # Tests whether AS boundaries are distributed differently from where
        # canonical splice junctions land, controlling for existing splice-site
        # preferences (e.g. β-strands are already slightly preferred by
        # canonical junctions).
        # -------------------------------------------------------------------
        if r1 is not None:
            pi_null_jct = r1["f_t"]
            print("\n-- Part C: AS boundaries vs. canonical junction null (from Analysis 1) --")
            print("  Null = empirical frequency of canonical splice junctions across motif-core elements")
            print("  rho > 1: more common at AS boundaries than at canonical splice sites")
            print("\nCanonical junction null (from Analysis 1):")
            for t in CATS4:
                print(f"  {LABELS4[t]:<16}: {pi_null_jct[t]:.4f}")

            if N_s3 > 0:
                f_t_s3j, rho_t_s3j, pv_s3j_raw, pv_s3j_bh, _ = chi_square_enrichment_4cat(
                    N_t_s3, N_s3, pi_null_jct)
                obs_s3j = np.array([N_t_s3[t] for t in CATS4], dtype=float)
                exp_s3j = np.array([N_s3 * pi_null_jct[t] for t in CATS4], dtype=float)
                chi2_s3j, p_gl_s3j = scipy_stats.chisquare(obs_s3j, f_exp=exp_s3j)
                print_table_4cat("Start positions (D_seq) vs. junction null", N_t_s3, N_s3,
                                 pi_null_jct, f_t_s3j, rho_t_s3j,
                                 pv_s3j_raw, pv_s3j_bh, chi2_s3j, p_gl_s3j)

            if N_e3 > 0:
                f_t_e3j, rho_t_e3j, pv_e3j_raw, pv_e3j_bh, _ = chi_square_enrichment_4cat(
                    N_t_e3, N_e3, pi_null_jct)
                obs_e3j = np.array([N_t_e3[t] for t in CATS4], dtype=float)
                exp_e3j = np.array([N_e3 * pi_null_jct[t] for t in CATS4], dtype=float)
                chi2_e3j, p_gl_e3j = scipy_stats.chisquare(obs_e3j, f_exp=exp_e3j)
                print_table_4cat("End positions (R_can-1) vs. junction null", N_t_e3, N_e3,
                                 pi_null_jct, f_t_e3j, rho_t_e3j,
                                 pv_e3j_raw, pv_e3j_bh, chi2_e3j, p_gl_e3j)

            if N_s3 > 0 and N_e3 > 0:
                N_pool3j = N_s3 + N_e3
                pv_p3j_raw = {}
                for t in CATS4:
                    N_t_p = N_t_s3[t] + N_t_e3[t]
                    E_t   = N_pool3j * pi_null_jct[t]
                    z     = (N_t_p - E_t) / np.sqrt(E_t) if pi_null_jct[t] > 0 else 0.0
                    pv_p3j_raw[t] = float(2 * scipy_stats.norm.sf(abs(z)))
                pv_p3j_bh = bh_correct(pv_p3j_raw)
                N_t_p3j = {t: N_t_s3[t] + N_t_e3[t] for t in CATS4}
                f_p3j   = {t: N_t_p3j[t] / N_pool3j for t in CATS4}
                rho_p3j = {t: f_p3j[t] / pi_null_jct[t] if pi_null_jct[t] > 0 else 0.0
                           for t in CATS4}
                print("\nPooled (start + end) vs. junction null:")
                for t in CATS4:
                    print(f"  {LABELS4[t]:<16}: N={N_t_p3j[t]}, f={f_p3j[t]:.3f}, "
                          f"rho={rho_p3j[t]:.3f}, BH p={pv_p3j_bh[t]:.4f} "
                          f"{sig_stars(pv_p3j_bh[t])}")

                _jct_ylabel = "AS boundaries / canonical junction null"
                plot_boundary_enrichment_4cat(
                    rho_t_s3j, rho_t_e3j, N_s3, N_e3, pi_null_jct,
                    pv_s3j_bh, pv_e3j_bh, pv_s3j_raw, pv_e3j_raw,
                    f"AS boundary enrichment vs. canonical junction null\n"
                    f"N_start = {N_s3}, N_end = {N_e3}  |  error bars = BH-adjusted CI\n* BH-adjusted p < 0.05",
                    OUTDIR / "as_junctions_vs_canonical_null.png",
                    ylabel=_jct_ylabel,
                    N_t_s=N_t_s3, N_t_e=N_t_e3)
                plot_pooled_4cat(
                    N_t_s3, N_t_e3, N_pool3j, pi_null_jct,
                    pv_p3j_bh, pv_p3j_raw,
                    f"AS boundary enrichment vs. canonical junction null (pooled)\n"
                    f"N = {N_pool3j}  |  error bars = BH-adjusted CI\n* BH-adjusted p < 0.05",
                    OUTDIR / "as_junctions_vs_canonical_null_pooled.png",
                    ylabel=_jct_ylabel,
                    show_counts=True)

    # -----------------------------------------------------------------------
    # Analysis 4 — VSP boundary enrichment (motif-core null, CATS4)
    # -----------------------------------------------------------------------
    print(SEPARATOR)
    print("  ANALYSIS 4 — VSP boundary enrichment (motif-core null)")
    print(SEPARATOR)

    boundaries = load_vsp_boundaries_strict(enst_isoforms, canonicals)
    n_vsps     = len(boundaries)
    n_covered  = len(set(b["uid"] for b in boundaries))
    print(f"\nVSP spans (Ensembl-matched, strict, diverge-checked): {n_vsps}  ({n_covered} proteins)")

    if n_vsps > 0:
        # Start — restrict to motif core
        counts_s4 = {c: 0 for c in CATS4}
        N_s4 = 0
        for b in boundaries:
            motifs = canonicals[b["uid"]]["motifs"]
            if not motifs:
                continue
            core_s = motifs[0]["beta_start"]
            core_e = motifs[-1]["alpha_end"]
            pos = b["can_start"]
            if core_s <= pos <= core_e:
                t = _tau5(pos, motifs)
                if t in CATS4:
                    counts_s4[t] += 1
                    N_s4 += 1

        # End — restrict to motif core
        counts_e4 = {c: 0 for c in CATS4}
        N_e4 = 0
        for b in boundaries:
            motifs = canonicals[b["uid"]]["motifs"]
            if not motifs:
                continue
            core_s = motifs[0]["beta_start"]
            core_e = motifs[-1]["alpha_end"]
            pos = b["can_end"]
            if core_s <= pos <= core_e:
                t = _tau5(pos, motifs)
                if t in CATS4:
                    counts_e4[t] += 1
                    N_e4 += 1

        print(f"  VSP starts within motif core: {N_s4},  ends within core: {N_e4}")

        if N_s4 > 0:
            f_t_s4, rho_t_s4, pv_s4_raw, pv_s4_bh, _ = chi_square_enrichment_4cat(
                counts_s4, N_s4, pi_null_core)
            obs_s4 = np.array([counts_s4[t] for t in CATS4], dtype=float)
            exp_s4 = np.array([N_s4 * pi_null_core[t] for t in CATS4], dtype=float)
            chi2_s4, p_gl_s4 = scipy_stats.chisquare(obs_s4, f_exp=exp_s4)
            print_table_4cat("VSP start positions (core)", counts_s4, N_s4, pi_null_core,
                             f_t_s4, rho_t_s4, pv_s4_raw, pv_s4_bh, chi2_s4, p_gl_s4)

        if N_e4 > 0:
            f_t_e4, rho_t_e4, pv_e4_raw, pv_e4_bh, _ = chi_square_enrichment_4cat(
                counts_e4, N_e4, pi_null_core)
            obs_e4 = np.array([counts_e4[t] for t in CATS4], dtype=float)
            exp_e4 = np.array([N_e4 * pi_null_core[t] for t in CATS4], dtype=float)
            chi2_e4, p_gl_e4 = scipy_stats.chisquare(obs_e4, f_exp=exp_e4)
            print_table_4cat("VSP end positions (core)", counts_e4, N_e4, pi_null_core,
                             f_t_e4, rho_t_e4, pv_e4_raw, pv_e4_bh, chi2_e4, p_gl_e4)

        if N_s4 > 0 and N_e4 > 0:
            N_pooled4 = N_s4 + N_e4
            counts_pool4 = {t: counts_s4[t] + counts_e4[t] for t in CATS4}
            pv_pool4_raw = {}
            for t in CATS4:
                E_t = N_pooled4 * pi_null_core[t]
                z   = (counts_pool4[t] - E_t) / np.sqrt(E_t) if pi_null_core[t] > 0 else 0.0
                pv_pool4_raw[t] = float(2 * scipy_stats.norm.sf(abs(z)))
            pv_pool4_bh = bh_correct(pv_pool4_raw)
            f_pool4  = {t: counts_pool4[t] / N_pooled4 for t in CATS4}
            rho_pool4 = {t: f_pool4[t] / pi_null_core[t] if pi_null_core[t] > 0 else 0.0
                         for t in CATS4}
            print("\nPooled (start + end):")
            for t in CATS4:
                print(f"  {LABELS4[t]:<16}: N={counts_pool4[t]}, f={f_pool4[t]:.3f}, "
                      f"rho={rho_pool4[t]:.3f}, BH p={pv_pool4_bh[t]:.4f} {sig_stars(pv_pool4_bh[t])}")

            plot_boundary_enrichment_4cat(
                rho_t_s4, rho_t_e4, N_s4, N_e4, pi_null_core,
                pv_s4_bh, pv_e4_bh, pv_s4_raw, pv_e4_raw,
                f"VSP boundary enrichment within motif core (strict filter)\n"
                f"N_start = {N_s4}, N_end = {N_e4}  |  error bars = BH-adjusted CI\n* BH-adjusted p < 0.05",
                OUTDIR / "vsp_boundary_enrichment.png")
            plot_pooled_4cat(
                counts_s4, counts_e4, N_pooled4, pi_null_core, pv_pool4_bh, pv_pool4_raw,
                f"VSP boundary enrichment\n"
                f"N = {N_pooled4}  |  error bars = BH-adjusted CI\n* BH-adjusted p < 0.05",
                OUTDIR / "vsp_boundary_pooled.png")

    # -----------------------------------------------------------------------
    # Analysis 5 — structural impact on barrel architecture
    # -----------------------------------------------------------------------
    print(SEPARATOR)
    print("  ANALYSIS 5 — Structural impact on barrel architecture (strict filter)")
    print(SEPARATOR)

    combined_map, combined_list = {}, []
    separate_map, separate_list = {}, []
    detailed_map = {}
    for iso in isoforms:
        cm = classify_combined(iso, canonicals)
        sm = classify_separate(iso, canonicals)
        combined_map[iso["isoform_id"]] = cm
        separate_map[iso["isoform_id"]] = sm
        combined_list.append(cm)
        separate_list.append(sm)
        detailed_map[iso["isoform_id"]] = classify_combined_detailed(iso, canonicals)

    intact_combined = [sum(1 for s in st if s == "intact") for st in combined_list]
    print(f"\n[Combined motif] Mean intact: {np.mean(intact_combined):.2f}  "
          f"Median: {np.median(intact_combined):.1f}")
    print("Intact-motif distribution:")
    for k, v in sorted(Counter(intact_combined).items()):
        print(f"  {k}: {v} isoforms ({100*v/len(intact_combined):.1f}%)")

    max_pos = max(len(st) for st in combined_list)
    print("\nPer-position disruption rate (combined):")
    for i in range(max_pos):
        relevant  = [st for st in combined_list if len(st) > i]
        disrupted = sum(1 for st in relevant if st[i] != "intact")
        rate = 100 * disrupted / len(relevant)
        print(f"  Position {i+1}: {disrupted}/{len(relevant)} ({rate:.1f}%)")

    n_fully_intact = sum(1 for st in combined_list if all(s == "intact" for s in st))
    print(f"\nFully intact barrel isoforms: {n_fully_intact} / {len(combined_list)} "
          f"({100*n_fully_intact/len(combined_list):.1f}%)")

    plot_disruption_combined(combined_list, OUTDIR / "as_domain_disruption.png")
    plot_isoform_heatmap(isoforms, combined_map, canonicals,
                         OUTDIR / "isoform_disruption_heatmap.png")
    plot_isoform_heatmap_by_gene(isoforms, combined_map, canonicals,
                                 OUTDIR / "isoform_disruption_heatmap_by_gene.png")
    # plot_isoform_heatmap_detailed(isoforms, detailed_map, canonicals,
    #                               OUTDIR / "isoform_disruption_heatmap_detailed.png")
    # plot_motif_heatmap(isoforms, combined_map, canonicals,
    #                    OUTDIR / "motif_disruption_heatmap.png")
    plot_gene_disruption_heatmap(isoforms, combined_map, canonicals,
                                 OUTDIR / "gene_disruption_count_heatmap.png")
    plot_gene_disruption_heatmap_split(isoforms, combined_map, canonicals,
                                       OUTDIR / "gene_disruption_count_heatmap_split.png")

    print(SEPARATOR)
    print("  Done.")
    print(SEPARATOR)


if __name__ == "__main__":
    main()
