#!/usr/bin/env python3
"""
VSP span boundary placement analysis (section 10 of Statistical-Analysis.md).

  10A. Element distribution of VSP START positions vs length-weighted null
  10B. Element distribution of VSP END positions vs length-weighted null
  10C. Start-end asymmetry (paired permutation)

Start position: domain-clipped to max(can_start, ds).
End position:   domain-clipped to min(can_end, de - 1).

Figures:
  figures/vsp_start_enrichment.png
  figures/vsp_end_enrichment.png
  figures/vsp_boundary_asymmetry.png
"""

import argparse
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config
from junction_utils import load_canonical_junctions


CATS5    = ["beta", "alpha", "inter", "loop", "flanking"]
LABELS5  = {"beta": "beta-strand", "alpha": "alpha-helix", "inter": "Inter-motif",
             "loop": "Loop (b->a)", "flanking": "Flanking"}
MDLABELS = {"beta": "β-strand", "alpha": "α-helix", "inter": "Inter-motif",
             "loop": "Loop (β→α)", "flanking": "Flanking"}
COLS5    = {"beta": "#4C72B0", "alpha": "#DD8452", "inter": "#55A868",
            "loop": "#8c8c8c", "flanking": "#C44E52"}


# ---------------------------------------------------------------------------
# Element assignment
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
# Data loading
# ---------------------------------------------------------------------------

def load_proteins(conn):
    enst_matched = set(load_canonical_junctions(conn).keys())
    rows = conn.execute("""
        SELECT vc.uniprot_id, vc.gene_name, vc.domain_start, vc.domain_end,
               vc.exon_annotations, vc.motif_annotations, ca.sequence
        FROM   view_canonical vc
        JOIN   canonical_analysis ca ON ca.uniprot_id = vc.uniprot_id
    """).fetchall()
    proteins = {}
    for uid, gene, ds, de, ea, ma, seq in rows:
        if uid not in enst_matched:
            continue
        motifs = json.loads(ma)
        exons     = json.loads(ea)
        junctions = [e["end"] for e in exons[:-1] if ds <= e["end"] < de]
        proteins[uid] = dict(uid=uid, gene=gene, ds=ds, de=de,
                             seq=seq, motifs=motifs, junctions=junctions)
    return proteins


def load_isoform_sequences(conn):
    rows = conn.execute(
        "SELECT isoform_id, sequence FROM isoforms WHERE sequence IS NOT NULL"
    ).fetchall()
    return {r[0]: r[1] for r in rows}


def load_vsp_boundaries(conn, proteins, iso_seqs):
    """
    For each VSP with domain overlap return domain-clipped start and end positions.
    can_start -> max(vs, ds);  can_end -> min(ve, de-1).

    Events are excluded when the isoform sequence does not actually diverge from
    the canonical before can_end (annotation artefacts / identical sequences).
    """
    rows = conn.execute("""
        SELECT uniprot_id, isoform_id, vsp_domain_events
        FROM   view_noncanonical
    """).fetchall()

    boundaries = []
    for uid, iso_id, vsp_json in rows:
        if uid not in proteins:
            continue
        if iso_id not in iso_seqs:
            continue
        can_seq = proteins[uid]["seq"]
        iso_seq = iso_seqs[iso_id]

        # Find first position where sequences diverge (1-indexed)
        diverge = None
        for i, (ca, ia) in enumerate(zip(can_seq, iso_seq), start=1):
            if ca != ia:
                diverge = i
                break

        vsps = json.loads(vsp_json)
        p = proteins[uid]
        ds, de = p["ds"], p["de"]
        for vsp in vsps:
            vs = vsp.get("can_start")
            ve = vsp.get("can_end")
            if vs is None or ve is None:
                continue
            if ve < ds or vs >= de:          # no domain overlap
                continue
            # Exclude events where isoform doesn't diverge within the VSP region
            if diverge is None or diverge >= ve:
                continue
            boundaries.append(dict(
                uid       = uid,
                can_start = max(vs, ds),
                can_end   = min(ve, de - 1),
            ))
    return boundaries


# ---------------------------------------------------------------------------
# Length-weighted null
# ---------------------------------------------------------------------------

def compute_length_null(proteins):
    counts = defaultdict(int)
    total  = 0
    for p in proteins.values():
        for r in range(p["ds"], p["de"]):
            counts[_tau5(r, p["motifs"])] += 1
            total += 1
    return {t: counts[t] / total if total > 0 else 0.0 for t in CATS5}


# ---------------------------------------------------------------------------
# VSP residue coverage (Ochoa-Leyva-style: which elements are inside spans)
# ---------------------------------------------------------------------------

def compute_vsp_residue_coverage(boundaries, proteins, pi_null):
    """
    For each VSP span [can_start, can_end], count every residue position and
    classify it by structural element type.  Multiple VSPs can cover the same
    residue; each is counted separately (per-VSP coverage, not per-protein).
    """
    counts = defaultdict(int)
    total  = 0
    for bnd in boundaries:
        uid = bnd["uid"]
        if uid not in proteins:
            continue
        motifs = proteins[uid]["motifs"]
        for pos in range(bnd["can_start"], bnd["can_end"] + 1):
            counts[_tau5(pos, motifs)] += 1
            total += 1
    if total == 0:
        return None
    f_t   = {t: counts[t] / total for t in CATS5}
    rho_t = {t: f_t[t] / pi_null[t] if pi_null[t] > 0 else 0.0 for t in CATS5}
    return dict(N=total, counts=counts, f_t=f_t, rho_t=rho_t)


# ---------------------------------------------------------------------------
# 10A / 10B: enrichment analysis
# ---------------------------------------------------------------------------

def compute_observed(boundaries, proteins, pi_null, endpoint):
    key = "can_start" if endpoint == "start" else "can_end"
    counts = defaultdict(int)
    for bnd in boundaries:
        counts[_tau5(bnd[key], proteins[bnd["uid"]]["motifs"])] += 1
    N     = len(boundaries)
    f_t   = {t: counts[t] / N for t in CATS5}
    rho_t = {t: f_t[t] / pi_null[t] if pi_null[t] > 0 else float("nan") for t in CATS5}
    return counts, f_t, rho_t


def chi_square_pvalues(boundaries, proteins, pi_null, endpoint):
    """
    Compute chi-square p-values for 10A or 10B enrichment.

    Returns (N_t, f_t, rho_t, pvals_raw, pvals_bh, z_scores) — all dicts
    keyed by element type (CATS5).
    """
    N_t, f_t, rho_t = compute_observed(boundaries, proteins, pi_null, endpoint)
    N = len(boundaries)

    z_scores  = {}
    pvals_raw = {}
    for t in CATS5:
        pi = pi_null[t]
        if pi <= 0:
            z_scores[t]  = float("nan")
            pvals_raw[t] = float("nan")
            continue
        E_t = N * pi
        z   = (N_t[t] - E_t) / np.sqrt(E_t)
        z_scores[t]  = float(z)
        pvals_raw[t] = float(2 * norm.sf(abs(z)))

    pvals_bh = bh_correct(pvals_raw, CATS5)
    return N_t, f_t, rho_t, pvals_raw, pvals_bh, z_scores


def bh_correct(pvals_dict, keys):
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


def _bh_ci_halfwidths(pvals_raw, pi_null, N, alpha=0.05):
    """
    CI half-widths in rho scale consistent with BH-corrected stars.

    Uses the test-consistent formula: half_width_t = z_t / sqrt(N * pi_t),
    where z_t = norm.ppf(1 - alpha_eff_t / 2) and
    alpha_eff_t = alpha * rank_t / m  (BH threshold for element t's rank).

    This ensures: CI excludes 1.0  <=>  BH p < alpha, so stars and error
    bars tell the same story.
    """
    m      = len(CATS5)
    order  = sorted(CATS5, key=lambda t: pvals_raw[t])
    ranks  = {t: r for r, t in enumerate(order, 1)}
    errs   = {}
    for t in CATS5:
        alpha_eff = alpha * ranks[t] / m
        z_t = norm.ppf(1 - alpha_eff / 2)
        E_t = N * pi_null[t]
        errs[t] = z_t / np.sqrt(E_t) if E_t > 0 else 0.0
    return errs


def plot_enrichment(rho_t, N, pi_null, endpoint_label, pvals_bh, pvals_raw, out):
    """
    Bar chart of enrichment ratios with BH-adjusted CI error bars centred at rho_t.
    Half-width = z_BH(t) / sqrt(E_t), where z_BH(t) uses the BH rank-adjusted
    critical value so that CI excludes 1.0 exactly when BH p < 0.05.
    """
    fig, ax = plt.subplots(figsize=(8, 4))
    x        = np.arange(len(CATS5))
    rho_vals = [rho_t[t] for t in CATS5]

    errs     = _bh_ci_halfwidths(pvals_raw, pi_null, N)
    err_lo   = [min(rho_vals[i], errs[t]) for i, t in enumerate(CATS5)]
    err_hi   = [errs[t] for t in CATS5]

    ax.bar(x, rho_vals, width=0.65, color=[COLS5[t] for t in CATS5], alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[err_lo, err_hi],
                fmt="none", color="black", capsize=5, lw=1.2, zorder=4,
                label="BH-adjusted CI")
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, t in enumerate(CATS5):
        p   = pvals_bh[t]
        sig = "**" if p < 0.01 else ("*" if p < 0.05 else "")
        if sig:
            ax.text(x[i], rho_vals[i] + err_hi[i] + 0.05,
                    sig, ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([MDLABELS[t] for t in CATS5], fontsize=9)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=10)
    ax.set_title(
        f"VSP {endpoint_label} position enrichment vs length-weighted null\n"
        f"$N$ = {N}  |  error bars = BH-adjusted CI;  ** BH p < 0.01,  * BH p < 0.05",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.18))
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


def plot_combined_enrichment(rho_s, rho_e, N, pi_null,
                             pvals_s_bh, pvals_e_bh,
                             pvals_s_raw, pvals_e_raw, out):
    """
    Grouped bar chart showing start and end enrichment ratios together,
    one pair of bars per structural element.
    """
    w   = 0.35
    x   = np.arange(len(CATS5))
    xs  = x - w / 2
    xe  = x + w / 2

    bh_s = _bh_ci_halfwidths(pvals_s_raw, pi_null, N)
    bh_e = _bh_ci_halfwidths(pvals_e_raw, pi_null, N)
    rho_s_vals = [rho_s[t] for t in CATS5]
    rho_e_vals = [rho_e[t] for t in CATS5]
    err_s_lo = [min(rho_s_vals[i], bh_s[t]) for i, t in enumerate(CATS5)]
    err_s_hi = [bh_s[t] for t in CATS5]
    err_e_lo = [min(rho_e_vals[i], bh_e[t]) for i, t in enumerate(CATS5)]
    err_e_hi = [bh_e[t] for t in CATS5]

    fig, ax = plt.subplots(figsize=(9, 4.5))

    ax.bar(xs, rho_s_vals, w, color=[COLS5[t] for t in CATS5],
           alpha=0.90, zorder=3, label="VSP start")
    ax.bar(xe, rho_e_vals, w, color=[COLS5[t] for t in CATS5],
           alpha=0.45, zorder=3, hatch="//", label="VSP end")

    ax.errorbar(xs, rho_s_vals, yerr=[err_s_lo, err_s_hi],
                fmt="none", color="black", capsize=4, lw=1.0, zorder=4)
    ax.errorbar(xe, rho_e_vals, yerr=[err_e_lo, err_e_hi],
                fmt="none", color="black", capsize=4, lw=1.0, zorder=4)

    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    def sig_label(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return ""

    for i, t in enumerate(CATS5):
        top_s = rho_s_vals[i] + err_s_hi[i]
        top_e = rho_e_vals[i] + err_e_hi[i]
        sl = sig_label(pvals_s_bh[t])
        el = sig_label(pvals_e_bh[t])
        if sl:
            ax.text(xs[i], top_s + 0.04, sl, ha="center", va="bottom",
                    fontsize=9, fontweight="bold")
        if el:
            ax.text(xe[i], top_e + 0.04, el, ha="center", va="bottom",
                    fontsize=9, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([MDLABELS[t] for t in CATS5], fontsize=9)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=10)
    ax.set_title(
        f"VSP boundary enrichment in TIM-barrel structural elements\n"
        f"$N$ = {N} VSPs  |  error bars = BH-adjusted CI;  *, **, *** = BH-adjusted $p$",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.14), ncol=2)

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


def plot_pooled_enrichment(N_t_s, N_t_e, N_total, pi_null, pvals_bh, pvals_raw, out):
    """
    Single bar chart treating start and end positions as one pooled set.
    N_total = 2 * number of VSPs; N_t = N_t_start + N_t_end per element.
    """
    N_t    = {t: N_t_s[t] + N_t_e[t] for t in CATS5}
    f_t    = {t: N_t[t] / N_total for t in CATS5}
    rho_t  = {t: f_t[t] / pi_null[t] if pi_null[t] > 0 else 0.0 for t in CATS5}

    x        = np.arange(len(CATS5))
    rho_vals = [rho_t[t] for t in CATS5]
    bh_hw    = _bh_ci_halfwidths(pvals_raw, pi_null, N_total)
    err_lo   = [min(rho_vals[i], bh_hw[t]) for i, t in enumerate(CATS5)]
    err_hi   = [bh_hw[t] for t in CATS5]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(x, rho_vals, width=0.6, color=[COLS5[t] for t in CATS5],
           alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[err_lo, err_hi],
                fmt="none", color="black", capsize=5, lw=1.2, zorder=4,
                label="BH-adjusted CI")
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    def sig_label(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return ""

    for i, t in enumerate(CATS5):
        sl = sig_label(pvals_bh[t])
        if sl:
            ax.text(x[i], rho_vals[i] + err_hi[i] + 0.04, sl,
                    ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([MDLABELS[t] for t in CATS5], fontsize=9)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=10)
    ax.set_title(
        "VSP boundary enrichment in TIM-barrel structural elements\n"
        f"(start + end pooled, $N$ = {N_total}  |  error bars = BH-adjusted CI)",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="upper left")

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    return N_t, f_t, rho_t


# ---------------------------------------------------------------------------
# VSP residue coverage figure
# ---------------------------------------------------------------------------

def plot_residue_coverage(cov, pi_null, pvals_bh, out):
    from scipy.stats import chi2 as _chi2
    x        = np.arange(len(CATS5))
    rho_vals = [cov["rho_t"][t] for t in CATS5]
    N        = cov["N"]

    lo_errs, hi_errs = [], []
    for t in CATS5:
        k   = cov["counts"][t]
        E_t = N * pi_null[t]
        if E_t > 0 and k > 0:
            lo = _chi2.ppf(0.025, 2 * k)       / (2 * E_t)
            hi = _chi2.ppf(0.975, 2 * (k + 1)) / (2 * E_t)
            lo_errs.append(max(cov["rho_t"][t] - lo, 0.0))
            hi_errs.append(max(hi - cov["rho_t"][t], 0.0))
        else:
            lo_errs.append(0.0)
            hi_errs.append(0.0)

    def sig_label(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return ""

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(x, rho_vals, width=0.6, color=[COLS5[t] for t in CATS5],
           alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[lo_errs, hi_errs],
                fmt="none", color="black", capsize=5, lw=1.2, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, t in enumerate(CATS5):
        sl = sig_label(pvals_bh[t])
        if sl:
            ax.text(x[i], rho_vals[i] + hi_errs[i] + 0.01, sl,
                    ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([MDLABELS[t] for t in CATS5], fontsize=9)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=10)
    ax.set_title(
        "Structural content of VSP spans\n"
        f"$N$ = {N} residue-positions  |  error bars = 95% CI (Poisson);  *** BH $p$ < 0.001",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# 10C: start-end asymmetry
# ---------------------------------------------------------------------------

def compute_asymmetry(boundaries, proteins):
    N         = len(boundaries)
    counts_s  = defaultdict(int)
    counts_e  = defaultdict(int)
    for bnd in boundaries:
        motifs = proteins[bnd["uid"]]["motifs"]
        counts_s[_tau5(bnd["can_start"], motifs)] += 1
        counts_e[_tau5(bnd["can_end"],   motifs)] += 1
    f_s  = {t: counts_s[t] / N for t in CATS5}
    f_e  = {t: counts_e[t] / N for t in CATS5}
    diff = {t: f_s[t] - f_e[t] for t in CATS5}
    L1   = sum(abs(diff[t]) for t in CATS5)
    return f_s, f_e, diff, L1


def permutation_asymmetry(boundaries, proteins, B, seed):
    N       = len(boundaries)
    rng     = np.random.default_rng(seed)
    cat_idx = {t: i for i, t in enumerate(CATS5)}

    bnd_elems = []
    for bnd in boundaries:
        motifs = proteins[bnd["uid"]]["motifs"]
        bnd_elems.append((
            cat_idx[_tau5(bnd["can_start"], motifs)],
            cat_idx[_tau5(bnd["can_end"],   motifs)],
        ))

    diff_perms = np.zeros((B, len(CATS5)))
    L1_perms   = np.zeros(B)
    for b in range(B):
        swap     = rng.random(N) < 0.5
        counts_s = np.zeros(len(CATS5), dtype=np.int64)
        counts_e = np.zeros(len(CATS5), dtype=np.int64)
        for i, (es, ee) in enumerate(bnd_elems):
            if swap[i]:
                counts_s[ee] += 1
                counts_e[es] += 1
            else:
                counts_s[es] += 1
                counts_e[ee] += 1
        d = (counts_s - counts_e) / N
        diff_perms[b] = d
        L1_perms[b]   = np.sum(np.abs(d))
    return diff_perms, L1_perms


def plot_asymmetry(f_s, f_e, diff_perms, L1_obs, L1_perms, pvals_bh, N, B, out):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    x  = np.arange(len(CATS5))
    w  = 0.35
    ax.bar(x - w/2, [f_s[t] for t in CATS5], w,
           color=[COLS5[t] for t in CATS5], alpha=0.85, label="VSP start", zorder=3)
    ax.bar(x + w/2, [f_e[t] for t in CATS5], w,
           color=[COLS5[t] for t in CATS5], alpha=0.45, label="VSP end",
           zorder=3, hatch="//")
    for i, t in enumerate(CATS5):
        p   = pvals_bh[t]
        sig = "**" if p < 0.01 else ("*" if p < 0.05 else "")
        if sig:
            ax.text(x[i], max(f_s[t], f_e[t]) + 0.02,
                    sig, ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([MDLABELS[t] for t in CATS5], fontsize=8)
    ax.set_ylabel("Fraction of VSP boundaries", fontsize=9)
    ax.set_title(f"Start vs end element distribution ($N$ = {N} VSPs)", fontsize=9)
    ax.legend(fontsize=8, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.18), ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax2 = axes[1]
    perm_p = (1 + np.sum(L1_perms >= L1_obs)) / (B + 1)
    lo, hi = np.percentile(L1_perms, [2.5, 97.5])
    ax2.hist(L1_perms, bins=40, color="#888888", alpha=0.6, edgecolor="white", zorder=2)
    ax2.axvspan(lo, hi, color="#888888", alpha=0.20, zorder=1, label="Null 95% interval")
    ax2.axvline(L1_obs, color="#C44E52", lw=2.0, zorder=3,
                label=f"Observed L1 = {L1_obs:.3f},  p = {perm_p:.4f}")
    ax2.set_xlabel("L1 norm of start-end difference", fontsize=9)
    ax2.set_ylabel("Permutation replicates", fontsize=9)
    ax2.set_title(f"Global asymmetry test (B = {B})", fontsize=9)
    ax2.legend(fontsize=8, frameon=False, loc="upper center",
               bbox_to_anchor=(0.5, -0.18), ncol=2)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Markdown helpers
# ---------------------------------------------------------------------------

def _sig(p):
    return "**" if p < 0.01 else ("*" if p < 0.05 else "ns")


def build_enrichment_md(N_t, f_t, pi_null, rho_t, pvals_raw, pvals_bh, z_scores, N, fig_path):
    rows = [
        f"| {MDLABELS[t]} | {N_t[t]} | {f_t[t]:.3f} | {pi_null[t]:.3f} | "
        f"{rho_t[t]:.3f} | {z_scores[t] ** 2:.3f} | "
        f"{pvals_raw[t]:.4f} | {pvals_bh[t]:.4f} | {_sig(pvals_bh[t])} |"
        for t in CATS5
    ]
    return (
        f"$N$ = {N} VSP boundary positions.\n\n"
        "| Element | $N_t$ | $f_t$ | $\\pi_t^0$ | $\\rho_t$ | "
        "$\\chi^2$ | Raw $p$ | BH $p$ | Sig |\n"
        "|---|---|---|---|---|---|---|---|---|\n"
        + "\n".join(rows)
        + "\n\nChi-square $p$-values are two-sided (z-score approximation); "
          "BH correction across 5 element types.\n\n"
        f"![]({fig_path})"
    )


def build_asymmetry_md(f_s, f_e, diff, diff_perms, L1_obs, L1_perms, pvals_bh, N, B, fig_path):
    perm_p = (1 + np.sum(L1_perms >= L1_obs)) / (B + 1)
    rows = [
        f"| {MDLABELS[t]} | {f_s[t]:.3f} | {f_e[t]:.3f} | {diff[t]:+.3f} | "
        f"{pvals_bh[t]:.4f} | {_sig(pvals_bh[t])} |"
        for t in CATS5
    ]
    reject = perm_p < 0.05
    verdict = (
        "$H_0$ **rejected**: start and end positions show significantly different "
        "element distributions."
        if reject else
        "$H_0$ **not rejected**: no significant overall asymmetry between "
        "start and end element distributions."
    )
    return (
        f"$N$ = {N} VSPs.\n\n"
        "| Element | $f_t^\\text{start}$ | $f_t^\\text{end}$ | Diff | BH $p$ | Sig |\n"
        "|---|---|---|---|---|---|\n"
        + "\n".join(rows)
        + f"\n\nGlobal L1 test: observed $L_1 = {L1_obs:.4f}$, "
          f"permutation $p = {perm_p:.4f}$ ($B = {B}$). " + verdict
        + f"\n\n![]({fig_path})"
    )


def build_interp_md(rho_s, pvals_s_bh, rho_e, pvals_e_bh, perm_p_global):
    sig_s = [t for t in CATS5 if pvals_s_bh[t] < 0.05]
    sig_e = [t for t in CATS5 if pvals_e_bh[t] < 0.05]
    lines = []

    if not sig_s and not sig_e:
        lines.append(
            "Neither VSP start nor end positions show significant enrichment or depletion "
            "relative to the length-weighted null after BH correction. "
            "VSP span boundaries are distributed across structural elements "
            "in rough proportion to element length, suggesting no strong structural "
            "preference for where AS-altered regions begin or end in the TIM-barrel domain."
        )
    else:
        if sig_s:
            parts = [
                f"{MDLABELS[t]} ($\\rho = {rho_s[t]:.3f}$, BH $p = {pvals_s_bh[t]:.4f}$)"
                for t in sig_s
            ]
            lines.append(
                f"VSP **start** positions show significant structural preference: "
                + "; ".join(parts) + "."
            )
        else:
            lines.append("VSP start positions show no significant structural preference.")
        if sig_e:
            parts = [
                f"{MDLABELS[t]} ($\\rho = {rho_e[t]:.3f}$, BH $p = {pvals_e_bh[t]:.4f}$)"
                for t in sig_e
            ]
            lines.append(
                f"VSP **end** positions show significant structural preference: "
                + "; ".join(parts) + "."
            )
        else:
            lines.append("VSP end positions show no significant structural preference.")

    if perm_p_global < 0.05:
        lines.append(
            "The global start-end asymmetry test is significant, indicating that VSP spans "
            "enter and exit the TIM-barrel domain at structurally distinct element types."
        )
    else:
        lines.append(
            "No significant global start-end asymmetry was detected: the element distributions "
            "of VSP start and end positions are not systematically different from each other."
        )

    lines.append(
        "These results complement §9 by characterising the structural context at the boundaries "
        "of AS-altered spans, rather than which canonical junctions fall within them. "
        "Together, §9 and §10 provide a two-level view: which exon boundaries are enclosed by "
        "VSP-defined AS events (§9), and where those events enter and exit the domain (§10)."
    )
    return "\n\n".join(lines)


# ---------------------------------------------------------------------------
# MD scaffold and update
# ---------------------------------------------------------------------------

SECTION10 = """
---

## §10. VSP span boundary placement in TIM-barrel structural elements

### Overview

§9 identified which canonical exon-intron junctions are enclosed within VSP-defined
AS-altered spans.  §10 asks the complementary question: **where in the domain structure
do VSP spans themselves begin and end?**  Each VSP feature marks the boundaries of the
canonical sequence region that is replaced or deleted in the alternatively spliced isoform;
these boundaries correspond to the actual splice points in protein-sequence space.

For each VSP with domain overlap, the domain-clipped start
$s_v = \\max(v_s,\\, d_s)$ and end $e_v = \\min(v_e,\\, d_e - 1)$ are extracted and
assigned to structural elements using the τ₅ classifier (§5).  The enrichment ratio

$$\\rho_t = \\frac{f_t}{\\pi_t^0}$$

compares the observed fraction $f_t$ of VSP boundaries in element type $t$ to the
length-weighted null $\\pi_t^0$ — the fraction of domain residues in element $t$ —
computed from the same 227 canonical proteins.  Statistical significance is assessed
with a chi-square z-score test: $z_t = (O_t - E_t)/\\sqrt{E_t}$ where $E_t = N \\pi_t^0$,
converted to a two-sided $p$-value and BH-corrected across the five element types.

---

### Dataset

| | |
|---|---|
| Canonical proteins (baseline) | 227 |
__DATASET__

---

### 10A.  Element distribution of VSP start positions

#### 10A Results

__10A__

---

### 10B.  Element distribution of VSP end positions

#### 10B Results

__10B__

---

### 10C.  Start-end asymmetry

For each of the $N$ VSPs, the structural elements at the start and end positions form a
paired comparison.  Under $H_0$: start and end positions are drawn from the same structural
distribution — no systematic asymmetry between where VSP spans enter and exit the domain.
The permutation null randomly and independently swaps start/end labels within each VSP,
preserving the overall marginal boundary-position distribution.

Per-element differences $f_t^\\text{start} - f_t^\\text{end}$ are tested against this null
with BH correction across five elements.  A global test uses the L1 norm
$\\sum_t |f_t^\\text{start} - f_t^\\text{end}|$ as the test statistic.

#### 10C Results

__10C__

---

### 10D.  Biological interpretation

__10D__

"""


def update_md(md_path, scaffold, dataset_md, r10a, r10b):
    text = Path(md_path).read_text(encoding="utf-8")

    # Remove existing §10 block if present
    marker = "\n---\n\n## §10."
    if marker in text:
        start = text.index(marker)
        ref_tag = "\n## References\n"
        if ref_tag in text:
            end = text.index(ref_tag, start)
            text = text[:start] + text[end:]
        else:
            text = text[:start]

    # Fill placeholders
    section = scaffold
    section = section.replace("__DATASET__", dataset_md)
    section = section.replace("__10A__",     r10a)
    section = section.replace("__10B__",     r10b)

    # Insert before References, or append
    ref_tag = "\n## References\n"
    if ref_tag in text:
        idx  = text.index(ref_tag)
        text = text[:idx] + section + text[idx:]
    else:
        text = text.rstrip() + "\n" + section

    Path(md_path).write_text(text, encoding="utf-8")
    print(f"Updated section 10 in {md_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",        default=None)
    parser.add_argument("--B",         type=int, default=2000,
                        help="Permutation replicates for §10C asymmetry test")
    parser.add_argument("--seed",      type=int, default=42,
                        help="RNG seed for §10C asymmetry permutation")
    parser.add_argument("--out-10a",      default="figures/vsp_start_enrichment.png")
    parser.add_argument("--out-10b",      default="figures/vsp_end_enrichment.png")
    parser.add_argument("--out-combined", default="figures/vsp_boundary_enrichment.png")
    parser.add_argument("--out-pooled",    default="figures/vsp_boundary_pooled.png")
    parser.add_argument("--out-coverage",  default="figures/vsp_residue_coverage.png")
    parser.add_argument("--md",        default="Statistical-Analysis.md")
    args = parser.parse_args()

    db_path    = args.db or get_config().db_path
    conn       = sqlite3.connect(db_path)
    proteins   = load_proteins(conn)
    iso_seqs   = load_isoform_sequences(conn)
    boundaries = load_vsp_boundaries(conn, proteins, iso_seqs)
    conn.close()

    n_prot    = len(proteins)
    n_vsps    = len(boundaries)
    n_covered = len(set(b["uid"] for b in boundaries))
    print(f"Canonical proteins: {n_prot}")
    print(f"VSP spans with domain overlap: {n_vsps}  ({n_covered} distinct proteins)")

    if n_vsps == 0:
        print("No VSP boundaries found -- nothing to analyse.")
        return

    print("Computing length-weighted null ...")
    pi_null = compute_length_null(proteins)
    print("   pi_t^0:", {t: f"{pi_null[t]:.3f}" for t in CATS5})

    # -- 10A: start positions ------------------------------------------------
    print("\n10A. VSP start positions (chi-square) ...")
    N_t_s, f_t_s, rho_t_s, pvals_s_raw, pvals_s_bh, z_scores_s = \
        chi_square_pvalues(boundaries, proteins, pi_null, "start")
    for t in CATS5:
        print(f"   {LABELS5[t]}: N={N_t_s[t]}, f={f_t_s[t]:.3f}, "
              f"pi={pi_null[t]:.3f}, rho={rho_t_s[t]:.3f}, "
              f"raw p={pvals_s_raw[t]:.4f}, BH p={pvals_s_bh[t]:.4f}")
    plot_enrichment(rho_t_s, n_vsps, pi_null, "start", pvals_s_bh, pvals_s_raw, args.out_10a)

    # -- 10B: end positions --------------------------------------------------
    print("\n10B. VSP end positions (chi-square) ...")
    N_t_e, f_t_e, rho_t_e, pvals_e_raw, pvals_e_bh, z_scores_e = \
        chi_square_pvalues(boundaries, proteins, pi_null, "end")
    for t in CATS5:
        print(f"   {LABELS5[t]}: N={N_t_e[t]}, f={f_t_e[t]:.3f}, "
              f"pi={pi_null[t]:.3f}, rho={rho_t_e[t]:.3f}, "
              f"raw p={pvals_e_raw[t]:.4f}, BH p={pvals_e_bh[t]:.4f}")
    plot_enrichment(rho_t_e, n_vsps, pi_null, "end", pvals_e_bh, pvals_e_raw, args.out_10b)

    # -- Combined start + end figure -----------------------------------------
    plot_combined_enrichment(rho_t_s, rho_t_e, n_vsps, pi_null,
                             pvals_s_bh, pvals_e_bh,
                             pvals_s_raw, pvals_e_raw, args.out_combined)

    # -- Pooled: start + end as one set --------------------------------------
    N_pooled = 2 * n_vsps
    pvals_pooled_raw = {}
    for t in CATS5:
        N_t_pool = N_t_s[t] + N_t_e[t]
        E_t = N_pooled * pi_null[t]
        z   = (N_t_pool - E_t) / np.sqrt(E_t) if pi_null[t] > 0 else 0.0
        from scipy.stats import norm as _norm
        pvals_pooled_raw[t] = float(2 * _norm.sf(abs(z)))
    pvals_pooled_bh = bh_correct(pvals_pooled_raw, CATS5)
    N_t_pool_d, _, rho_t_pool_d = plot_pooled_enrichment(
        N_t_s, N_t_e, N_pooled, pi_null, pvals_pooled_bh, pvals_pooled_raw, args.out_pooled
    )
    print("\nPooled boundary enrichment:")
    for t in CATS5:
        print(f"   {LABELS5[t]}: N={N_t_pool_d[t]}, rho={rho_t_pool_d[t]:.3f}, "
              f"BH p={pvals_pooled_bh[t]:.4f}")

    # -- VSP residue coverage (Ochoa-Leyva comparison) -----------------------
    from scipy.stats import norm as _norm2
    cov = compute_vsp_residue_coverage(boundaries, proteins, pi_null)
    if cov:
        cov_pvals_raw = {}
        for t in CATS5:
            E_t = cov["N"] * pi_null[t]
            z   = (cov["counts"][t] - E_t) / np.sqrt(E_t) if pi_null[t] > 0 else 0.0
            cov_pvals_raw[t] = float(2 * _norm2.sf(abs(z)))
        cov_pvals_bh = bh_correct(cov_pvals_raw, CATS5)
        print(f"\nVSP residue coverage (Ochoa-Leyva-style)  —  N = {cov['N']} residue-positions")
        print(f"  {'Element':<16}  {'Count':>7}  {'f_t':>7}  {'pi_t':>7}  {'rho':>6}  {'p_raw':>7}  {'p_BH':>7}  Sig")
        print("  " + "-"*16 + "  " + "  ".join(["-"*7]*5 + ["-"*7, "---"]))
        for t in CATS5:
            sig = "***" if cov_pvals_bh[t] < 0.001 else ("**" if cov_pvals_bh[t] < 0.01 else ("*" if cov_pvals_bh[t] < 0.05 else "ns"))
            print(f"  {LABELS5[t]:<16}  {cov['counts'][t]:>7}  {cov['f_t'][t]:>7.4f}  "
                  f"{pi_null[t]:>7.4f}  {cov['rho_t'][t]:>6.3f}  "
                  f"{cov_pvals_raw[t]:>7.4f}  {cov_pvals_bh[t]:>7.4f}  {sig}")
        plot_residue_coverage(cov, pi_null, cov_pvals_bh, args.out_coverage)

    # -- Build markdown ------------------------------------------------------
    dataset_md = (
        f"| VSP spans with domain overlap | {n_vsps} |\n"
        f"| Distinct proteins covered | {n_covered} |"
    )
    r10a_md = build_enrichment_md(N_t_s, f_t_s, pi_null, rho_t_s,
                                   pvals_s_raw, pvals_s_bh, z_scores_s,
                                   n_vsps, args.out_10a)
    r10b_md = build_enrichment_md(N_t_e, f_t_e, pi_null, rho_t_e,
                                   pvals_e_raw, pvals_e_bh, z_scores_e,
                                   n_vsps, args.out_10b)
    md_path = Path(args.md)
    if md_path.exists():
        update_md(str(md_path), SECTION10, dataset_md, r10a_md, r10b_md)
    else:
        print(f"  [note] {md_path} not found -- markdown not updated.")


if __name__ == "__main__":
    main()
