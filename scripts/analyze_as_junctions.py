#!/usr/bin/env python3
"""
AS junction placement analysis (§8 of Statistical-Analysis.md).

  8A. Element-level distribution of AS-affected canonical junctions vs baseline
  8B. Cross-protein positional distribution (KS test vs canonical baseline)
  8C. Within-protein hotspot reuse across multiple AS isoforms

Output figures:
  figures/as_element_enrichment.png
  figures/as_positional.png
  figures/as_hotspot.png

Usage:
    python scripts/analyze_as_junctions.py
    python scripts/analyze_as_junctions.py --B 5000
    python scripts/analyze_as_junctions.py --full-only
"""

import argparse
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config


# ---------------------------------------------------------------------------
# Category definitions (identical to analyze_junction_enrichment.py)
# ---------------------------------------------------------------------------

CATS5    = ["beta", "alpha", "inter", "loop", "flanking"]
LABELS5  = {"beta": "β-strand", "alpha": "α-helix", "inter": "Inter-motif",
            "loop": "Loop (β→α)", "flanking": "Flanking"}
ALABELS5 = {"beta": "beta-strand", "alpha": "alpha-helix", "inter": "Inter-motif",
            "loop": "Loop (b->a)", "flanking": "Flanking"}
COLS5    = {"beta": "#4C72B0", "alpha": "#DD8452", "inter": "#55A868",
            "loop": "#8c8c8c", "flanking": "#C44E52"}


# ---------------------------------------------------------------------------
# Element assignment — τ_5, identical to analyze_junction_enrichment.py
# ---------------------------------------------------------------------------

def _tau5(pos, motifs):
    """5-category label for pos; first-match, lower motif index wins."""
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
# Load canonical proteins (same filter as other analysis scripts)
# ---------------------------------------------------------------------------

def load_proteins(conn, full_only=False):
    """Return dict uid -> protein record (only proteins with >=1 domain-internal junction)."""
    rows = conn.execute("""
        SELECT uniprot_id, gene_name, domain_start, domain_end,
               exon_annotations, motif_annotations
        FROM   canonical_analysis
        WHERE  exon_annotations  IS NOT NULL
          AND  motif_annotations IS NOT NULL
          AND  domain_start      IS NOT NULL
          AND  domain_end        IS NOT NULL
    """).fetchall()
    proteins = {}
    for uid, gene, ds, de, ea, ma in rows:
        motifs = json.loads(ma)
        if full_only and len(motifs) != 8:
            continue
        exons     = json.loads(ea)
        junctions = [e["end"] for e in exons[:-1] if ds <= e["end"] < de]
        if not junctions:
            continue
        proteins[uid] = dict(uid=uid, gene=gene, ds=ds, de=de,
                             motifs=motifs, n_motifs=len(motifs),
                             junctions=junctions)
    return proteins


# ---------------------------------------------------------------------------
# Load AS isoforms
# ---------------------------------------------------------------------------

def load_as_isoforms(conn, proteins):
    """
    For each entry in affected_isoforms whose canonical protein is in proteins,
    compute J_a^AS = canonical junctions inside any VSP span (union over all VSPs).
    Only isoforms with at least one AS junction are returned.
    """
    rows = conn.execute("""
        SELECT uniprot_id, isoform_id, vsp_domain_events
        FROM   affected_isoforms
        WHERE  vsp_domain_events IS NOT NULL
    """).fetchall()

    isoforms = []
    for uid, iso_id, vsp_json in rows:
        if uid not in proteins:
            continue
        vsps = json.loads(vsp_json)
        if not vsps:
            continue
        p = proteins[uid]
        as_junctions = set()
        for vsp in vsps:
            # can_start/can_end: VSP canonical span (full canonical sequence coordinates).
            # overlap_start/overlap_end in the same record refer to the domain-overlap
            # subinterval — do not use them here.
            vs = vsp.get("can_start")
            ve = vsp.get("can_end")
            if vs is None or ve is None:
                continue
            for j in p["junctions"]:
                if vs <= j <= ve:
                    as_junctions.add(j)
        as_junctions = sorted(as_junctions)
        if not as_junctions:
            continue
        isoforms.append(dict(uid=uid, isoform_id=iso_id,
                             as_junctions=as_junctions,
                             n_as=len(as_junctions)))
    return isoforms


# ---------------------------------------------------------------------------
# 8A — Element-level distribution
# ---------------------------------------------------------------------------

def compute_8a(proteins, isoforms):
    """
    Returns N_t_AS, f_t_AS, rho_t_AS, baseline f_t, canonical N, AS N.
    """
    N_total = 0
    N_t = defaultdict(int)
    for p in proteins.values():
        for j in p["junctions"]:
            N_t[_tau5(j, p["motifs"])] += 1
            N_total += 1
    f_t = {t: N_t[t] / N_total for t in CATS5}

    N_AS_total = 0
    N_t_AS = defaultdict(int)
    for iso in isoforms:
        p = proteins[iso["uid"]]
        for j in iso["as_junctions"]:
            N_t_AS[_tau5(j, p["motifs"])] += 1
            N_AS_total += 1

    f_t_AS  = {t: N_t_AS[t] / N_AS_total for t in CATS5} if N_AS_total else {}
    rho_t_AS = {t: (f_t_AS[t] / f_t[t]) if f_t.get(t, 0) > 0 else float("nan")
                for t in CATS5}
    return N_t_AS, f_t_AS, rho_t_AS, f_t, N_total, N_AS_total


def permutation_8a(proteins, isoforms, B, seed):
    """
    For each replicate, draw |J_a^AS| junctions uniformly from J_{p(a)} without
    replacement per isoform.  Returns rho_perms of shape (B, 5).
    """
    N_total = 0
    N_t = defaultdict(int)
    for p in proteins.values():
        for j in p["junctions"]:
            N_t[_tau5(j, p["motifs"])] += 1
            N_total += 1
    f_t = {t: N_t[t] / N_total for t in CATS5}
    cat_idx = {t: i for i, t in enumerate(CATS5)}

    # Cache per-protein type lists for speed
    uid_to_types = {uid: [_tau5(j, p["motifs"]) for j in p["junctions"]]
                    for uid, p in proteins.items()}

    iso_specs  = [(iso["uid"], iso["n_as"]) for iso in isoforms]
    N_AS_total = sum(k for _, k in iso_specs)

    rng = np.random.default_rng(seed)
    rho_perms = np.zeros((B, len(CATS5)))
    for b in range(B):
        counts = np.zeros(len(CATS5), dtype=np.int64)
        for uid, k in iso_specs:
            types_list = uid_to_types[uid]
            k_draw = min(k, len(types_list))
            for idx in rng.choice(len(types_list), size=k_draw, replace=False):
                counts[cat_idx[types_list[idx]]] += 1
        f_perm = counts / N_AS_total
        for ci, t in enumerate(CATS5):
            rho_perms[b, ci] = (f_perm[ci] / f_t[t]) if f_t[t] > 0 else float("nan")
    return rho_perms


def plot_8a(rho_t_AS, rho_perms, pvals_bh, n_iso, N_AS, out):
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(CATS5))
    rho_vals = [rho_t_AS.get(t, 1.0) for t in CATS5]
    null_lo  = [np.percentile(rho_perms[:, i],  2.5) for i in range(len(CATS5))]
    null_hi  = [np.percentile(rho_perms[:, i], 97.5) for i in range(len(CATS5))]
    null_mid = [(lo + hi) / 2 for lo, hi in zip(null_lo, null_hi)]
    null_err = [(hi - lo) / 2 for lo, hi in zip(null_lo, null_hi)]

    ax.bar(x, rho_vals, width=0.65, color=[COLS5[t] for t in CATS5],
           alpha=0.85, zorder=3)
    ax.errorbar(x, null_mid, yerr=null_err,
                fmt="none", color="black", capsize=5, lw=1.2, zorder=4,
                label="Null 95% interval")
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, t in enumerate(CATS5):
        p = pvals_bh.get(t, float("nan"))
        sig = ("**" if p < 0.01 else ("*" if p < 0.05 else "")) if not np.isnan(p) else ""
        if sig:
            ax.text(x[i], max(rho_vals[i], null_hi[i]) + 0.04,
                    sig, ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[t] for t in CATS5], fontsize=9)
    ax.set_ylabel(r"$\rho_t^{AS}$ (AS / canonical fraction)", fontsize=10)
    ax.set_title(
        f"AS junction element enrichment relative to canonical baseline\n"
        f"{n_iso} isoforms, $N^{{AS}}$ = {N_AS}  |  "
        "error bars = null 95% interval;  ** BH p < 0.01,  * BH p < 0.05",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# 8B — Cross-protein positional distribution
# ---------------------------------------------------------------------------

def normalized_positions_as(proteins, isoforms):
    xs = []
    for iso in isoforms:
        p = proteins[iso["uid"]]
        ds, de = p["ds"], p["de"]
        dlen = de - ds
        for j in iso["as_junctions"]:
            xs.append((j - ds) / dlen)
    return np.array(xs)


def normalized_positions_canonical(proteins):
    xs = []
    for p in proteins.values():
        ds, de = p["ds"], p["de"]
        dlen = de - ds
        for j in p["junctions"]:
            xs.append((j - ds) / dlen)
    return np.array(xs)


def ks_two_sample_vs_fixed(x_test, x_ref_sorted):
    """D = sup|F_test(x) - F_ref(x)| with F_ref kept as the empirical CDF of x_ref."""
    n_t = len(x_test)
    n_r = len(x_ref_sorted)
    if n_t == 0 or n_r == 0:
        return 0.0
    combined = np.sort(np.concatenate([x_test, x_ref_sorted]))
    f_test = np.searchsorted(np.sort(x_test),  combined, side="right") / n_t
    f_ref  = np.searchsorted(x_ref_sorted,      combined, side="right") / n_r
    return float(np.max(np.abs(f_test - f_ref)))


def permutation_8b(proteins, isoforms, B, seed, x_ref_sorted):
    """
    Same draw logic as 8A; return per-replicate KS statistics vs x_ref_sorted.
    """
    uid_to_info = {}
    for uid, p in proteins.items():
        ds, de = p["ds"], p["de"]
        jarr   = np.array(p["junctions"], dtype=float)
        uid_to_info[uid] = (jarr, ds, de - ds)

    iso_specs = [(iso["uid"], iso["n_as"]) for iso in isoforms]
    rng = np.random.default_rng(seed)
    D_perms = np.zeros(B)
    for b in range(B):
        xs = []
        for uid, k in iso_specs:
            jarr, ds, dlen = uid_to_info[uid]
            k_draw = min(k, len(jarr))
            drawn  = rng.choice(len(jarr), size=k_draw, replace=False)
            xs.extend((jarr[drawn] - ds) / dlen)
        D_perms[b] = ks_two_sample_vs_fixed(np.array(xs), x_ref_sorted)
    return D_perms


def mean_element_spans(proteins):
    spans = defaultdict(list)
    for p in proteins.values():
        ds, de = p["ds"], p["de"]
        dlen = de - ds
        if dlen <= 0:
            continue
        for i, m in enumerate(p["motifs"]):
            n = m["motif"]
            def r(x, ds=ds, dlen=dlen): return (x - ds) / dlen
            spans[f"beta_{n}"].append( (r(m["beta_start"]),      r(m["beta_end"])) )
            spans[f"loop_{n}"].append( (r(m["beta_end"] + 1),    r(m["alpha_start"] - 1)) )
            spans[f"alpha_{n}"].append((r(m["alpha_start"]),      r(m["alpha_end"])) )
            if i + 1 < len(p["motifs"]):
                nxt = p["motifs"][i + 1]
                spans[f"inter_{n}_{n+1}"].append(
                    (r(m["alpha_end"] + 1), r(nxt["beta_start"] - 1)))
    return {k: (np.mean([s for s, e in v]), np.mean([e for s, e in v]))
            for k, v in spans.items() if v}


def plot_8b(x_as, x_const, D_obs, perm_p, spans, n_iso, B, out):
    bw   = "scott"
    grid = np.linspace(0, 1, 500)
    BG   = {"beta":  ("#4C72B0", 0.18), "alpha": ("#DD8452", 0.18),
            "loop":  ("#eeeeee", 0.80), "inter": ("#eeeeee", 0.80)}

    kde_as    = stats.gaussian_kde(x_as,    bw_method=bw)(grid) if len(x_as) > 1 else np.ones(500)
    kde_const = stats.gaussian_kde(x_const, bw_method=bw)(grid)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    for elem, (s, e) in spans.items():
        if e <= s:
            continue
        prefix = elem.split("_")[0]
        col, alpha = BG.get(prefix, ("#eeeeee", 0.5))
        ax.axvspan(s, e, color=col, alpha=alpha, zorder=0)

    ax.plot(grid, kde_const, color="#4C72B0", lw=1.2, ls="--", zorder=3,
            label=f"Canonical junctions (N = {len(x_const)})")
    ax.plot(grid, kde_as,    color="black",   lw=1.8, zorder=4,
            label=f"AS-affected junctions (N = {len(x_as)})")

    bg_patches = [
        mpatches.Patch(facecolor="#4C72B0", alpha=0.40, label="β-strand region (mean)"),
        mpatches.Patch(facecolor="#DD8452", alpha=0.40, label="α-helix region (mean)"),
        mpatches.Patch(facecolor="#cccccc", alpha=0.80, label="Loop / inter-motif (mean)"),
    ]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + bg_patches, fontsize=8, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3)

    ax.set_xlabel("Normalised domain position", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_title(
        f"Positional distribution: AS-affected vs canonical junctions\n"
        f"{n_iso} isoforms  |  KS D = {D_obs:.4f},  permutation p = {perm_p:.4f}  (B = {B})",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# 8C — Within-protein hotspot reuse
# ---------------------------------------------------------------------------

def compute_8c(proteins, isoforms):
    """
    Returns H_bar, multi-isoform protein dict, per-protein stats,
    total hotspot count, total unique junction count.
    """
    by_protein = defaultdict(list)
    for iso in isoforms:
        by_protein[iso["uid"]].append(iso["as_junctions"])

    multi = {uid: jsets for uid, jsets in by_protein.items() if len(jsets) >= 2}

    total_hotspot = 0
    total_unique  = 0
    per_protein   = {}
    for uid, jsets in multi.items():
        usage = defaultdict(int)
        for jlist in jsets:
            for j in jlist:
                usage[j] += 1
        union_j   = set(usage)
        hotspot_j = {j for j, u in usage.items() if u >= 2}
        H_p = len(hotspot_j) / len(union_j) if union_j else 0.0
        per_protein[uid] = dict(H_p=H_p, hotspot=len(hotspot_j),
                                unique=len(union_j), k=len(jsets))
        total_hotspot += len(hotspot_j)
        total_unique  += len(union_j)

    H_bar = total_hotspot / total_unique if total_unique > 0 else 0.0
    return H_bar, multi, per_protein, total_hotspot, total_unique


def permutation_8c(proteins, isoforms, B, seed):
    """
    Permutation null for H_bar: each isoform draws |J_a^AS| junctions from its
    canonical pool without replacement.
    """
    by_protein = defaultdict(list)
    for iso in isoforms:
        by_protein[iso["uid"]].append(iso["n_as"])

    multi = {uid: ks for uid, ks in by_protein.items() if len(ks) >= 2}
    if not multi:
        return np.zeros(B)

    uid_to_pool = {uid: proteins[uid]["junctions"] for uid in multi}
    rng = np.random.default_rng(seed)
    H_bar_perms = np.zeros(B)
    for b in range(B):
        total_hotspot = 0
        total_unique  = 0
        for uid, ks in multi.items():
            pool   = uid_to_pool[uid]
            n_pool = len(pool)
            usage  = defaultdict(int)
            for k in ks:
                k_draw = min(k, n_pool)
                for idx in rng.choice(n_pool, size=k_draw, replace=False):
                    usage[pool[idx]] += 1
            union_j   = set(usage)
            hotspot_j = {j for j, u in usage.items() if u >= 2}
            total_hotspot += len(hotspot_j)
            total_unique  += len(union_j)
        H_bar_perms[b] = total_hotspot / total_unique if total_unique > 0 else 0.0
    return H_bar_perms


def plot_8c(H_bar, H_bar_perms, perm_p, n_multi, total_hotspot, total_unique, B, out):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(H_bar_perms, bins=40, color="#888888", alpha=0.6,
            edgecolor="white", zorder=2, label="Permutation null")
    lo, hi = np.percentile(H_bar_perms, [2.5, 97.5])
    ax.axvspan(lo, hi, color="#888888", alpha=0.20, zorder=1, label="Null 95% interval")
    ax.axvline(H_bar, color="#C44E52", lw=2.0, zorder=3,
               label=f"Observed H-bar = {H_bar:.3f}")
    ax.set_xlabel(r"Hotspot fraction $\bar{H}$", fontsize=10)
    ax.set_ylabel("Permutation replicates", fontsize=10)
    ax.set_title(
        f"Within-protein AS hotspot reuse  ({n_multi} proteins with k>=2 isoforms)\n"
        f"{total_hotspot} hotspot / {total_unique} unique junctions  |  "
        f"permutation p = {perm_p:.4f}  (B = {B})",
        fontsize=9,
    )
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# BH correction
# ---------------------------------------------------------------------------

def bh_correct(pvals_dict, keys):
    pvals = np.array([pvals_dict[k] for k in keys])
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.zeros(m)
    for rank, idx in enumerate(order, 1):
        adj[idx] = pvals[idx] * m / rank
    for j in range(m - 2, -1, -1):
        adj[order[j]] = min(adj[order[j]], adj[order[j + 1]])
    adj = np.minimum(adj, 1.0)
    return {k: float(adj[i]) for i, k in enumerate(keys)}


# ---------------------------------------------------------------------------
# Update §8 in Statistical-Analysis.md
# ---------------------------------------------------------------------------

def _replace_block(text, start_tag, end_tag, new_content):
    si = text.find(start_tag)
    if si == -1:
        return text
    content_start = si + len(start_tag)
    if end_tag is None:
        return text[:content_start] + "\n" + new_content + "\n"
    ei = text.find(end_tag, content_start)
    if ei == -1:
        return text[:content_start] + "\n" + new_content + "\n"
    return text[:content_start] + "\n" + new_content + "\n" + text[ei:]


def update_section9(md_path, content_9a, content_9b, content_9c):
    text = Path(md_path).read_text(encoding="utf-8")
    text = _replace_block(text, "#### 9A Results\n", "\n---\n\n### 9B.", content_9a)
    text = _replace_block(text, "#### 9B Results\n", "\n---\n\n### 9C.", content_9b)
    text = _replace_block(text, "#### 9C Results\n", "\n---\n\n### 9D.", content_9c)
    Path(md_path).write_text(text, encoding="utf-8")
    print(f"Updated §9 in {md_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",        default=None)
    parser.add_argument("--B",         type=int, default=2000)
    parser.add_argument("--seed",      type=int, default=42)
    parser.add_argument("--full-only", action="store_true")
    parser.add_argument("--out-8a",    default="figures/as_element_enrichment.png")
    parser.add_argument("--out-8b",    default="figures/as_positional.png")
    parser.add_argument("--out-8c",    default="figures/as_hotspot.png")
    parser.add_argument("--md",        default="Statistical-Analysis.md")
    args = parser.parse_args()

    db_path  = args.db or get_config().db_path
    conn     = sqlite3.connect(db_path)
    proteins = load_proteins(conn, full_only=args.full_only)
    isoforms = load_as_isoforms(conn, proteins)
    conn.close()

    n_prot  = len(proteins)
    n_iso   = len(isoforms)
    N_AS    = sum(iso["n_as"] for iso in isoforms)
    n_genes = len(set(iso["uid"] for iso in isoforms))
    print(f"Loaded {n_prot} canonical proteins.")
    print(f"AS isoforms with >=1 domain junction in VSP span: {n_iso} "
          f"({n_genes} distinct canonical proteins), {N_AS} AS-affected junction instances.")

    if n_iso == 0:
        print("No AS isoforms found -- nothing to analyse.")
        return

    # ── 8A. Element distribution ─────────────────────────────────────────────
    print(f"\n8A. Element enrichment (B = {args.B}) ...")
    N_t_AS, f_t_AS, rho_t_AS, f_t, N_total, N_AS_total = compute_8a(proteins, isoforms)
    for t in CATS5:
        print(f"   {ALABELS5[t]}: N_t_AS={N_t_AS[t]}, "
              f"f_t_AS={f_t_AS.get(t,0):.3f}, f_t={f_t[t]:.3f}, "
              f"rho={rho_t_AS.get(t,float('nan')):.3f}")

    print("   Permutation ...")
    rho_perms = permutation_8a(proteins, isoforms, args.B, args.seed)

    pvals_8a = {}
    for ci, t in enumerate(CATS5):
        rho_obs = rho_t_AS.get(t, 1.0)
        p_val   = (1 + np.sum(np.abs(rho_perms[:, ci] - 1) >= abs(rho_obs - 1))) / (args.B + 1)
        pvals_8a[t] = float(p_val)
    pvals_8a_bh = bh_correct(pvals_8a, CATS5)
    for t in CATS5:
        print(f"   {ALABELS5[t]}: raw p={pvals_8a[t]:.4f}, BH p={pvals_8a_bh[t]:.4f}")

    plot_8a(rho_t_AS, rho_perms, pvals_8a_bh, n_iso, N_AS_total, args.out_8a)

    # ── 8B. Positional distribution ──────────────────────────────────────────
    print(f"\n8B. Positional distribution (B = {args.B}) ...")
    x_as    = normalized_positions_as(proteins, isoforms)
    x_const = normalized_positions_canonical(proteins)
    x_const_sorted = np.sort(x_const)
    D_obs   = ks_two_sample_vs_fixed(x_as, x_const_sorted)
    print(f"   KS D = {D_obs:.4f}")

    D_perms   = permutation_8b(proteins, isoforms, args.B, args.seed, x_const_sorted)
    perm_p_8b = (1 + np.sum(D_perms >= D_obs)) / (args.B + 1)
    print(f"   Permutation p = {perm_p_8b:.4f}")

    spans = mean_element_spans(proteins)
    plot_8b(x_as, x_const, D_obs, perm_p_8b, spans, n_iso, args.B, args.out_8b)

    # ── 8C. Hotspot reuse ────────────────────────────────────────────────────
    print(f"\n8C. Hotspot reuse (B = {args.B}) ...")
    H_bar, multi, per_protein, total_hotspot, total_unique = compute_8c(proteins, isoforms)
    n_multi = len(multi)
    print(f"   Multi-isoform proteins (k>=2): {n_multi}")
    print(f"   H_bar = {H_bar:.4f}  ({total_hotspot} hotspot / {total_unique} unique junctions)")

    H_bar_perms = permutation_8c(proteins, isoforms, args.B, args.seed)
    perm_p_8c   = (1 + np.sum(H_bar_perms >= H_bar)) / (args.B + 1)
    print(f"   Permutation p = {perm_p_8c:.4f}")

    plot_8c(H_bar, H_bar_perms, perm_p_8c, n_multi, total_hotspot, total_unique,
            args.B, args.out_8c)

    # ── Build markdown content ────────────────────────────────────────────────
    def sig(p):
        return "**" if p < 0.01 else ("*" if p < 0.05 else "ns")

    null_lo_r = {t: np.percentile(rho_perms[:, ci], 2.5)  for ci, t in enumerate(CATS5)}
    null_hi_r = {t: np.percentile(rho_perms[:, ci], 97.5) for ci, t in enumerate(CATS5)}

    rows_8a = [
        f"| {LABELS5[t]} | {N_t_AS[t]} | {f_t_AS.get(t,0):.3f} | {f_t[t]:.3f} | "
        f"{rho_t_AS.get(t, float('nan')):.3f} | "
        f"[{null_lo_r[t]:.3f}, {null_hi_r[t]:.3f}] | "
        f"{pvals_8a[t]:.4f} | {pvals_8a_bh[t]:.4f} | {sig(pvals_8a_bh[t])} |"
        for t in CATS5
    ]

    content_8a_md = (
        f"Dataset: $|\\mathcal{{A}}|$ = {n_iso} isoforms across {n_genes} canonical proteins; "
        f"$N^{{AS}}$ = {N_AS_total} AS-affected junction instances.\n\n"
        "| Element | $N_t^{AS}$ | $f_t^{AS}$ | $f_t$ | $\\rho_t^{AS}$ | "
        "Null 95% interval | Raw $p$ | BH $p$ | Sig |\n"
        "|---|---|---|---|---|---|---|---|---|\n"
        + "\n".join(rows_8a)
        + "\n\nPermutation $p$-values are two-sided; BH correction across all 5 element types.\n\n"
        "![AS junction element enrichment](figures/as_element_enrichment.png)"
    )

    reject_8b = perm_p_8b < 0.05
    content_8b_md = (
        f"$N^{{AS}}$ = {len(x_as)} normalised positions from {n_iso} isoforms; "
        f"canonical baseline $N$ = {len(x_const)}.\n\n"
        "| Statistic | Value |\n"
        "|---|---|\n"
        f"| KS statistic $D_N^{{AS}}$ | {D_obs:.4f} |\n"
        f"| Permutation $p$ ($B = {args.B}$) | {perm_p_8b:.4f} |\n\n"
        f"$H_0^{{AS}}$ **{'rejected' if reject_8b else 'not rejected'}** "
        f"(permutation $p = {perm_p_8b:.4f}$, $\\alpha = 0.05$). "
        + (
            "The AS-affected junction positional distribution is significantly different "
            "from the canonical baseline."
            if reject_8b else
            "The AS-affected junction positional distribution is not significantly different "
            "from the canonical baseline; AS-affected canonical junctions do not cluster at "
            "different domain positions than the full canonical pool."
        )
        + "\n\n![AS junction positional distribution](figures/as_positional.png)"
    )

    null_H_lo, null_H_hi = np.percentile(H_bar_perms, [2.5, 97.5])
    reject_8c = perm_p_8c < 0.05
    content_8c_md = (
        f"Multi-isoform proteins ($k_p \\ge 2$): $|\\mathcal{{P}}^{{(2)}}|$ = {n_multi}.\n\n"
        "| Statistic | Value |\n"
        "|---|---|\n"
        f"| Proteins with $k_p \\ge 2$ | {n_multi} |\n"
        f"| Total unique AS-affected junctions | {total_unique} |\n"
        f"| Hotspot junctions ($u_p(j) \\ge 2$) | {total_hotspot} |\n"
        f"| Observed $\\bar{{H}}$ | {H_bar:.4f} |\n"
        f"| Null 95% interval | [{null_H_lo:.4f}, {null_H_hi:.4f}] |\n"
        f"| Permutation $p$ ($B = {args.B}$) | {perm_p_8c:.4f} |\n\n"
        f"$H_0^{{AS}}$ **{'rejected' if reject_8c else 'not rejected'}** "
        f"(permutation $p = {perm_p_8c:.4f}$, $\\alpha = 0.05$). "
        + (
            "The observed hotspot fraction significantly exceeds the permutation null, indicating "
            "that multiple AS isoforms of the same protein non-randomly overlap the same canonical "
            "junctions — consistent with within-protein AS hotspots."
            if reject_8c else
            "The observed hotspot fraction is within the permutation null band; "
            "within-protein AS hotspot reuse is not significantly above chance given the "
            "canonical junction pool sizes."
        )
        + "\n\n![Within-protein hotspot reuse](figures/as_hotspot.png)"
    )

    md_path = Path(args.md)
    if md_path.exists():
        update_section9(str(md_path), content_8a_md, content_8b_md, content_8c_md)
    else:
        print(f"  [note] {md_path} not found — markdown not updated.")


if __name__ == "__main__":
    main()
