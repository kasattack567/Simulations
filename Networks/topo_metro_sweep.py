"""
Metro-span sweep on the 6 real optical topologies — coverage-first, balanced.

For each real topology rescaled to spans 10-200 km, we compute DV vs CV coverage
and average key rate (routed along the real edges, each protocol using its own
reach: DV ~275 km, CV ~40 km).

CRITERION (stated explicitly): full coverage is treated as a hard requirement.
A network that cannot connect every user pair has failed regardless of key rate,
so the rate comparison is only meaningful while a protocol maintains 100%
coverage. We therefore report each topology's CV FULL-COVERAGE LIMIT — the span
beyond which CV can no longer connect the whole network — and compare rate only
up to that limit.

This framing is even-handed: coverage-first FAVOURS DV (reach ~275 km keeps DV
at 100% far longer), while rate FAVOURS CV (higher per-link rate within reach).
Neither dominates; the figure shows both regimes.

CV detection: heterodyne is the CV representative in the summary panels
(homodyne tracks it within a few percent, per the point-to-point results);
homodyne coverage is also computed so any divergence near the limit is visible.

Three panels:
  1) Coverage vs span (all 6; CV solid, DV dashed) — DV's reach advantage.
  2) CV/DV rate ratio vs span (all 6), greyed beyond each CV coverage limit —
     CV's rate advantage, shown only where CV is a valid full-coverage network.
  3) CV full-coverage span limit per topology (bar) — the headline number.

Usage:
    python topo_metro_sweep.py
"""
import argparse
import os
import heapq
import numpy as np
import matplotlib.pyplot as plt

from net_common import link_rate, CV_REACH_KM, DV_REACH_KM
from topo_loader import TOPOLOGY_FILES, load_topology, edge_lengths

REACH = {"dv": DV_REACH_KM, "cv_het": CV_REACH_KM, "cv_hom": CV_REACH_KM}

TCOLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]


def widest_paths_from(src, n, adj):
    best = [0.0] * n
    best[src] = float("inf")
    pq = [(-float("inf"), src)]
    while pq:
        negc, u = heapq.heappop(pq)
        c = -negc
        if c < best[u]:
            continue
        for v, w in adj[u]:
            nc = min(c, w)
            if nc > best[v]:
                best[v] = nc
                heapq.heappush(pq, (-nc, v))
    return best


def rate_and_coverage(n_nodes, edges, elen, protocol):
    reach = REACH[protocol]
    adj = [[] for _ in range(n_nodes)]
    for (a, b), L in zip(edges, elen):
        if L <= reach:
            r = link_rate(L, protocol)
            if r > 1e-12:
                adj[a].append((b, r)); adj[b].append((a, r))
    total, connected = 0.0, 0
    n_pairs = n_nodes * (n_nodes - 1) // 2
    for s in range(n_nodes):
        caps = widest_paths_from(s, n_nodes, adj)
        for t in range(s + 1, n_nodes):
            if caps[t] not in (0.0, float("inf")):
                total += caps[t]; connected += 1
    avg = total / n_pairs if n_pairs else 0.0
    cov = connected / n_pairs if n_pairs else 0.0
    return avg, cov


def cv_coverage_limit(spans, cv_cov):
    """Largest span at which CV still has full (100%) coverage."""
    full = [s for s, c in zip(spans, cv_cov) if c >= 0.999]
    return max(full) if full else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--spans", type=float, nargs="+",
                   default=list(range(10, 301, 5)))
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    spans = np.array(args.spans, dtype=float)
    names = list(TOPOLOGY_FILES.keys())

    data = {}
    print("CV full-coverage limit per topology "
          "(span beyond which CV cannot connect the whole network):\n")
    for name in names:
        topo = load_topology(name)
        n = len(topo["node_ids"]); edges = topo["edges"]
        dv_r, cv_r, dv_c, cv_c, hom_c = [], [], [], [], []
        for span in spans:
            elen = edge_lengths(topo, "metro", span)
            a_dv, c_dv = rate_and_coverage(n, edges, elen, "dv")
            a_cv, c_cv = rate_and_coverage(n, edges, elen, "cv_het")
            _, c_hom = rate_and_coverage(n, edges, elen, "cv_hom")
            dv_r.append(a_dv); cv_r.append(a_cv)
            dv_c.append(c_dv); cv_c.append(c_cv); hom_c.append(c_hom)
        data[name] = dict(dv_r=np.array(dv_r), cv_r=np.array(cv_r),
                          dv_c=np.array(dv_c), cv_c=np.array(cv_c),
                          hom_c=np.array(hom_c), n=n)
        lim = cv_coverage_limit(spans, cv_c)
        data[name]["limit"] = lim
        print(f"  {name:12} ({n:3d} nodes): CV full coverage up to "
              f"{'~%.0f km' % lim if lim else '<10 km'}")

    plt.rcParams.update({"font.family": "serif", "font.size": 11})
    fig = plt.figure(figsize=(18, 6.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.9], wspace=0.28)
    axC = fig.add_subplot(gs[0]); axR = fig.add_subplot(gs[1]); axB = fig.add_subplot(gs[2])

    # Panel 1: coverage vs span
    for name, col in zip(names, TCOLORS):
        d = data[name]
        axC.plot(spans, d["cv_c"] * 100, "-", color=col, lw=2, label=name)
        axC.plot(spans, d["dv_c"] * 100, "--", color=col, lw=1.2, alpha=0.6)
    axC.axvspan(10, 40, color="#7cc47f", alpha=0.10, zorder=0)
    axC.set_xlabel("Network span (km)")
    axC.set_ylabel("Coverage (%)")
    axC.set_ylim(0, 105)
    axC.set_title("Coverage vs span\n(solid = CV, dashed = DV)", fontsize=12)
    axC.grid(True, alpha=0.3)
    axC.legend(fontsize=8, ncol=2, loc="lower left")

    # Panel 2: CV/DV rate ratio, greyed past CV coverage limit
    for name, col in zip(names, TCOLORS):
        d = data[name]
        ratio = np.divide(d["cv_r"], d["dv_r"],
                          out=np.zeros_like(d["cv_r"]), where=d["dv_r"] > 0)
        lim = d["limit"]
        if lim is not None:
            valid = spans <= lim
            axR.plot(spans[valid], ratio[valid], "-", color=col, lw=2, label=name)
            axR.plot(spans[~valid], ratio[~valid], ":", color=col, lw=1, alpha=0.35)
        else:
            axR.plot(spans, ratio, ":", color=col, lw=1, alpha=0.35, label=name)
    axR.axhline(1.0, color="0.3", lw=1.2, ls="-")
    axR.text(spans[-1], 1.0, " CV=DV", fontsize=8, va="bottom", ha="right", color="0.3")
    axR.axvspan(10, 40, color="#7cc47f", alpha=0.10, zorder=0)
    axR.set_yscale("log")
    axR.set_xlabel("Network span (km)")
    axR.set_ylabel("CV / DV  average rate ratio")
    axR.set_title("Rate advantage (CV/DV)\nsolid = CV at full coverage, dotted = CV degraded",
                  fontsize=12)
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=8, ncol=2, loc="upper right")

    # Panel 3: CV full-coverage span limit per topology (bar)
    limits = [data[n]["limit"] if data[n]["limit"] else 0 for n in names]
    order = np.argsort(limits)
    bnames = [names[i] for i in order]
    blims = [limits[i] for i in order]
    bcols = [TCOLORS[i] for i in order]
    axB.barh(range(len(bnames)), blims, color=bcols)
    axB.set_yticks(range(len(bnames))); axB.set_yticklabels(bnames, fontsize=9)
    axB.axvspan(10, 40, color="#7cc47f", alpha=0.12, zorder=0)
    axB.set_xlabel("CV full-coverage span limit (km)")
    axB.set_title("How far CV stays a\nfull-coverage network", fontsize=12)
    axB.grid(True, axis="x", alpha=0.3)
    for i, v in enumerate(blims):
        axB.text(v, i, f" {v:.0f}", va="center", fontsize=9)

    fig.suptitle("CV vs DV on 6 real topologies - coverage limits and rate advantage",
                 fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or "topo_metro_sweep.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()