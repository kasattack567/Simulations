"""
topo_direct.py -- Direct-link (no relays) per-topology coverage and key-rate.

For each of the six selected topologies at REAL scale, compute per-protocol:
  - link coverage : fraction of edges with rate > 0 directly
  - pair coverage : fraction of node pairs with a positive-rate path
                    (widest-path / max-min bottleneck over live edges)
  - avg pair rate : SURVIVORSHIP-CORRECTED mean over all pairs, with
                    unreachable pairs counted as zero (do NOT average only
                    over connected pairs).

Two unit views produced from the same underlying rates:
  - bits per channel use  (protocol physics)
  - bits per second       (deployment view, DV*1GHz vs CV*100MHz)

USAGE:
    python Networks/topo_direct.py

Outputs:
    topo_direct_summary.csv      -- one row per (topology, protocol)
    topo_direct_coverage.png     -- link + pair coverage bar chart
    topo_direct_rate.png         -- avg pair rate (both unit views)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import heapq

from topo_loader import TOPOLOGY_FILES, load_topology, edge_lengths
from net_common import link_rate

DV_CLOCK_HZ = 1e9      # commercial-grade DV (Clavis XGR class)
CV_CLOCK_HZ = 100e6    # QOSST-class high-rate CV
CLOCK = {"dv": DV_CLOCK_HZ, "cv_het": CV_CLOCK_HZ}

PROTOS = ("dv", "cv_het")
PROTO_LABEL = {"dv": "DV (decoy BB84)", "cv_het": "CV (GG02 het.)"}
PROTO_COLOR = {"dv": "tab:blue", "cv_het": "tab:red"}


def widest_path_rates(n_nodes, edges, erates):
    """All-pairs bottleneck rates via modified Dijkstra (max-min).
    Returns an (n,n) matrix; 0.0 where no positive-rate path exists."""
    adj = [[] for _ in range(n_nodes)]
    for (u, v), r in zip(edges, erates):
        if r > 0.0:
            adj[u].append((v, r))
            adj[v].append((u, r))
    pair = np.zeros((n_nodes, n_nodes))
    for s in range(n_nodes):
        best = np.zeros(n_nodes)
        best[s] = np.inf
        heap = [(-np.inf, s)]
        while heap:
            neg_w, u = heapq.heappop(heap)
            w = -neg_w
            if w < best[u]:
                continue
            for v, r in adj[u]:
                cand = min(w, r)
                if cand > best[v]:
                    best[v] = cand
                    heapq.heappush(heap, (-cand, v))
        best[s] = 0.0
        pair[s, :] = np.where(np.isfinite(best), best, 0.0)
    return pair


def analyse(name, topo, proto):
    """Direct-link (no relays) metrics for one (topology, protocol)."""
    elens = edge_lengths(topo, "real")
    erates = np.array([link_rate(float(L), proto) for L in elens])
    n = len(topo["node_ids"])
    edges = topo["edges"]

    link_alive = np.count_nonzero(erates > 0.0)
    link_cov = link_alive / len(edges) if edges else 0.0

    pair = widest_path_rates(n, edges, erates)
    iu = np.triu_indices(n, k=1)
    w = pair[iu]
    n_pairs = len(w)
    pair_cov = float(np.count_nonzero(w > 0.0)) / n_pairs if n_pairs else 0.0
    # SURVIVORSHIP-CORRECTED: mean over ALL pairs, failed pairs count as zero
    avg_rate_channel = float(w.sum()) / n_pairs if n_pairs else 0.0
    avg_rate_bps = avg_rate_channel * CLOCK[proto]

    return dict(
        topology=name, protocol=proto, nodes=n, edges=len(edges),
        min_link_km=float(elens.min()), max_link_km=float(elens.max()),
        median_link_km=float(np.median(elens)),
        link_coverage=link_cov, links_alive=int(link_alive),
        pair_coverage=pair_cov,
        avg_pair_rate_channel=avg_rate_channel,
        avg_pair_rate_bps=avg_rate_bps,
    )


def main():
    rows = []
    for name in TOPOLOGY_FILES:
        topo = load_topology(name)
        for proto in PROTOS:
            r = analyse(name, topo, proto)
            rows.append(r)
            print(f"{name:>10} {proto:>7} | "
                  f"link_cov {r['link_coverage']*100:5.1f}% "
                  f"({r['links_alive']}/{r['edges']}) | "
                  f"pair_cov {r['pair_coverage']*100:5.1f}% | "
                  f"rate {r['avg_pair_rate_channel']:.3e} b/ch, "
                  f"{r['avg_pair_rate_bps']:.3e} b/s")

    df = pd.DataFrame(rows)
    df.to_csv("topo_direct_summary.csv", index=False)
    print(f"\nSaved: topo_direct_summary.csv")

    # coverage figure
    names = list(TOPOLOGY_FILES.keys())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    x = np.arange(len(names))
    w = 0.35
    for i, proto in enumerate(PROTOS):
        cov_link = [df[(df.topology == n) & (df.protocol == proto)]
                    ["link_coverage"].iloc[0] * 100 for n in names]
        cov_pair = [df[(df.topology == n) & (df.protocol == proto)]
                    ["pair_coverage"].iloc[0] * 100 for n in names]
        ax1.bar(x + (i - 0.5) * w, cov_link, w,
                color=PROTO_COLOR[proto], label=PROTO_LABEL[proto], alpha=0.85)
        ax2.bar(x + (i - 0.5) * w, cov_pair, w,
                color=PROTO_COLOR[proto], label=PROTO_LABEL[proto], alpha=0.85)
    for ax, ttl in [(ax1, "Direct link coverage"),
                    (ax2, "Pair coverage (widest-path)")]:
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right", fontsize=9)
        ax.set_ylabel("Coverage (%)")
        ax.set_ylim(0, 105)
        ax.set_title(ttl)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis="y")
    fig.suptitle("Direct-link (no-relay) coverage across six topologies "
                 "(real scale)", fontsize=12)
    plt.tight_layout()
    plt.savefig("topo_direct_coverage.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved: topo_direct_coverage.png")

    # rate figure (two unit views)
    fig, (axc, axs) = plt.subplots(1, 2, figsize=(12, 4.5))
    for i, proto in enumerate(PROTOS):
        r_ch = [df[(df.topology == n) & (df.protocol == proto)]
                ["avg_pair_rate_channel"].iloc[0] for n in names]
        r_bs = [df[(df.topology == n) & (df.protocol == proto)]
                ["avg_pair_rate_bps"].iloc[0] for n in names]
        axc.bar(x + (i - 0.5) * w, np.maximum(r_ch, 1e-20), w,
                color=PROTO_COLOR[proto], label=PROTO_LABEL[proto], alpha=0.85)
        axs.bar(x + (i - 0.5) * w, np.maximum(r_bs, 1e-20), w,
                color=PROTO_COLOR[proto], label=PROTO_LABEL[proto], alpha=0.85)
    for ax, ylabel, ttl in [
        (axc, "Avg pair rate (bits/channel use)", "Protocol physics view"),
        (axs, "Avg pair rate (bits/s)", "Deployment view (with clocks)")]:
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right", fontsize=9)
        ax.set_ylabel(ylabel)
        ax.set_title(ttl)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis="y", which="both")
    fig.suptitle("Direct-link (no-relay) average pair rate — survivorship "
                 "corrected", fontsize=12)
    plt.tight_layout()
    plt.savefig("topo_direct_rate.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved: topo_direct_rate.png")


if __name__ == "__main__":
    main()
