"""
topo_direct.py -- Direct-link (no relays) per-topology coverage and key-rate.

For each of the six selected topologies at REAL scale, compute per-protocol:
  - link coverage : fraction of edges with rate > 0 directly
  - pair coverage : fraction of node pairs with a positive-rate path
                    (widest-path / max-min bottleneck over live edges)
  - avg pair rate : SURVIVORSHIP-CORRECTED mean over all pairs, with
                    unreachable pairs counted as zero (do NOT average only
                    over connected pairs).

Rates are reported in BITS/S at a MATCHED 1 GHz clock for both DV and CV. This
figure pins that clock locally (see CLOCK_HZ below) rather than reading
net_common's configurable value, because a matched-clock comparison is what this
figure is for: it isolates the protocol physics, so any DV/CV separation here is
intrinsic to the protocols and not to the electronics. It is NOT the deployed
view — fielded CV symbol rates are ~100 MHz. The CSV retains the raw
bits/channel-use column alongside the bits/s one for traceability.

ONE FIGURE, TWO PANELS SIDE BY SIDE:
  (left)  coverage — link coverage (solid bars) and pair coverage (hatched),
          for both protocols. Link coverage is per-edge; pair coverage is
          per node pair after widest-path routing, so a single dead edge that
          severs the graph drops it far faster than it drops link coverage.
  (right) avg pair rate in bits/s, survivorship corrected.

USAGE:
    python Networks/topo_direct.py

Outputs:
    topo_direct_summary.csv         -- one row per (topology, protocol)
    topo_direct_coverage_rate.png   -- the two-panel figure
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import heapq

from topo_loader import TOPOLOGY_FILES, load_topology, edge_lengths
from net_common import (link_rate, is_reachable, reachable_threshold,
                        REACHABLE_BPS)

# MATCHED CLOCK, PINNED. Both protocols are clocked at 1 GHz for this figure and
# this value is deliberately NOT taken from net_common.clock_hz: net_common's
# clocks are configurable (QKD_CV_CLOCK_HZ) so that other scripts can be switched
# to the deployed view, whereas this figure is defined as the matched-clock
# comparison and must not silently change when that env var is set.
CLOCK_HZ = 1e9

PROTOS = ("dv", "cv_het")
PROTO_LABEL = {"dv": "DV (decoy BB84)", "cv_het": "CV (GG02 het.)"}
PROTO_COLOR = {"dv": "tab:blue", "cv_het": "tab:red"}


def widest_path_rates(n_nodes, edges, erates, floor=0.0):
    """All-pairs bottleneck rates via modified Dijkstra (max-min).
    Edges below `floor` (bits/channel use) are not traversable.
    Returns an (n,n) matrix; 0.0 where no reachable path exists."""
    adj = [[] for _ in range(n_nodes)]
    for (u, v), r in zip(edges, erates):
        if r >= floor:
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

    link_alive = int(np.count_nonzero(is_reachable(erates, proto)))
    link_cov = link_alive / len(edges) if edges else 0.0

    pair = widest_path_rates(n, edges, erates,
                             floor=reachable_threshold(proto))
    iu = np.triu_indices(n, k=1)
    w = pair[iu]
    n_pairs = len(w)
    pair_cov = float(np.count_nonzero(w > 0.0)) / n_pairs if n_pairs else 0.0
    # w is already floored: widest_path_rates dropped sub-floor edges, and a
    # bottleneck is a min over surviving edges, so any positive w clears it.
    # SURVIVORSHIP-CORRECTED: mean over ALL pairs, failed pairs count as zero
    avg_rate_channel = float(w.sum()) / n_pairs if n_pairs else 0.0
    avg_rate_bps = avg_rate_channel * CLOCK_HZ   # same clock for both protocols

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

    # ============================================================
    # ONE FIGURE, TWO PANELS: coverage (left) | key rate (right)
    # No sharey — the panels carry different quantities and different scales
    # (percent on a linear axis vs bits/s on a log axis).
    # ============================================================
    names = list(TOPOLOGY_FILES.keys())
    fig, (axC, axR) = plt.subplots(1, 2, figsize=(14, 5.2))
    x = np.arange(len(names))

    def _col(proto, field):
        return [df[(df.topology == n) & (df.protocol == proto)][field].iloc[0]
                for n in names]

    # ---- left: coverage. Four bars per topology: {DV, CV} x {link, pair}.
    # Link vs pair is encoded by hatching rather than by colour, so protocol
    # stays readable as colour throughout the figure.
    wc = 0.2
    offsets = {("dv", "link"): -1.5, ("dv", "pair"): -0.5,
               ("cv_het", "link"): 0.5, ("cv_het", "pair"): 1.5}
    for proto in PROTOS:
        for kind, field in (("link", "link_coverage"), ("pair", "pair_coverage")):
            axC.bar(x + offsets[(proto, kind)] * wc,
                    np.array(_col(proto, field)) * 100, wc,
                    color=PROTO_COLOR[proto], alpha=0.85,
                    hatch="" if kind == "link" else "///",
                    edgecolor="white", linewidth=0.4)
    axC.set_ylabel("Coverage (%)")
    # Headroom above 100% so the legend sits clear of the bars rather than on
    # top of them; ticks stop at 100 so the axis still reads as a percentage.
    axC.set_ylim(0, 138)
    axC.set_yticks(range(0, 101, 20))
    axC.set_title("Coverage: direct links and routed pairs")
    axC.grid(True, alpha=0.3, axis="y")
    cov_handles = (
        [Patch(facecolor=PROTO_COLOR[p], alpha=0.85, label=PROTO_LABEL[p])
         for p in PROTOS]
        + [Patch(facecolor="0.75", edgecolor="white", label="Link coverage"),
           Patch(facecolor="0.75", edgecolor="white", hatch="///",
                 label="Pair coverage (widest-path)")])
    axC.legend(handles=cov_handles, fontsize=8, ncol=2, loc="upper center",
               framealpha=0.9, borderaxespad=0.4)

    # ---- right: average pair rate, bits/s at the matched 1 GHz clock.
    # Log axis: a zero average cannot be drawn, so it is clamped to the floor
    # below and shows as an absent bar. Absent = zero, not missing data.
    wr = 0.35
    for i, proto in enumerate(PROTOS):
        axR.bar(x + (i - 0.5) * wr, np.maximum(_col(proto, "avg_pair_rate_bps"), 1e-20),
                wr, color=PROTO_COLOR[proto], label=PROTO_LABEL[proto], alpha=0.85)
    axR.set_yscale("log")
    axR.set_ylabel("Average key rate per pair (bits / s)")
    axR.set_title("Key rate: average over all pairs")
    axR.grid(True, alpha=0.3, axis="y", which="both")
    axR.legend(fontsize=9)

    for ax in (axC, axR):
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right", fontsize=9)

    fig.suptitle("Direct-link (no-relay) comparison across six real topologies "
                 "(real scale) — rate is survivorship corrected", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    plt.savefig("topo_direct_coverage_rate.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved: topo_direct_coverage_rate.png")


if __name__ == "__main__":
    main()