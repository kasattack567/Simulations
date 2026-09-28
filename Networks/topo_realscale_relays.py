"""
Real-scale relay analysis — DV vs CV on the 6 real optical backbones.

At true backbone scale (links 100s-1000s of km), neither protocol can carry most
edges directly. Real trusted-node backbones solve
this by placing trusted repeater nodes ALONG each too-long fibre route. We model
this as a MIDPOINT-CHAIN per edge: an edge of length L under a protocol of reach
R is split into ceil(L/R) equal sub-segments by ceil(L/R)-1 evenly-spaced trusted
relays, so every sub-segment is within reach. This matches how deployed
trusted-node backbones (e.g. Beijing-Shanghai) are actually built.

With every edge relayed to viability, BOTH protocols reach 100% edge coverage, so
we compare on two axes:
  - COST: total trusted relays each protocol needs.
  - PERFORMANCE: average end-to-end key rate per pair in BITS/S once relayed,
          routing across the topology at the bottleneck (widest-path) rate. Each
          relayed edge returns at its weakest sub-segment's rate. (This figure
          previously reported bits per channel use with no clock applied at all,
          which was not comparable with the rest of the suite.)

RELAY SIZING (--criterion). Default is 'reach': a hop may be as long as the
distance at which that protocol's key rate falls to the REACHABILITY FLOOR
(net_common.REACHABLE_BPS, 1 kbit/s), so relays are inserted only where a link
would otherwise fail. Note this is NOT the zero-rate distance: sizing hops at
279 km for DV would produce segments carrying a few hundred bit/s, which the
floor then rejects, leaving dead edges in a graph the script claims is fully
covered. Both criteria are therefore derived against a bits/s goal and differ
only in what that goal is.

The alternative, --criterion rate, sizes hops to deliver a service key rate
(--target-bps, default 10 Mbps). It exists because the cost model needs a hop
budget at a usable rate. Be aware of what it does to the comparison: with DV and
CV on the same clock both face the same per-channel-use threshold, and 10 Mbps
sits near the rate at which the two curves cross, so both spans land a few km
apart and the relay counts come out nearly equal — not because the protocols are
equally capable at distance, but because the criterion never looks at the region
where they differ. Use it for costing, not for the reach comparison.

The criterion is applied identically to both protocols either way.

Headline: the relay-count gap (CV needs N x more trusted nodes than DV to light up
the same real backbone), and whether CV's per-link rate advantage survives the
much heavier relaying.

Each trusted relay is also a SECURITY LIABILITY (holds key material in the clear),
so the relay count is a trust-surface cost, not just hardware.

Usage:
    python topo_realscale_relays.py
"""
import argparse
import os
import heapq
import math
import numpy as np
import matplotlib.pyplot as plt

from net_common import (link_rate, span_km, to_bps, is_reachable,
                        REACHABLE_BPS, TARGET_BPS)
from topo_loader import TOPOLOGY_FILES, load_topology, edge_lengths


def edge_relays_and_rate(L_km, protocol, criterion=None, target_bps=None):
    """For an edge of length L, return (n_relays, end_to_end_rate) under the
    midpoint-chain model: split into ceil(L/span) equal sub-segments; the edge
    rate is the bottleneck (equal segments => the single segment rate).

    The hop span comes from span_km(protocol, criterion), so it is protocol- AND
    criterion-dependent. Excess noise is evaluated per hop (link_rate is called on
    seg_len, not on L), which is required under Wang Eq. 12 since xi_r ~ 1/(eta*T).
    """
    reach = span_km(protocol, criterion, target_bps)
    if L_km <= 0 or reach <= 0:
        return 0, 0.0
    n_seg = max(1, math.ceil(L_km / reach))
    seg_len = L_km / n_seg
    seg_rate = link_rate(seg_len, protocol)
    return (n_seg - 1), seg_rate   # n_relays, bottleneck rate


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


def analyse(topo, protocol, criterion=None, target_bps=None):
    """Return dict: total_relays, avg_rate in BITS/S (all pairs), coverage, over
    the real topology with midpoint-chain relays on every edge.

    Bottleneck routing is done per channel use (a min over hops is invariant
    under a common positive scaling, so the widest path is the same either way);
    the single clock conversion is applied once to the final average."""
    n = len(topo["node_ids"])
    edges = topo["edges"]
    elen = edge_lengths(topo, "real")
    adj = [[] for _ in range(n)]
    total_relays = 0
    for (a, b), L in zip(edges, elen):
        nr, rate = edge_relays_and_rate(L, protocol, criterion, target_bps)
        total_relays += nr
        if is_reachable(rate, protocol):
            adj[a].append((b, rate)); adj[b].append((a, rate))
    total, connected = 0.0, 0
    n_pairs = n * (n - 1) // 2
    for s in range(n):
        caps = widest_paths_from(s, n, adj)
        for t in range(s + 1, n):
            if caps[t] not in (0.0, float("inf")):
                total += caps[t]; connected += 1
    return dict(total_relays=total_relays,
                avg_rate=to_bps(total / n_pairs, protocol) if n_pairs else 0.0,
                coverage=connected / n_pairs if n_pairs else 0.0,
                n=n, n_edges=len(edges))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cv", choices=["heterodyne", "homodyne"], default="heterodyne")
    p.add_argument("--criterion", choices=["reach", "rate"], default="reach",
                   help="hop sizing: 'reach' = hops sized at the reachability "
                        "floor, 1 kbit/s (default); 'rate' = hops sized to hold "
                        "--target-bps")
    p.add_argument("--target-bps", type=float, default=TARGET_BPS,
                   help="service key rate, used only with --criterion rate "
                        "(default 10 Mbps)")
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()
    cv_proto = "cv_hom" if args.cv == "homodyne" else "cv_het"
    names = list(TOPOLOGY_FILES.keys())

    res = {}
    spans = {q: span_km(q, args.criterion, args.target_bps)
             for q in ("dv", cv_proto)}
    print(f"Real-scale midpoint-chain relays. CV={args.cv}. "
          f"criterion={args.criterion}"
          + (f" @ {args.target_bps/1e6:g} Mbps" if args.criterion == "rate" else ""))
    print(f"  hop span: DV {spans['dv']:.1f} km, CV {spans[cv_proto]:.1f} km\n")
    print(f"{'Topology':12} {'nodes':>5} {'edges':>5} "
          f"{'DV relays':>9} {'CV relays':>9} {'ratio':>6} "
          f"{'DV bit/s':>10} {'CV bit/s':>10}")
    zero_relay = []
    for name in names:
        topo = load_topology(name)
        dv = analyse(topo, "dv", args.criterion, args.target_bps)
        cv = analyse(topo, cv_proto, args.criterion, args.target_bps)
        res[name] = (dv, cv)
        ratio = (cv["total_relays"] / dv["total_relays"]
                 if dv["total_relays"] else float("inf"))
        if dv["total_relays"] == 0 or cv["total_relays"] == 0:
            zero_relay.append(name)
        print(f"{name:12} {dv['n']:5d} {dv['n_edges']:5d} "
              f"{dv['total_relays']:9d} {cv['total_relays']:9d} "
              f"{ratio:6.1f} {dv['avg_rate']:10.2e} {cv['avg_rate']:10.2e}")

    if zero_relay:
        print(f"\n  [note] a protocol needs zero relays on: {', '.join(zero_relay)}"
              "\n         the CV/DV ratio is undefined there (division by zero).")

    # ============================================================
    # PLOT: relay count (left) + avg rate (right), DV vs CV per topology
    # ============================================================
    plt.rcParams.update({"font.family": "serif", "font.size": 11})
    fig, (axK, axR) = plt.subplots(1, 2, figsize=(16, 6))
    x = np.arange(len(names)); w = 0.38
    DV_C, CV_C = "#1f4e9c", "#c0392b"

    dv_relays = [res[n][0]["total_relays"] for n in names]
    cv_relays = [res[n][1]["total_relays"] for n in names]
    axK.bar(x - w/2, dv_relays, w, color=DV_C, label="DV — decoy BB84")
    axK.bar(x + w/2, cv_relays, w, color=CV_C, label=f"CV — {args.cv}")
    axK.set_xticks(x); axK.set_xticklabels(names, rotation=30, ha="right", fontsize=9)
    axK.set_ylabel("Trusted relays needed for full coverage")
    axK.set_yscale("log")
    axK.set_title("Cost: trusted relays to light up the real backbone", fontsize=12)
    axK.grid(True, axis="y", which="both", alpha=0.3)
    axK.legend(fontsize=9)
    for i, (d, c) in enumerate(zip(dv_relays, cv_relays)):
        if d > 0: axK.text(i - w/2, d, str(d), ha="center", va="bottom", fontsize=8)
        if c > 0: axK.text(i + w/2, c, str(c), ha="center", va="bottom", fontsize=8)

    dv_rate = [res[n][0]["avg_rate"] for n in names]
    cv_rate = [res[n][1]["avg_rate"] for n in names]
    axR.bar(x - w/2, dv_rate, w, color=DV_C, label="DV — decoy BB84")
    axR.bar(x + w/2, cv_rate, w, color=CV_C, label=f"CV — {args.cv}")
    axR.set_xticks(x); axR.set_xticklabels(names, rotation=30, ha="right", fontsize=9)
    axR.set_ylabel("Average key rate per pair, relayed (bits / s)")
    axR.set_yscale("log")
    axR.set_title("Performance: average key rate once relayed to full coverage",
                  fontsize=12)
    axR.grid(True, axis="y", which="both", alpha=0.3)
    axR.legend(fontsize=9)

    fig.suptitle("Real-scale trusted-relay comparison: DV vs CV on 6 real "
                 "optical backbones", fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or f"topo_realscale_relays_{args.cv}.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()