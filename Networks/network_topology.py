"""
Topology comparison — DV vs CV across canonical network topologies.

For each topology (ring, star, tree, mesh), sweep the network SIZE (area side
length) and plot the AVERAGE key rate per user pair vs size, comparing DV and CV.

Why average (not total): different topologies have very different link counts
(mesh N(N-1)/2, ring N, star/tree N-1), so total rate would mostly measure "how
many links" rather than per-connection quality. Average per pair (over ALL pairs,
unreachable = 0) normalises this and keeps the topologies comparable, while still
penalising a topology that cannot connect distant pairs.

Routing model: only a topology's EDGES are physical QKD links. A key between two
non-adjacent users is relayed along the topology through trusted intermediate
nodes; its end-to-end rate is the BOTTLENECK (min) link rate along the path
(widest-path / max-min routing). Mesh = every pair is a direct edge.

Node placement per topology (in an area x area square):
  - ring : nodes evenly spaced on a circle inscribed in the area
  - star : one central hub + spokes to nodes on a surrounding circle
  - tree : hub + two-level branching (hub -> branch nodes -> leaves)
  - mesh : nodes uniformly at random (all pairs direct)

Usage:
    python network_topology.py
    python network_topology.py --n 12 --runs 5 --cv homodyne
"""
import argparse
import os
import heapq
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

from net_common import (link_rate, to_bps, is_reachable, REACHABLE_BPS,
                        CV_REACH_KM, DV_REACH_KM)

# NOTE ON EDGE VIABILITY. This script places NO relays: it measures how far each
# protocol's *unaided* reach carries a topology. There is therefore no artificial
# reach cutoff — link_rate() already returns exactly 0.0 beyond a protocol's
# reach, so the physics decides. (Gating edges with a relay-SPACING span, e.g.
# net_common.span_km(..., 'rate'), would be wrong here: that quantity answers
# "how far apart may relays sit", not "does this edge carry a key at all", and
# would kill DV edges beyond 65 km that in fact work out to 279 km.)
#
# For a service-provisioning view instead of a connectivity view, use --min-rate:
# an edge is then viable only if it meets that rate, and "coverage" reads as
# "% of pairs meeting the service target" rather than "% of pairs with any key".
REACH = {"dv": DV_REACH_KM, "cv_het": CV_REACH_KM, "cv_hom": CV_REACH_KM}


# ============================================================
# TOPOLOGY CONSTRUCTION — returns (positions, edge_list)
# ============================================================
def make_ring(n, area):
    r = area / 2.0
    cx = cy = area / 2.0
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pos = np.column_stack([cx + r * np.cos(ang), cy + r * np.sin(ang)])
    edges = [(i, (i + 1) % n) for i in range(n)]     # ring edges
    return pos, edges


def make_star(n, area):
    # hub at centre, n-1 users on a surrounding circle
    r = area / 2.0
    cx = cy = area / 2.0
    hub = np.array([[cx, cy]])
    ang = np.linspace(0, 2 * np.pi, n - 1, endpoint=False)
    spokes = np.column_stack([cx + r * np.cos(ang), cy + r * np.sin(ang)])
    pos = np.vstack([hub, spokes])
    edges = [(0, k) for k in range(1, n)]            # hub-to-spoke
    return pos, edges


def make_tree(n, area):
    # simple 2-level tree: hub (0) -> branches -> leaves, laid out radially
    cx = cy = area / 2.0
    n_branch = max(2, int(round(np.sqrt(n - 1))))
    pos = [[cx, cy]]
    edges = []
    remaining = n - 1
    ang_b = np.linspace(0, 2 * np.pi, n_branch, endpoint=False)
    r1 = area / 4.0
    r2 = area / 2.0
    branch_ids = []
    for a in ang_b:
        if remaining <= 0:
            break
        bid = len(pos)
        pos.append([cx + r1 * np.cos(a), cy + r1 * np.sin(a)])
        edges.append((0, bid)); branch_ids.append((bid, a)); remaining -= 1
    # distribute leaves among branches
    li = 0
    while remaining > 0 and branch_ids:
        bid, a = branch_ids[li % len(branch_ids)]
        jitter = 0.4 * (li // len(branch_ids) + 1) * (-1) ** li
        lid = len(pos)
        pos.append([cx + r2 * np.cos(a + 0.2 * jitter),
                    cy + r2 * np.sin(a + 0.2 * jitter)])
        edges.append((bid, lid)); remaining -= 1; li += 1
    return np.array(pos, float), edges


def make_mesh(n, area, seed):
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, area, size=(n, 2))
    edges = list(combinations(range(n), 2))          # all pairs direct
    return pos, edges


# ============================================================
# ROUTING — average key rate per pair over the topology
# ============================================================
def edge_rate(pos, a, b, protocol, min_rate=0.0):
    """Link rate for edge (a,b) under `protocol`, in BITS/S.

    No artificial reach cutoff: link_rate() is already exactly 0.0 beyond a
    protocol's reach. `min_rate` is a service threshold in BITS/S (it used to be
    in bits/channel use — the units changed with the switch to a bits/s output,
    so any --min-rate value from an earlier run must be rescaled by the clock).
    A positive min_rate turns coverage into "% of pairs meeting the target"."""
    d = float(np.hypot(*(pos[a] - pos[b])))
    r_ch = link_rate(d, protocol)
    if not is_reachable(r_ch, protocol):
        return 0.0
    r = to_bps(r_ch, protocol)
    return r if r >= min_rate else 0.0


def widest_paths_from(src, n, adj):
    """Max-min (widest path) capacities from src to all nodes."""
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


def avg_rate_per_pair(pos, edges, protocol, min_rate=0.0):
    """Average end-to-end key rate over ALL user pairs (unreachable = 0),
    routing along the topology edges with bottleneck (widest-path) rates."""
    n = len(pos)
    adj = [[] for _ in range(n)]
    for a, b in edges:
        r = edge_rate(pos, a, b, protocol, min_rate)
        if r > 0.0:      # edge_rate already returns 0.0 below the floor
            adj[a].append((b, r)); adj[b].append((a, r))
    total = 0.0
    n_pairs = n * (n - 1) // 2
    for s in range(n):
        caps = widest_paths_from(s, n, adj)
        for t in range(s + 1, n):
            if caps[t] not in (0.0, float("inf")):
                total += caps[t]
    return total / n_pairs if n_pairs else 0.0


def rate_and_coverage(pos, edges, protocol, min_rate=0.0):
    """Return (avg_rate_over_all_pairs, coverage_fraction).
    A pair is 'covered' if a routed path connects it (bottleneck rate > 0)."""
    n = len(pos)
    adj = [[] for _ in range(n)]
    for a, b in edges:
        r = edge_rate(pos, a, b, protocol, min_rate)
        if r > 0.0:      # edge_rate already returns 0.0 below the floor
            adj[a].append((b, r)); adj[b].append((a, r))
    total = 0.0
    connected = 0
    n_pairs = n * (n - 1) // 2
    for s in range(n):
        caps = widest_paths_from(s, n, adj)
        for t in range(s + 1, n):
            if caps[t] not in (0.0, float("inf")):
                total += caps[t]
                connected += 1
    avg = total / n_pairs if n_pairs else 0.0
    cov = connected / n_pairs if n_pairs else 0.0
    return avg, cov


# ============================================================
# SWEEP + PLOT
# ============================================================
TOPOS = [("Ring", make_ring), ("Star", make_star),
         ("Tree", make_tree), ("Mesh", make_mesh)]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=12, help="Nodes per topology")
    p.add_argument("--runs", type=int, default=5,
                   help="Layouts averaged (mesh is random; ring/star/tree fixed)")
    p.add_argument("--areas", type=float, nargs="+",
                   default=[10, 50, 100, 150, 175, 200, 250, 300, 350, 400])
    p.add_argument("--cv", choices=["heterodyne", "homodyne"], default="heterodyne")
    p.add_argument("--min-rate", type=float, default=0.0,
                   help="EXTRA service threshold in bits/s, on top of the "
                        "net_common reachability floor. 0 (default) = the floor "
                        "alone.")
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    cv_proto = "cv_hom" if args.cv == "homodyne" else "cv_het"
    areas = np.array(args.areas, dtype=float)

    plt.rcParams.update({"font.family": "serif", "font.size": 10})
    # 2 rows (rate, coverage) x 4 cols (topologies)
    fig, axes = plt.subplots(2, 4, figsize=(20, 9), sharex=True)
    # dotted vertical + "-> 0" in the top row marks where a protocol's average
    # rate becomes exactly zero (no connected pairs), which a log axis cannot show.

    view = (f"connectivity (>= {REACHABLE_BPS:g} bits/s floor)"
            if args.min_rate <= 0
            else f"service (>= {args.min_rate:g} bits/s)")
    print(f"N={args.n} nodes, CV={args.cv}. Reach: DV {DV_REACH_KM:.0f} km, "
          f"CV {CV_REACH_KM:.0f} km. View: {view}\n")
    for col, (tname, builder) in enumerate(TOPOS):
        ax_rate = axes[0, col]
        ax_cov = axes[1, col]
        dv_r, cv_r, dv_c, cv_c = [], [], [], []
        for area in areas:
            dvr, cvr, dvc, cvc = [], [], [], []
            runs = args.runs if tname == "Mesh" else 1
            for r in range(runs):
                if tname == "Mesh":
                    pos, edges = builder(args.n, area, seed=r)
                else:
                    pos, edges = builder(args.n, area)
                a_dv, c_dv = rate_and_coverage(pos, edges, "dv", args.min_rate)
                a_cv, c_cv = rate_and_coverage(pos, edges, cv_proto, args.min_rate)
                dvr.append(a_dv); cvr.append(a_cv); dvc.append(c_dv); cvc.append(c_cv)
            dv_r.append(np.mean(dvr)); cv_r.append(np.mean(cvr))
            dv_c.append(np.mean(dvc)); cv_c.append(np.mean(cvc))
        dv_r, cv_r = np.array(dv_r), np.array(cv_r)
        dv_c, cv_c = np.array(dv_c), np.array(cv_c)

        # --- rate panel (top) ---
        # A rate of exactly 0 cannot be drawn on a log axis, so a curve that
        # reaches zero simply STOPS. That reads as missing data when it actually
        # means "no connected pairs". Mark the extinction point explicitly.
        def _plot_rate(ax, x, y, color, **kw):
            y = np.array(y, float)
            mask = y > 0
            ax.plot(x[mask], y[mask], color=color, **kw)
            if mask.any() and not mask.all():
                i = int(np.where(mask)[0][-1])       # last surviving point
                if i + 1 < len(x):
                    ax.axvline(x[i + 1], color=color, ls=":", lw=1.4, alpha=0.8,
                               zorder=1)
                    ax.annotate("→ 0", xy=(x[i + 1], y[mask][-1]),
                                xytext=(4, 0), textcoords="offset points",
                                color=color, fontsize=8, va="center")
        _plot_rate(ax_rate, areas, dv_r, "#1f4e9c", marker="o", ls="-",
                   label="DV — decoy BB84", ms=5)
        _plot_rate(ax_rate, areas, cv_r, "#c0392b", marker="s", ls="-",
                   label=f"CV — {args.cv}", ms=5)
        ax_rate.axvspan(10, 40, color="#7cc47f", alpha=0.12, zorder=0,
                        label="Metro-scale span")
        ax_rate.set_yscale("log")
        ax_rate.set_title(f"{tname}", fontsize=13, weight="bold")
        ax_rate.grid(True, which="both", alpha=0.3)
        if col == 0:
            ax_rate.set_ylabel("Average key rate per pair\n(bits / s)")
        ax_rate.legend(fontsize=8, loc="lower left")

        # --- coverage panel (bottom) ---
        # DV drawn dashed + thicker so it stays visible where it overlaps CV at 100%
        ax_cov.plot(areas, dv_c * 100, "--o", color="#1f4e9c", label="DV",
                    ms=5, lw=2.5, zorder=3)
        ax_cov.plot(areas, cv_c * 100, "-s", color="#c0392b", label="CV",
                    ms=5, lw=1.5, zorder=2)
        ax_cov.axvspan(10, 40, color="#7cc47f", alpha=0.12, zorder=0)
        ax_cov.set_ylim(0, 105)
        ax_cov.grid(True, alpha=0.3)
        ax_cov.set_xlabel("Area side (km)")
        if col == 0:
            ax_cov.set_ylabel("Coverage\n(% of pairs connected)")
        ax_cov.legend(fontsize=8, loc="center right")

        print(f"{tname:5s}: CV rate {cv_r[0]:.2e}->{cv_r[-1]:.2e}  "
              f"cov {cv_c[0]*100:.0f}%->{cv_c[-1]*100:.0f}%   "
              f"DV rate {dv_r[0]:.2e}->{dv_r[-1]:.2e}  "
              f"cov {dv_c[0]*100:.0f}%->{dv_c[-1]*100:.0f}%")

    fig.suptitle(f"Topology comparison: DV vs CV — rate (top) and coverage (bottom)   "
                 f"(N={args.n} nodes, routed along topology)",
                 fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or f"network_topology_all4_N{args.n}_{args.cv}.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()