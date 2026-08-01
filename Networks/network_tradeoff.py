"""
Cost-performance tradeoff — DV vs CV once relays give BOTH full coverage.

Builds on the min-K result: at each area size, place the MINIMUM number of
trusted relays that bring CV to 100% coverage (k-means placement, dedicated
relays, users are not transit hops). With both protocols now at full coverage,
coverage is no longer the differentiator, so we compare on RATE, with the relay
count as the COST.

Crucially, a relayed CV pair does NOT keep its direct rate: the end-to-end rate
is the BOTTLENECK (widest-path max-min) over the hops along its route, and
relayed pairs are exactly the long ones whose hops sit near CV's reach limit, so
they return at LOW rate. This is computed honestly here (max-min routing), not
assumed.

Two panels vs area:
  (left)  total network key rate in BITS/S: DV (all direct) vs CV (direct + relayed)
  (right) minimum trusted relays CV needs (DV needs zero)

Together: "CV delivers X total rate for N trusted nodes; DV delivers Y for none."

Usage:
    python network_tradeoff.py
    python network_tradeoff.py --n 20 --runs 10 --cv homodyne
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import heapq

from net_common import (place_users, all_pairs, pair_distance_km, link_rate,
                        CV_REACH_KM, span_km, to_bps, is_reachable)


# ---------------- k-means (self-contained) ----------------
def kmeans(points, k, seed=0, iters=50):
    rng = np.random.default_rng(seed)
    if k <= 0:
        return np.empty((0, 2))
    if k >= len(points):
        return points.copy()
    centres = [points[rng.integers(len(points))]]
    for _ in range(1, k):
        d2 = np.min([np.sum((points - c) ** 2, axis=1) for c in centres], axis=0)
        probs = d2 / d2.sum() if d2.sum() > 0 else None
        idx = rng.choice(len(points), p=probs)
        centres.append(points[idx])
    centres = np.array(centres)
    for _ in range(iters):
        dists = np.array([np.sum((points - c) ** 2, axis=1) for c in centres])
        labels = np.argmin(dists, axis=0)
        new = np.array([points[labels == m].mean(axis=0) if np.any(labels == m)
                        else centres[m] for m in range(k)])
        if np.allclose(new, centres):
            break
        centres = new
    return centres


# ---------------- coverage check (dedicated relays; users not transit) --------
def full_coverage(users, relays, reach_km):
    n_u = len(users); n_r = len(relays)
    if n_r:
        radj = [[] for _ in range(n_r)]
        for a in range(n_r):
            for b in range(a + 1, n_r):
                if np.hypot(*(relays[a] - relays[b])) <= reach_km:
                    radj[a].append(b); radj[b].append(a)
        rcomp = [-1] * n_r; c = 0
        for s in range(n_r):
            if rcomp[s] == -1:
                st = [s]; rcomp[s] = c
                while st:
                    x = st.pop()
                    for y in radj[x]:
                        if rcomp[y] == -1:
                            rcomp[y] = c; st.append(y)
                c += 1
        user_rcomps = [{rcomp[k] for k in range(n_r)
                        if np.hypot(*(users[i] - relays[k])) <= reach_km}
                       for i in range(n_u)]
    else:
        user_rcomps = [set() for _ in range(n_u)]
    for i, j in combinations(range(n_u), 2):
        if np.hypot(*(users[i] - users[j])) <= reach_km:
            continue
        if len(user_rcomps[i] & user_rcomps[j]) == 0:
            return False
    return True


def min_relays_positions(users, reach_km, kmax, seed):
    """Return (K, relay_positions) for the minimum K giving full coverage."""
    if full_coverage(users, np.empty((0, 2)), reach_km):
        return 0, np.empty((0, 2))
    for k in range(1, kmax + 1):
        relays = kmeans(users, k, seed=seed)
        if full_coverage(users, relays, reach_km):
            return k, relays
    return kmax + 1, kmeans(users, kmax, seed=seed)


# ---------------- widest-path (max-min) routing for end-to-end rate ----------
def relayed_pair_rate(i, j, users, relays, protocol, reach_km):
    """End-to-end rate for user pair (i,j): direct if reachable, else the
    max-min (widest-path) bottleneck rate through the users+relays graph where
    ONLY relays are transit nodes (users i, j are endpoints only)."""
    d_ij = np.hypot(*(users[i] - users[j]))
    direct = link_rate(d_ij, protocol) if d_ij <= reach_km else 0.0
    if len(relays) == 0:
        return direct

    # nodes: 0..n_u-1 users (but only i and j usable as endpoints), then relays
    n_u = len(users)
    relay_ids = list(range(n_u, n_u + len(relays)))
    coords = np.vstack([users, relays])

    # build feasible edges (<= reach): i<->relays, j<->relays, relay<->relay,
    # and the direct i<->j. Users other than i,j are NOT nodes in the path graph.
    usable = {i, j, *relay_ids}
    adj = {u: [] for u in usable}
    def maybe_edge(a, b):
        d = np.hypot(*(coords[a] - coords[b]))
        if d <= reach_km:
            r = link_rate(d, protocol)
            if is_reachable(r, protocol):
                adj[a].append((b, r)); adj[b].append((a, r))
    # i,j to relays
    for rid in relay_ids:
        maybe_edge(i, rid); maybe_edge(j, rid)
    # relay-relay
    for a_idx in range(len(relay_ids)):
        for b_idx in range(a_idx + 1, len(relay_ids)):
            maybe_edge(relay_ids[a_idx], relay_ids[b_idx])
    # direct i-j
    if direct > 0:
        adj[i].append((j, direct)); adj[j].append((i, direct))

    # widest path from i to j: maximise the minimum edge along the path
    best = {u: 0.0 for u in usable}
    best[i] = float("inf")
    pq = [(-float("inf"), i)]
    while pq:
        negcap, u = heapq.heappop(pq)
        cap = -negcap
        if cap < best[u]:
            continue
        if u == j:
            return best[j] if best[j] != float("inf") else direct
        for v, w in adj[u]:
            nc = min(cap, w)
            if nc > best[v]:
                best[v] = nc
                heapq.heappush(pq, (-nc, v))
    return best[j] if best[j] not in (0.0, float("inf")) else direct


def network_total_relayed(users, relays, protocol, reach_km):
    """Total end-to-end key rate in BITS/S over all user pairs, with relay
    routing. Summed per channel use, converted once at the protocol's clock."""
    total = 0.0
    for i, j in combinations(range(len(users)), 2):
        total += relayed_pair_rate(i, j, users, relays, protocol, reach_km)
    return to_bps(total, protocol)


def dv_total_direct(users, reach_km):
    """DV total in BITS/S over all pairs (direct; DV reaches all pairs)."""
    total = 0.0
    for i, j in all_pairs(len(users)):
        total += link_rate(pair_distance_km(users, i, j), "dv")
    return to_bps(total, "dv")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=20)
    p.add_argument("--runs", type=int, default=10)
    p.add_argument("--areas", type=float, nargs="+",
                   default=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    p.add_argument("--kmax", type=int, default=25)
    p.add_argument("--cv", choices=["heterodyne", "homodyne"], default="heterodyne")
    p.add_argument("--reach", type=float, default=span_km("cv_het"),
                   help="CV hop span, km (default: net_common sizing criterion)")
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    cv_proto = "cv_hom" if args.cv == "homodyne" else "cv_het"
    areas = np.array(args.areas, dtype=float)
    dv_tot = np.zeros(len(areas)); dv_std = np.zeros(len(areas))
    cv_tot = np.zeros(len(areas)); cv_std = np.zeros(len(areas))
    k_mean = np.zeros(len(areas)); k_std = np.zeros(len(areas))

    print(f"N={args.n}, {args.runs} layouts/area, CV={args.cv}, reach {args.reach:.0f} km\n")
    for a_idx, area in enumerate(areas):
        dvs, cvs, ks = [], [], []
        for r in range(args.runs):
            users = place_users(args.n, area_km=area, seed=9000 * a_idx + r)
            k, relays = min_relays_positions(users, args.reach, args.kmax,
                                             seed=9000 * a_idx + r)
            dvs.append(dv_total_direct(users, args.reach))
            cvs.append(network_total_relayed(users, relays, cv_proto, args.reach))
            ks.append(k)
        dv_tot[a_idx], dv_std[a_idx] = np.mean(dvs), np.std(dvs)
        cv_tot[a_idx], cv_std[a_idx] = np.mean(cvs), np.std(cvs)
        k_mean[a_idx], k_std[a_idx] = np.mean(ks), np.std(ks)
        print(f"area {area:4.0f} km:  DV total {np.mean(dvs):.3e}   "
              f"CV total {np.mean(cvs):.3e}   relays {np.mean(ks):.1f}   "
              f"CV/DV {np.mean(cvs)/np.mean(dvs):.2f}")

    # ============================================================
    plt.rcParams.update({"font.family": "serif", "font.size": 12})
    fig, (axR, axK) = plt.subplots(1, 2, figsize=(15, 6))

    # Rates are already bits/s (converted in the total_* functions above).
    axR.plot(areas, dv_tot, "-o", color="#1f4e9c",
             label="DV — decoy BB84", markersize=5)
    axR.fill_between(areas, np.maximum(dv_tot - dv_std, 1e-12),
                     dv_tot + dv_std, color="#1f4e9c", alpha=0.15)
    axR.plot(areas, cv_tot, "-s", color="#c0392b",
             label=f"CV — GG02 {args.cv} (relayed)", markersize=5)
    axR.fill_between(areas, np.maximum(cv_tot - cv_std, 1e-12),
                     cv_tot + cv_std, color="#c0392b", alpha=0.15)
    axR.axvspan(10, 40, color="#7cc47f", alpha=0.12, zorder=0,
                label="Deployed area band")
    axR.set_xlabel("Area side length (km)")
    axR.set_ylabel("Total network key rate (bits / s)")
    axR.set_title("Performance: total key rate (both at 100% coverage)")
    axR.set_yscale("log")
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=10)

    axK.plot(areas, np.zeros_like(areas), "-o", color="#1f4e9c",
             label="DV (needs none)", markersize=5)
    axK.plot(areas, k_mean, "-s", color="#c0392b", label="CV trusted relays",
             markersize=5)
    axK.fill_between(areas, np.maximum(k_mean - k_std, 0), k_mean + k_std,
                     color="#c0392b", alpha=0.15)
    axK.axvspan(10, 40, color="#7cc47f", alpha=0.12, zorder=0)
    axK.set_xlabel("Area side length (km)")
    axK.set_ylabel("Minimum trusted relays for full coverage")
    axK.set_title("Cost: trusted nodes CV requires")
    axK.grid(True, alpha=0.3)
    axK.legend(fontsize=10)
    axK.set_ylim(bottom=-0.5)

    fig.suptitle(f"Cost vs performance: DV vs CV at full coverage   "
                 f"(N={args.n} users, {args.runs} layouts/area)",
                 fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or f"network_tradeoff_N{args.n}_{args.cv}_bps.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()