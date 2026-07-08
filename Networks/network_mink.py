"""
Network analysis — minimum trusted nodes for full CV coverage vs metro area.

THE HEADLINE COST CURVE.

CV's ~40 km reach means that as a metro area grows, more and more user pairs
exceed a direct CV link and must be routed through SHARED trusted relays. This
script finds, for each area size, the MINIMUM number of trusted relays K needed
so that EVERY user pair is connected (directly or via a relay path where every
hop is under CV reach). DV needs zero relays anywhere in the metro range
(reach ~275 km), so the contrast is the cost of choosing CV:

    "CV needs N trusted nodes to cover a metro of this size; DV needs none."

Each trusted node is both a hardware cost and a SECURITY LIABILITY (it holds key
material in the clear), so N is a meaningful deployment cost.

Model:
  - Users placed uniformly at random; K relays placed by k-means on user
    positions (relays at cluster centres — standard, matches the peer's method).
  - Build a graph on users+relays; connect any two nodes whose separation is
    under CV reach (a feasible QKD hop). A user pair is "covered" if a path
    exists (direct or multi-hop through relays); end-to-end security is the
    trusted-node chain.
  - Minimum K = smallest relay count giving 100% pair coverage.
  - Averaged over random layouts; plotted vs area with a +/- 1 std band.

Usage:
    python network_mink.py
    python network_mink.py --n 20 --runs 10 --kmax 15
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

from net_common import place_users, CV_REACH_KM


# ------------------------------------------------------------------
# lightweight k-means (no sklearn dependency)
# ------------------------------------------------------------------
def kmeans(points, k, seed=0, iters=50):
    """Return k cluster-centre positions for `points` (n,2)."""
    rng = np.random.default_rng(seed)
    if k >= len(points):
        return points.copy()
    # k-means++ style init: pick first centre at random, then spread out
    centres = [points[rng.integers(len(points))]]
    for _ in range(1, k):
        d2 = np.min([np.sum((points - c) ** 2, axis=1) for c in centres], axis=0)
        probs = d2 / d2.sum() if d2.sum() > 0 else None
        idx = rng.choice(len(points), p=probs)
        centres.append(points[idx])
    centres = np.array(centres)
    for _ in range(iters):
        # assign
        dists = np.array([np.sum((points - c) ** 2, axis=1) for c in centres])
        labels = np.argmin(dists, axis=0)
        # update
        new = np.array([points[labels == m].mean(axis=0) if np.any(labels == m)
                        else centres[m] for m in range(k)])
        if np.allclose(new, centres):
            break
        centres = new
    return centres


# ------------------------------------------------------------------
# coverage via relay reachability
# ------------------------------------------------------------------
def full_coverage(users, relays, reach_km):
    """True if EVERY user pair is connectable, under the dedicated-relay model:
    a pair (u,v) is connected iff either
      (a) they are within reach directly (a direct QKD link), OR
      (b) both can reach the relay backbone and the relays they reach are in the
          same relay-connected component.
    USER NODES ARE NOT TRANSIT HOPS — a user never relays for a third party
    (this preserves the assumption that users need not trust one another; only
    the operator's dedicated relays are trusted intermediaries).
    """
    n_u = len(users)
    n_r = len(relays)

    # relay-relay adjacency + relay connected components
    if n_r:
        radj = [[] for _ in range(n_r)]
        for a in range(n_r):
            for b in range(a + 1, n_r):
                if np.hypot(*(relays[a] - relays[b])) <= reach_km:
                    radj[a].append(b); radj[b].append(a)
        rcomp = [-1] * n_r; c = 0
        for s in range(n_r):
            if rcomp[s] == -1:
                stack = [s]; rcomp[s] = c
                while stack:
                    x = stack.pop()
                    for y in radj[x]:
                        if rcomp[y] == -1:
                            rcomp[y] = c; stack.append(y)
                c += 1
        # which relay components each user can reach (user within reach of a relay)
        user_rcomps = []
        for i in range(n_u):
            comps = {rcomp[k] for k in range(n_r)
                     if np.hypot(*(users[i] - relays[k])) <= reach_km}
            user_rcomps.append(comps)
    else:
        user_rcomps = [set() for _ in range(n_u)]

    # check every user pair
    for i, j in combinations(range(n_u), 2):
        direct = np.hypot(*(users[i] - users[j])) <= reach_km
        via_relay = len(user_rcomps[i] & user_rcomps[j]) > 0
        if not (direct or via_relay):
            return False
    return True


def min_relays(users, reach_km, kmax, seed):
    """Smallest K (0..kmax) of k-means relays giving full user coverage."""
    if full_coverage(users, np.empty((0, 2)), reach_km):
        return 0
    for k in range(1, kmax + 1):
        relays = kmeans(users, k, seed=seed)
        if full_coverage(users, relays, reach_km):
            return k
    return kmax + 1   # not achieved within kmax (flag)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=20, help="Number of users")
    p.add_argument("--runs", type=int, default=10, help="Random layouts per area")
    p.add_argument("--areas", type=float, nargs="+",
                   default=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100],)
    p.add_argument("--kmax", type=int, default=20, help="Max relays to try")
    p.add_argument("--reach", type=float, default=CV_REACH_KM,
                   help="CV reach (km); hops must be under this")
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    areas = np.array(args.areas, dtype=float)
    mean_k = np.zeros(len(areas))
    std_k = np.zeros(len(areas))

    print(f"N={args.n} users, {args.runs} layouts/area, CV reach {args.reach:.0f} km, "
          f"k-means relays\n")
    for a_idx, area in enumerate(areas):
        ks = []
        for r in range(args.runs):
            users = place_users(args.n, area_km=area, seed=7000 * a_idx + r)
            k = min_relays(users, args.reach, args.kmax, seed=7000 * a_idx + r)
            ks.append(k)
        mean_k[a_idx] = np.mean(ks)
        std_k[a_idx] = np.std(ks)
        print(f"area {area:4.0f} km:  min relays = {np.mean(ks):4.1f} "
              f"(std {np.std(ks):.1f}, range {min(ks)}-{max(ks)})")

    # ============================================================
    # PLOT — the cost curve
    # ============================================================
    plt.rcParams.update({"font.family": "serif", "font.size": 12})
    fig, ax = plt.subplots(figsize=(10, 6.2))

    # DV needs zero relays across the whole metro range
    ax.plot(areas, np.zeros_like(areas), "-o", color="#1f4e9c",
            label="DV — decoy BB84 (reach ~275 km)", markersize=5)
    # CV min relays
    ax.plot(areas, mean_k, "-s", color="#c0392b",
            label=f"CV — GG02 (reach ~{args.reach:.0f} km)", markersize=5)
    ax.fill_between(areas, np.maximum(mean_k - std_k, 0), mean_k + std_k,
                    color="#c0392b", alpha=0.15)

    # shade the realistic deployed-metro band (Madrid 1.9-33 km, Tokyo ~45 km)
    ax.axvspan(10, 40, color="#7cc47f", alpha=0.12, zorder=0,
               label="Deployed metro band (~10-40 km)")

    ax.set_xlabel("Metro area side length (km)")
    ax.set_ylabel("Minimum trusted relays for full coverage")
    ax.set_title(f"Trusted-node cost of CV vs metro area   "
                 f"(N={args.n} users, {args.runs} layouts/area)",
                 fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11, loc="upper left")
    ax.set_ylim(bottom=-0.5)

    fig.tight_layout()
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or f"network_mink_N{args.n}_runs{args.runs}.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()
