"""
Network sweep 2 — how DV and CV scale with the NUMBER of users, at fixed area.

Fixes the metro area at 40 km (the upper end of the realistic deployed-metro
band; Madrid links 1.9-33 km, Tokyo ~45 km) and sweeps the number of users from
small up to 40. Why 40: the largest deployed metro QKD network is Hefei with
40 user nodes (Chen et al., npj QI 7, 134, 2021); Jinan reached ~56. So 40 users
is a realistic ceiling for a metropolitan QKD network.

For each user count, USER_RUNS random layouts are generated and the network's
average key rate (over ALL pairs, failed links = 0) is averaged across them.
The +/- 1 std band reflects run-to-run variation in random placement (the
engine itself is deterministic).

y-axis: average key rate per user pair (bits / s). Rates come out of the
engines per channel use and are converted with net_common.to_bps at the
protocol's clock; DV and CV are both clocked at 1 GHz (matched).

Usage:
    python network_users.py
    python network_users.py --area 40 --runs 5 --umax 40
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from net_common import (place_users, all_pairs, pair_distance_km, link_rate,
                        to_bps, is_reachable)

PROTOCOLS = [("dv", "DV — decoy BB84", "#1f4e9c", "-o"),
             ("cv_het", "CV — heterodyne", "#c0392b", "-s"),
             ("cv_hom", "CV — homodyne", "#7d2d8c", "--D")]


def avg_rate_all_pairs(users, protocol):
    """Average key rate in BITS/S over ALL pairs (failed links count as zero)."""
    n = len(users)
    n_pairs = n * (n - 1) // 2
    if n_pairs == 0:
        return 0.0
    total = 0.0
    for i, j in all_pairs(n):
        r = link_rate(pair_distance_km(users, i, j), protocol)
        if is_reachable(r, protocol):
            total += r
    return to_bps(total / n_pairs, protocol)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--area", type=float, default=40.0,
                   help="Fixed metro area side (km); default 40 (deployed-metro ceiling)")
    p.add_argument("--runs", type=int, default=5, help="Random layouts per user count")
    p.add_argument("--umin", type=int, default=3, help="Min users")
    p.add_argument("--umax", type=int, default=40, help="Max users (Hefei = 40)")
    p.add_argument("--ustep", type=int, default=2, help="User-count step")
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    user_counts = list(range(args.umin, args.umax + 1, args.ustep))
    results = {proto: {"mean": np.zeros(len(user_counts)),
                       "std": np.zeros(len(user_counts))}
               for proto, *_ in PROTOCOLS}

    print(f"Area {args.area:.0f} km fixed, {args.runs} layouts per user count, "
          f"users {args.umin}..{args.umax}\n")
    for u_idx, n_users in enumerate(user_counts):
        seeds = range(5000 * u_idx, 5000 * u_idx + args.runs)
        layouts = [place_users(n_users, area_km=args.area, seed=s) for s in seeds]
        line = f"N={n_users:3d}:"
        for proto, label, *_ in PROTOCOLS:
            vals = [avg_rate_all_pairs(u, proto) for u in layouts]
            results[proto]["mean"][u_idx] = np.mean(vals)
            results[proto]["std"][u_idx] = np.std(vals)
            line += f"  {proto}={np.mean(vals):.3e}"
        print(line)

    # ============================================================
    # PLOT
    # ============================================================
    plt.rcParams.update({"font.family": "serif", "font.size": 12})
    fig, ax = plt.subplots(figsize=(10, 6.2))

    x = np.array(user_counts)
    for proto, label, color, style in PROTOCOLS:
        m = results[proto]["mean"]
        s = results[proto]["std"]
        ax.plot(x, m, style, color=color, label=label, markersize=5)
        ax.fill_between(x, np.maximum(m - s, 1e-12), m + s, color=color, alpha=0.15)

    ax.set_xlabel("Number of users in the network")
    ax.set_ylabel("Average key rate per pair  (bits / s)")
    ax.set_yscale("log")
    ax.set_title(f"DV vs CV — average key rate vs network size "
                 f"(fixed {args.area:.0f}x{args.area:.0f} km metro area)",
                 fontsize=13)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=11)

    # mark the Hefei deployed benchmark
    if args.umax >= 40 >= args.umin:
        ax.axvline(40, color="0.5", ls=":", lw=1.2)
        ax.text(40, ax.get_ylim()[1], " Hefei (40 users)", rotation=90,
                va="top", ha="right", color="0.4", fontsize=9)

    fig.tight_layout()
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or f"network_users_area{args.area:.0f}_runs{args.runs}.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()