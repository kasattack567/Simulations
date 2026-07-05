"""
Network sweep — how DV and CV perform as the metro area grows.

For each area size (10..80 km, 10 km steps), generate SWEEP_RUNS random user
layouts (fixed N users), compute each layout's network metrics, and average.
Because the user PLACEMENT is random (the engine itself is deterministic), there
is genuine run-to-run spread here — shown as a +/- 1 std band. This is unlike the
point-to-point / sensitivity plots, which have no randomness and no error bars.

Two panels:
  (left)  coverage  = fraction of user pairs that can share a key
  (right) total key rate = sum of direct-link rates over ALL pairs
                            (failed links contribute zero — no survivorship bias)

Curves: DV, CV-heterodyne, CV-homodyne. Direct links only (no relays).

Usage:
    python network_sweep.py
    python network_sweep.py --n 30 --runs 20
    python network_sweep.py --areas 10 20 30 40 50 60 70 80
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from net_common import place_users, all_pairs, pair_distance_km, link_rate

PROTOCOLS = [("dv", "DV — decoy BB84", "#1f4e9c", "-o"),
             ("cv_het", "CV — heterodyne", "#c0392b", "-s"),
             ("cv_hom", "CV — homodyne", "#7d2d8c", "--D")]


def network_metrics(users, protocol):
    """Coverage and total key rate (over all pairs) for one layout."""
    n = len(users)
    n_pairs = n * (n - 1) // 2
    total, connected = 0.0, 0
    for i, j in all_pairs(n):
        r = link_rate(pair_distance_km(users, i, j), protocol)
        if r > 1e-9:
            total += r
            connected += 1
    return connected / n_pairs, total   # coverage, total rate


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=20, help="Number of users (fixed)")
    p.add_argument("--runs", type=int, default=10, help="Random layouts per area")
    p.add_argument("--areas", type=float, nargs="+",
                   default=[10, 20, 30, 40, 50, 60, 70, 80],
                   help="Metro area sizes (km) to sweep")
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    areas = np.array(args.areas, dtype=float)
    # results[proto] = dict(cov_mean, cov_std, tot_mean, tot_std) arrays over areas
    results = {proto: {k: np.zeros(len(areas)) for k in
                       ("cov_mean", "cov_std", "tot_mean", "tot_std")}
               for proto, *_ in PROTOCOLS}

    print(f"N={args.n} users, {args.runs} random layouts per area\n")
    for a_idx, area in enumerate(areas):
        # one set of seeds per area, shared across protocols (same layouts)
        seeds = range(1000 * a_idx, 1000 * a_idx + args.runs)
        layouts = [place_users(args.n, area_km=area, seed=s) for s in seeds]
        line = f"area {area:4.0f} km:"
        for proto, label, *_ in PROTOCOLS:
            covs, tots = [], []
            for users in layouts:
                c, t = network_metrics(users, proto)
                covs.append(c); tots.append(t)
            results[proto]["cov_mean"][a_idx] = np.mean(covs)
            results[proto]["cov_std"][a_idx]  = np.std(covs)
            results[proto]["tot_mean"][a_idx] = np.mean(tots)
            results[proto]["tot_std"][a_idx]  = np.std(tots)
            line += f"  {proto} cov={np.mean(covs)*100:3.0f}%"
        print(line)

    # ============================================================
    # PLOT — two panels
    # ============================================================
    plt.rcParams.update({"font.family": "serif", "font.size": 12})
    fig, (axC, axR) = plt.subplots(1, 2, figsize=(15, 6))

    for proto, label, color, style in PROTOCOLS:
        r = results[proto]
        ls = style[:-1] if style[-1] in "osD^" else style
        marker = style[-1] if style[-1] in "osD^" else None
        # coverage panel
        axC.plot(areas, r["cov_mean"] * 100, style, color=color, label=label,
                 markersize=5)
        axC.fill_between(areas, (r["cov_mean"] - r["cov_std"]) * 100,
                         (r["cov_mean"] + r["cov_std"]) * 100,
                         color=color, alpha=0.15)
        # total-rate panel
        axR.plot(areas, r["tot_mean"], style, color=color, label=label,
                 markersize=5)
        axR.fill_between(areas, np.maximum(r["tot_mean"] - r["tot_std"], 1e-12),
                         r["tot_mean"] + r["tot_std"],
                         color=color, alpha=0.15)

    axC.set_xlabel("Metro area side length (km)")
    axC.set_ylabel("Coverage (% of user pairs with a key)")
    axC.set_title("Network coverage vs metro area")
    axC.set_ylim(0, 105)
    axC.grid(True, alpha=0.3)
    axC.legend(fontsize=10)

    axR.set_xlabel("Metro area side length (km)")
    axR.set_ylabel("Total network key rate (bits / channel use)")
    axR.set_title("Total key rate vs metro area")
    axR.set_yscale("log")
    axR.grid(True, alpha=0.3, which="both")
    axR.legend(fontsize=10)

    fig.suptitle(
        f"DV vs CV over random metro networks   "
        f"(N={args.n} users, {args.runs} layouts/area, direct links)",
        fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or f"network_sweep_N{args.n}_runs{args.runs}.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()
