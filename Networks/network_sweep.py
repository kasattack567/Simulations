"""
Network sweep — how DV and CV perform as the area grows.

For each area size, generate --runs random user layouts (fixed N users), compute
each layout's network metrics, and average. Because the user PLACEMENT is random
(the engine itself is deterministic), there is genuine run-to-run spread here —
shown as a +/- 1 std band. This is unlike the point-to-point / sensitivity plots,
which have no randomness and no error bars.

Two panels:
  (left)  coverage  = fraction of user pairs that can share a key
  (right) total key rate in BITS/S = sum of direct-link rates over ALL pairs
                            (failed links contribute zero — no survivorship bias)

Curves: DV, CV-heterodyne, CV-homodyne. Direct links only (no relays).

REACHABILITY FLOOR. A pair counts as connected only if its rate clears
net_common.REACHABLE_BPS. The default is 256 bit/s — one AES-256 key per second,
the smallest floor with a stated operational meaning (Nweke et al.,
arXiv:2306.15031, define the required rate as key length x refresh rate). The
floor is NOT neutral between the protocols: CV's rate-distance curve falls
vertically at its limit because eps_b/(eta*T) diverges, so CV coverage barely
moves with it, while DV's tail decays gently and DV coverage moves more. Raising
the floor is close to a DV-only handicap. The coverage panel labels the floor it
used, so a figure cannot be read without it. Sweep it per run:

    QKD_REACHABLE_BPS=20   python network_sweep.py   # Jouguet CV field test
    QKD_REACHABLE_BPS=256  python network_sweep.py   # default, 1 AES-256 key/s
    QKD_REACHABLE_BPS=1e4  python network_sweep.py   # carrier-grade

Usage:
    python network_sweep.py
    python network_sweep.py --n 30 --runs 20
    python network_sweep.py --areas 10 20 30 40 50 60 70 80

Layouts are seeded per area (1000*area_index + run), so the same --n and --runs
reproduce the same geometry exactly. A rerun after a parameter change therefore
isolates the parameter change, with the user placement held fixed.
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from net_common import (place_users, all_pairs, pair_distance_km, link_rate,
                        to_bps, is_reachable, REACHABLE_BPS, CLOCK_LABEL)

PROTOCOLS = [("dv", "DV — decoy BB84", "#1f4e9c", "-o"),
             ("cv_het", "CV — heterodyne", "#c0392b", "-s"),
             ("cv_hom", "CV — homodyne", "#7d2d8c", "--D")]

# OUTPUT UNITS: bits/s throughout. Clocks come from net_common (single source of
# truth) and are MATCHED at 1 GHz for DV and CV. The old --units flag is gone:
# with matched clocks the bits/channel-use view is this figure rescaled by a
# single common factor, so it carried no information the bits/s view lacks.
# Run with QKD_CV_CLOCK_HZ=100e6 for the deployed-clock version.


def _fmt_floor(bps):
    """Format the floor for a label without forcing it into kbit/s.

    The previous version divided by 1e3 unconditionally, which rendered the
    256 bit/s default as '0.256 kbit/s'.
    """
    if bps >= 1e6:
        return f"{bps/1e6:g} Mbit/s"
    if bps >= 1e3:
        return f"{bps/1e3:g} kbit/s"
    return f"{bps:g} bit/s"


def network_metrics(users, protocol):
    """Coverage and total key rate in BITS/S (over all pairs) for one layout."""
    n = len(users)
    n_pairs = n * (n - 1) // 2
    total, connected = 0.0, 0
    for i, j in all_pairs(n):
        r = link_rate(pair_distance_km(users, i, j), protocol)
        if is_reachable(r, protocol):
            total += r
            connected += 1
    return connected / n_pairs, to_bps(total, protocol)   # coverage, bits/s


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=20, help="Number of users (fixed)")
    p.add_argument("--runs", type=int, default=10, help="Random layouts per area")
    p.add_argument("--areas", type=float, nargs="+",
                   default=[10, 50, 100, 150, 200, 250],
                   help="Area sizes (km) to sweep")
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    areas = np.array(args.areas, dtype=float)
    # results[proto] = dict(cov_mean, cov_std, tot_mean, tot_std) arrays over areas
    results = {proto: {k: np.zeros(len(areas)) for k in
                       ("cov_mean", "cov_std", "tot_mean", "tot_std")}
               for proto, *_ in PROTOCOLS}

    print(f"N={args.n} users, {args.runs} random layouts per area")
    print(f"reachability floor: {_fmt_floor(REACHABLE_BPS)} "
          f"({REACHABLE_BPS/256:.2f} AES-256 keys/s)   clocks: {CLOCK_LABEL}\n")
    for a_idx, area in enumerate(areas):
        # one set of seeds per area, shared across protocols (same layouts)
        seeds = range(1000 * a_idx, 1000 * a_idx + args.runs)
        layouts = [place_users(args.n, area_km=area, seed=s) for s in seeds]
        line = f"area {area:4.0f} km:"
        for proto, label, *_ in PROTOCOLS:
            covs, tots = [], []
            for users in layouts:
                c, t = network_metrics(users, proto)
                covs.append(c)
                tots.append(t)
            results[proto]["cov_mean"][a_idx] = np.mean(covs)
            results[proto]["cov_std"][a_idx] = np.std(covs)
            results[proto]["tot_mean"][a_idx] = np.mean(tots)
            results[proto]["tot_std"][a_idx] = np.std(tots)
            line += f"  {proto} cov={np.mean(covs)*100:3.0f}%"
        print(line)

    # ============================================================
    # PLOT — two panels
    # ============================================================
    plt.rcParams.update({"font.family": "serif", "font.size": 12})
    fig, (axC, axR) = plt.subplots(1, 2, figsize=(15, 6))

    for proto, label, color, style in PROTOCOLS:
        r = results[proto]
        # coverage panel
        axC.plot(areas, r["cov_mean"] * 100, style, color=color, label=label,
                 markersize=5)
        axC.fill_between(areas, (r["cov_mean"] - r["cov_std"]) * 100,
                         (r["cov_mean"] + r["cov_std"]) * 100,
                         color=color, alpha=0.15)
        # total-rate panel — already in bits/s (converted in network_metrics)
        axR.plot(areas, r["tot_mean"], style, color=color, label=label,
                 markersize=5)
        axR.fill_between(areas,
                         np.maximum(r["tot_mean"] - r["tot_std"], 1e-12),
                         r["tot_mean"] + r["tot_std"],
                         color=color, alpha=0.15)

    axC.set_xlabel("Area side length (km)")
    axC.set_ylabel(f"Coverage (% of pairs at $\\geq$ {_fmt_floor(REACHABLE_BPS)})")
    axC.set_title("Network coverage vs area")
    axC.set_ylim(0, 105)
    axC.grid(True, alpha=0.3)
    axC.legend(fontsize=10)

    axR.set_xlabel("Area side length (km)")
    axR.set_ylabel("Total network key rate (bits / s)")
    axR.set_title("Total key rate vs area")
    axR.set_yscale("log")
    axR.grid(True, alpha=0.3, which="both")
    # Limits from the plotted data. The per-pair floor is not drawn here: this
    # panel sums over pairs, so the floor is not a threshold for these curves,
    # and including it stretched the log axis over several empty decades.
    lows = [r["tot_mean"][r["tot_mean"] > 0].min() for r in results.values()
            if np.any(r["tot_mean"] > 0)]
    highs = [(r["tot_mean"] + r["tot_std"]).max() for r in results.values()]
    if lows:
        axR.set_ylim(min(lows) / 3, max(highs) * 3)
    axR.legend(fontsize=10)

    fig.suptitle(
        f"DV vs CV over random networks   "
        f"(N={args.n} users, {args.runs} layouts/area, direct links, "
        f"{_fmt_floor(REACHABLE_BPS)} floor)",
        fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    # Floor goes in the filename: a 256 bit/s run and a 1 kbit/s run are different
    # figures and must not silently overwrite each other.
    default = (f"network_sweep_N{args.n}_runs{args.runs}"
               f"_floor{REACHABLE_BPS:.0f}bps_bps.png")
    path = os.path.join(outdir, args.save or default)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()