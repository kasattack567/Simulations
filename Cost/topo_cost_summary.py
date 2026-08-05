"""
topo_cost_summary.py — the cost model as one figure.

THE PROBLEM THIS SOLVES. The model has two free assumptions that nobody can
pin down from the literature: the service target c_min (how much key each link
must carry) and the trusted-node site cost c_site. Plotting a curve against one
of them means fixing the other at an invented value, which is exactly the kind
of unsupported constant a reviewer will challenge. So neither is fixed here.

The figure maps the break-even system-cost ratio R* over the whole (c_min,
c_site) plane, one panel per clock regime. The reader supplies their own
assumptions and reads the answer off the map. Contour lines are labelled with
R* directly.

Reading: R* is the break-even system-cost ratio. For R > R*, the CV network is
the cheaper build. The R* = 1 contour, where present, is parity — CV wins there
even at equal system prices.

Topology is collapsed to the median over the six backbone networks; it shifts
R* by only 1-33% (typically ~10%), against a factor of ~5 for the clock regime.
Run with --spread to also print the per-topology spread.
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cost_common import CostModel, build_tables, TOPOLOGY_FILES

plt.rcParams.update({
    "font.family": "serif", "font.size": 11, "axes.titlesize": 11.5,
    "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "axes.linewidth": 0.9, "figure.dpi": 120, "savefig.dpi": 200,
    "savefig.bbox": "tight",
})

LEVELS = [1, 1.5, 2, 3, 4, 6, 8, 10, 14, 20]


def grid(model, c_vals, s_vals):
    Z = np.empty((len(s_vals), len(c_vals)))
    for i, s in enumerate(s_vals):
        for j, c in enumerate(c_vals):
            Z[i, j] = np.median([model.break_even(n, c, s)
                                 for n in TOPOLOGY_FILES])
    return Z


def report_spread(regimes, c_vals, s_vals):
    """Print the per-topology spread of R* over the six networks.

    Spread is reported as the half-range about the median, in percent:
    100 * (max - min) / (2 * median). This is what backs the claim that the
    six topologies agree to within about +/-10% except at the extreme corner.
    """
    cs_show = [1e6, 1e7, 1e8]
    ss_show = [0.0, 5.0, 20.0]
    print("\n" + "=" * 78)
    print("PER-TOPOLOGY SPREAD of R*  (min / median / max over the six networks)")
    print("=" * 78)
    for _title, model, key in regimes:
        print(f"\n{key} clocks")
        print(f"{'c_min':>7} {'c_site':>7} | {'min':>6} {'med':>6} {'max':>6}"
              f" {'spread':>8} | per-topology")
        print("-" * 82)
        for c in cs_show:
            for s in ss_show:
                vals = {n: model.break_even(n, c, s) for n in TOPOLOGY_FILES}
                arr = np.array(list(vals.values()))
                lo, med, hi = arr.min(), np.median(arr), arr.max()
                spread_pct = 100.0 * (hi - lo) / (2.0 * med) if med else 0.0
                detail = "  ".join(
                    f"{n.split('_')[-1][:4]}:{v:.1f}" for n, v in vals.items())
                print(f"{c/1e6:>6.0f}M {s:>7.1f} | "
                      f"{lo:>6.2f} {med:>6.2f} {hi:>6.2f} {spread_pct:>6.1f}% "
                      f"| {detail}")

    # worst-case spread over the full grid, with where it occurs
    worst, where = 0.0, None
    for _title, model, key in regimes:
        for c in c_vals:
            for s in s_vals:
                arr = np.array([model.break_even(n, c, s)
                                for n in TOPOLOGY_FILES])
                med = np.median(arr)
                if med:
                    sp = 100.0 * (arr.max() - arr.min()) / (2.0 * med)
                    if sp > worst:
                        worst, where = sp, (key, c, s)
    k, c, s = where
    print(f"\nWorst-case spread over the full grid: {worst:.1f}% "
          f"(half-range / median) at {k} clocks, "
          f"c_min={c/1e6:.0f} Mbps, c_site={s:.1f}.")
    print("Elsewhere the six topologies agree far more tightly, so the median "
          "is a faithful summary.")


def plot_slices(regimes, outdir):
    """Six per-topology R* curves at two site-cost cuts, matched clocks.

    Supporting figure for the median heatmap: shows the six topologies track
    closely (so the median is faithful), and where the spread opens up.
    """
    model = dict((k, m) for _t, m, k in regimes)["matched"]
    c_vals = np.logspace(6, 8, 40)
    cuts = [(0.0, r"(a) No site cost ($c_{site}=0$)"),
            (20.0, r"(b) High site cost ($c_{site}=20$)")]
    colors = plt.cm.tab10(np.linspace(0, 1, len(TOPOLOGY_FILES)))

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), sharey=True)
    for ax, (s, title) in zip(axes, cuts):
        for (name, _), col in zip(TOPOLOGY_FILES.items(), colors):
            R = [model.break_even(name, c, s) for c in c_vals]
            ax.plot(c_vals, R, color=col, lw=1.6, label=name)
        ax.axhline(1.0, color="0.4", ls="--", lw=0.9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Service target per link, $c_{min}$ (bit/s)")
        ax.set_title(title)
    axes[0].set_ylabel("Break-even cost ratio, $R^{*}$")
    axes[0].legend(fontsize=8.5, ncol=2, frameon=False, loc="upper right")
    fig.suptitle("Per-topology break-even ratio (matched clocks): "
                 "the six networks track closely", fontsize=12.5, y=1.02)
    fig.text(0.5, -0.04,
             "Each line is one of the six backbone topologies. At $c_{site}=0$ the curves are nearly coincident, so "
             "the median heatmap is a faithful summary; the spread widens only at high site cost, driven by SAGO, the "
             "one topology whose DV network needs no relays.",
             ha="center", va="top", fontsize=8.5, style="italic",
             color="0.35", wrap=True)
    path = os.path.join(outdir, "cost_slices.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved: {path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nc", type=int, default=13, help="c_min grid points")
    p.add_argument("--ns", type=int, default=21, help="c_site grid points")
    p.add_argument("--smax", type=float, default=20.0)
    p.add_argument("--outdir", type=str, default="figures_cost")
    p.add_argument("--spread", action="store_true",
                   help="also print the per-topology break-even spread")
    p.add_argument("--slices", action="store_true",
                   help="also save the per-topology slice figure (cost_slices.png)")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    dep_t, cvk = build_tables("heterodyne", dv_hz=1e9, cv_hz=100e6)
    mat_t, _ = build_tables("heterodyne", dv_hz=1e9, cv_hz=1e9,
                            dv_table=dep_t["dv"])
    regimes = [("(a) Deployed clocks: DV 1 GHz, CV 100 MHz",
                CostModel(dep_t, cvk), "deployed"),
               ("(b) Matched clocks: both 1 GHz",
                CostModel(mat_t, cvk), "matched")]

    c_vals = np.logspace(4, 8, args.nc)
    s_vals = np.linspace(0, args.smax, args.ns)
    C, S = np.meshgrid(c_vals, s_vals)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), sharey=True)
    rows = []
    Zs = [grid(m, c_vals, s_vals) for _, m, _ in regimes]
    vmin = min(Z.min() for Z in Zs)
    vmax = max(Z.max() for Z in Zs)
    norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

    for ax, (title, model, key), Z in zip(axes, regimes, Zs):
        im = ax.pcolormesh(C, S, Z, norm=norm, cmap="RdYlBu_r",
                           shading="gouraud")
        cs = ax.contour(C, S, Z, levels=LEVELS, colors="0.15",
                        linewidths=0.7)
        ax.clabel(cs, inline=True, fontsize=9, fmt=lambda v: f"{v:g}")
        ax.set_xscale("log")
        ax.set_xlabel("Service target per link, $c_{min}$ (bit/s)")
        ax.set_title(title)
        for i, s in enumerate(s_vals):
            for j, c in enumerate(c_vals):
                rows.append(dict(regime=key, cmin_bps=c, site_cost=s,
                                 break_even_R=round(Z[i, j], 3)))
    axes[0].set_ylabel("Trusted-node site cost, $c_{site}$\n(CV-system units)")

    cb = fig.colorbar(im, ax=axes, fraction=0.035, pad=0.02,
                      ticks=[1, 2, 3, 5, 8, 12, 20])
    cb.ax.set_yticklabels([f"{v:g}" for v in [1, 2, 3, 5, 8, 12, 20]])
    cb.set_label("Break-even cost ratio, $R^{*}$")

    fig.suptitle("Break-even system-cost ratio for CV relative to DV QKD "
                 "networks", fontsize=12.5, y=1.0)
    fig.text(0.5, -0.06,
             "Contours give the factor by which one DV system must exceed one CV system in cost before the CV network "
             "is the cheaper build. Neither free assumption is fixed: the service target varies along the horizontal "
             "axis and the trusted-node site cost along the vertical, so a reader may read the result off the map "
             "under their own assumptions. Median over the six backbone topologies.",
             ha="center", va="top", fontsize=8.5, style="italic", color="0.35",
             wrap=True)
    path = os.path.join(args.outdir, "cost_summary.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved: {path}")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.outdir, "table_summary.csv"), index=False)

    # readable extract at round values
    print("\nBreak-even R* (median over the six topologies)\n")
    cs_show = [1e6, 1e7, 1e8]
    ss_show = [0, 5, 10, 20]
    for _, _, key in regimes:
        print(f"{key} clocks")
        print(f"{'c_site':>8} | " + " | ".join(f"{c/1e6:>6.0f} Mbps"
                                               for c in cs_show))
        print("-" * 42)
        for s in ss_show:
            si = int(np.argmin(np.abs(s_vals - s)))
            cells = []
            for c in cs_show:
                ci = int(np.argmin(np.abs(c_vals - c)))
                z = df[(df.regime == key)].pivot(index="site_cost",
                                                 columns="cmin_bps",
                                                 values="break_even_R")
                cells.append(f"{z.iloc[si, ci]:>11.2f}")
            print(f"{s_vals[si]:>8.0f} | " + " | ".join(cells))
        print()
    print("For R above the quoted ratio, the CV network is the cheaper build.")

    if args.spread:
        report_spread(regimes, c_vals, s_vals)

    if args.slices:
        plot_slices(regimes, args.outdir)


if __name__ == "__main__":
    main()