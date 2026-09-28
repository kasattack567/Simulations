"""
topo_cost.py — Karavias-style cost comparison of DV vs CV on the six real
optical backbones, run in BOTH clock regimes and plotted side by side.

The clock assumption dominates the cost conclusion, so it is treated as the
primary axis of the study rather than a footnote. Every figure has two panels:

    left   DEPLOYED clocks  — DV 1 GHz, CV 100 MHz (commercial hardware today)
    right  MATCHED clocks   — both at the DV clock (isolates protocol physics
                              from the hardware clock)

Read the matched panel with care: at matched clocks a common bit/s target maps
to the SAME bits/channel-use threshold for both protocols, near where the two
rate curves cross, so the comparison is evaluated on top of the crossover and
is correspondingly sensitive (this is the degeneracy net_common warns about).
The deployed panel is the one that answers a procurement question; the matched
panel is a sensitivity check.

FIGURES (default output: ./figures_cost/)
  cost_fd_sites.png            RESULT. Fractional difference
                               FD = (C_DV - C_CV)/C_CV vs the DV/CV system cost
                               ratio R, at site cost s > 0, one curve per
                               topology, zero crossings marked.
                               Kept (rather than the s = 0 version) because at
                               s = 0 the identity FD = R/R* - 1 holds exactly,
                               so that figure carries nothing beyond the six
                               break-even numbers already in Table 2. At s > 0
                               the curves depart from it by up to ~0.24, and
                               the curvature is physical: DV re-optimises
                               toward fewer, longer hops as it gets pricier.
  cost_breakeven_vs_site.png   RESULT / headline. Break-even R* vs site cost s.
                               Log y so both regimes stay legible on one scale.
  cost_fom_distance_rate.png   METHODOLOGY. The distance-rate product l*r(l),
                               which sets the high-target asymptote and
                               explains the mechanism behind the break-even.

TABLES
  table1_provisioning.csv      What the model built: edges, fibre, machines and
                               trusted nodes per topology per protocol, in both
                               regimes, at baseline c_min.
  table2_breakeven.csv         What the model concluded: break-even R* vs
                               service target, both regimes.
  tables_cost.md               Both of the above as markdown, ready to paste.

  (--all-figures additionally emits the s = 0 FD panel and the R* vs c_min
   panel; both are redundant with Table 2 and are off by default.)

USAGE
  python topo_cost.py
  python topo_cost.py --cmin 1e8 --site 10
  python topo_cost.py --cv homodyne --all-figures
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cost_common as cc
from cost_common import (CostModel, build_tables, provision_network,
                         REACHABLE_BPS, TOPOLOGY_FILES, PROTO_LABEL)
from topo_loader import edge_lengths

plt.rcParams.update({
    "font.family": "serif", "font.size": 11, "axes.titlesize": 12,
    "axes.labelsize": 11, "legend.fontsize": 8.5, "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5, "axes.linewidth": 0.9, "lines.linewidth": 1.8,
    "figure.dpi": 120, "savefig.dpi": 200, "savefig.bbox": "tight",
})
DV_C, CV_C = "#1f4e9c", "#c0392b"
TOPO_COLORS = ["#1f4e9c", "#c0392b", "#2e8b57", "#8c564b", "#7d2d8c", "#e08214"]


def _fmt_hz(f):
    return f"{f/1e9:g} GHz" if f >= 1e9 else f"{f/1e6:g} MHz"


def _fmt_bps(b):
    return f"{b/1e6:g} Mbit/s" if b >= 1e6 else f"{b/1e3:g} kbit/s"


def _save(fig, outdir, fname):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, fname)
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved: {path}")


def _caption(fig, text):
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=8,
             style="italic", color="0.35", wrap=True)


# ============================================================
# A regime bundles one clock pair with its rate tables and cost model.
# ============================================================
class Regime:
    def __init__(self, label, dv_hz, cv_hz, cv_detection, dv_table=None):
        self.label = label
        self.dv_hz, self.cv_hz = float(dv_hz), float(cv_hz)
        self.tables, self.cv_key = build_tables(
            cv_detection, verbose=True, dv_hz=dv_hz, cv_hz=cv_hz,
            dv_table=dv_table)
        self.model = CostModel(self.tables, self.cv_key)

    @property
    def panel_title(self):
        if self.dv_hz == self.cv_hz:
            return f"Matched clocks — both {_fmt_hz(self.dv_hz)}"
        return (f"Deployed clocks — DV {_fmt_hz(self.dv_hz)}, "
                f"CV {_fmt_hz(self.cv_hz)}")

    def span(self, cmin):
        return (self.tables["dv"].l_at_rate(cmin),
                self.tables[self.cv_key].l_at_rate(cmin))


def _panels(figsize):
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    return fig, axes


# ============================================================
# FIGURE 1 (result): FD vs R, paired panels
# ============================================================
def fig_fd(regimes, cmin, s, rmin, rmax, outdir, band=None):
    R = np.linspace(rmin, rmax, 240)
    fig, axes = _panels((12.6, 5.4))
    stars = {}
    lo = hi = 0.0
    for ax, rg in zip(axes, regimes):
        st = {}
        for name, color in zip(TOPOLOGY_FILES, TOPO_COLORS):
            fd = np.array([rg.model.fd(name, cmin, r, s) for r in R])
            rstar = rg.model.break_even(name, cmin, s)
            st[name] = rstar
            lo, hi = min(lo, fd.min()), max(hi, fd.max())
            ax.plot(R, fd, color=color, label=f"{name}  ($R^*$={rstar:.2f})")
            if rmin <= rstar <= rmax:
                ax.plot([rstar], [0.0], "o", color=color, ms=5, zorder=5)
        if band:
            ax.axvspan(band[0], band[1], color="#7cc47f", alpha=0.15, zorder=0)
        ax.axhline(0.0, color="0.3", lw=1.0)
        ax.set_xlim(rmin, rmax)
        ax.set_xlabel("$R = c_{sys}^{DV} / c_{sys}^{CV}$")
        ax.set_title(rg.panel_title)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", framealpha=0.92, edgecolor="0.7",
                  ncol=1, fontsize=8)
        stars[rg.label] = st
    for ax in axes:
        ax.set_ylim(lo * 1.05, hi * 1.05)
        ax.text(0.985, 0.965, "CV cheaper", transform=ax.transAxes,
                ha="right", va="top", fontsize=9, color="0.35")
        ax.text(0.985, 0.03, "DV cheaper", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=9, color="0.35")
    axes[0].set_ylabel("Fractional cost difference  "
                       "$(C_{DV}-C_{CV})\\,/\\,C_{CV}$")
    site_txt = ("hardware only ($c_{site}=0$)" if s == 0 else
                f"site cost $c_{{site}}={s:g}$ CV-system units")
    fig.suptitle(f"Network cost comparison at {_fmt_bps(cmin)} per link — "
                 f"{site_txt}", fontsize=13, y=0.99)
    tail = (" At $c_{site}=0$ the identity FD = R/R*-1 holds exactly, so this "
            "panel adds nothing to the break-even table." if s == 0 else
            " Curvature is DV re-optimising toward fewer, longer hops as it "
            "gets pricier; panels share a y-axis.")
    _caption(fig,
             "Each edge carries c_min: n hops of L/n km, ceil(c_min / r(L/n)) parallel systems per hop, n-1 trusted-node "
             "sites; n minimises n*m*c_sys + (n-1)*c_site per protocol. Zero crossings R* are the break-even cost "
             f"ratios: for R > R* the CV network is the cheaper build. Hops must clear the {REACHABLE_BPS:.0f} bit/s "
             "validity floor." + tail)
    _save(fig, outdir, f"cost_fd_{'hw' if s == 0 else 'sites'}.png")
    return stars


# ============================================================
# FIGURE 2 (headline): break-even R* vs site cost, paired panels
# ============================================================
def fig_breakeven_vs_site(regimes, cmin, s_grid, outdir):
    rows = []
    fig, axes = _panels((12.6, 5.2))
    for ax, rg in zip(axes, regimes):
        for name, color in zip(TOPOLOGY_FILES, TOPO_COLORS):
            rs = [rg.model.break_even(name, cmin, s) for s in s_grid]
            ax.plot(s_grid, rs, "-o", ms=4, color=color, label=name)
            for s, r in zip(s_grid, rs):
                rows.append(dict(regime=rg.label, topology=name,
                                 site_cost_cv_units=s, break_even_R=r))
        ax.axhline(1.0, color="0.45", ls="--", lw=1.0)
        ax.set_yscale("log")
        ax.set_ylim(0.55, 13)
        ax.set_yticks([0.6, 1, 2, 3, 5, 8, 12])
        ax.get_yaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, p: f"{v:g}"))
        ax.set_xlabel("Trusted-node site cost  $c_{site}$  (CV-system units)")
        ax.set_title(rg.panel_title)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(ncol=2, framealpha=0.92, edgecolor="0.7", fontsize=8)
    axes[0].set_ylabel("Break-even system-cost ratio  $R^*$")
    fig.suptitle("How expensive may a DV system be before the CV network wins? "
                 f" ({_fmt_bps(cmin)} per link)", fontsize=13, y=0.99)
    _caption(fig,
             "For R > R* the CV network is cheaper, so a larger R* is a harder case for CV. Dashed line marks parity "
             "(R* = 1); below it CV wins even at equal system prices. Site costs weigh against CV because it builds "
             "coverage from more, shorter hops and carries ~3x the sites. Log y-axis, shared across panels, so the "
             "gap between regimes and the shape within each stay readable. Both protocols re-optimise their relay "
             "spacing at every point.")
    _save(fig, outdir, "cost_breakeven_vs_site.png")
    return pd.DataFrame(rows)


# ============================================================
# FIGURE 3 (methodology): distance-rate product, paired panels
# ============================================================
def fig_fom(regimes, cmin, outdir):
    fig, axes = _panels((12.6, 5.0))
    asyms = {}
    for ax, rg in zip(axes, regimes):
        marks = {}
        for proto, color in (("dv", DV_C), (rg.cv_key, CV_C)):
            L, prod, lstar, pstar = rg.tables[proto].distance_rate_product()
            ax.plot(L, prod, color=color, label=PROTO_LABEL[proto])
            ax.plot([lstar], [pstar], "o", color=color, ms=5)
            marks[proto] = (lstar, pstar)
        a = marks["dv"][1] / marks[rg.cv_key][1]
        asyms[rg.label] = (a, marks)
        ax.set_xlabel("Hop length  $\\ell$  (km)")
        ax.set_yscale("log")
        ax.set_ylim(1e4, 1e10)
        ax.set_title(f"{rg.panel_title}\nasymptotic $R^* \\to$ {a:.1f}",
                     fontsize=11)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="lower left", fontsize=8.5)
    axes[0].set_ylabel("$\\ell \\cdot r(\\ell)$   (bit$\\cdot$km / s per system)")
    fig.suptitle("Hardware figure of merit per system", fontsize=13, y=0.99)
    _caption(fig,
             "In the continuum limit (targets large against per-system rates, no site costs) an edge of length L needs "
             "~ L*c_min / max[l*r(l)] systems, so the break-even ratio tends to max_DV / max_CV. Matching the clocks "
             "lifts the CV curve by the clock ratio and inverts the ordering, which is why the two regimes reach "
             "opposite asymptotes. At moderate targets provisioning sits in the one-system-per-hop regime and the "
             "break-even follows the ratio of the longest hops sustaining c_min instead.")
    _save(fig, outdir, "cost_fom_distance_rate.png")
    return asyms


# ============================================================
# TABLE 1 — provisioning, both regimes
# ============================================================
def table1_provisioning(regimes, cmin, outdir):
    rows = []
    for rg in regimes:
        for name in TOPOLOGY_FILES:
            topo = rg.model.topos[name]
            fkm = float(np.sum(edge_lengths(topo, "real")))
            rec = dict(regime=rg.label, topology=name,
                       edges=len(edge_lengths(topo, "real")),
                       fibre_km=round(fkm))
            for proto, pre in (("dv", "DV"), (rg.cv_key, "CV")):
                pn = provision_network(topo, rg.tables[proto], cmin, 0.0)
                rec[f"{pre}_machines"] = pn["systems"]
                rec[f"{pre}_trusted_nodes"] = pn["relays"]
                rec[f"{pre}_mean_hop_km"] = round(
                    float(np.mean([e["hop_km"] for e in pn["per_edge"]])), 1)
            rec["machines_CV_over_DV"] = round(
                rec["CV_machines"] / rec["DV_machines"], 2)
            rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, "table1_provisioning.csv"), index=False)
    return df


# ============================================================
# TABLE 2 — break-even vs service target, both regimes
# ============================================================
def table2_breakeven(regimes, cmin_grid, outdir, s=0.0):
    rows = []
    for rg in regimes:
        for c in cmin_grid:
            rec = dict(regime=rg.label, cmin_bps=c, site_cost=s)
            for name in TOPOLOGY_FILES:
                rec[name] = round(rg.model.break_even(name, c, s), 2)
            vals = [rec[n] for n in TOPOLOGY_FILES]
            rec["range"] = f"{min(vals):.2f}–{max(vals):.2f}"
            rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, "table2_breakeven.csv"), index=False)
    return df


def _to_md(df):
    """Markdown table without the optional `tabulate` dependency."""
    cols = [str(c) for c in df.columns]
    body = [[("" if pd.isna(v) else str(v)) for v in row]
            for row in df.itertuples(index=False)]
    num = [all(_isnum(r[i]) for r in body) if body else False
           for i in range(len(cols))]
    w = [max(len(cols[i]), *(len(r[i]) for r in body)) if body else len(cols[i])
         for i in range(len(cols))]
    def line(cells):
        return "| " + " | ".join(
            c.rjust(w[i]) if num[i] else c.ljust(w[i])
            for i, c in enumerate(cells)) + " |"
    sep = "|" + "|".join(("-" * (w[i] + 1)) + (":" if num[i] else " ")
                         for i in range(len(cols))) + "|"
    return "\n".join([line(cols), sep] + [line(r) for r in body])


def _isnum(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def write_markdown(df1, df2, cmin, s, outdir, regimes):
    p = os.path.join(outdir, "tables_cost.md")
    with open(p, "w") as f:
        f.write(f"# Cost model tables (c_min = {_fmt_bps(cmin)} per link)\n\n")
        f.write("## Table 1 — provisioning at hardware-optimal spacing "
                "(c_site = 0)\n\n")
        f.write(_to_md(df1))
        f.write("\n\n## Table 2 — break-even system-cost ratio R* vs service "
                f"target (c_site = {s:g})\n\n")
        f.write("For R > R*, the CV network is the cheaper build.\n\n")
        f.write(_to_md(df2))
        f.write("\n\nSpans at each target (longest hop sustaining c_min):\n\n")
        for rg in regimes:
            d, c = rg.span(cmin)
            f.write(f"- {rg.panel_title}: DV {d:.1f} km, CV {c:.1f} km "
                    f"(ratio {d/c:.2f})\n")
    print(f"Saved: {p}")


# ============================================================
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cmin", type=float, default=10e6)
    p.add_argument("--site", type=float, default=5.0)
    p.add_argument("--site-grid", type=float, nargs="+",
                   default=[0, 1, 2, 5, 10, 20])
    p.add_argument("--cmin-grid", type=float, nargs="+",
                   default=[1e6, 3e6, 1e7, 3e7, 1e8])
    p.add_argument("--cv", choices=["heterodyne", "homodyne"],
                   default="heterodyne")
    p.add_argument("--dv-clock", type=float, default=1e9)
    p.add_argument("--cv-deployed-clock", type=float, default=100e6)
    p.add_argument("--rmin", type=float, default=0.1)
    p.add_argument("--rmax", type=float, default=10.0)
    p.add_argument("--band", type=float, nargs=2, default=None,
                   metavar=("LO", "HI"))
    p.add_argument("--all-figures", action="store_true",
                   help="also emit the two figures that are redundant with "
                        "Table 2 (s=0 FD panel, R* vs c_min panel)")
    p.add_argument("--outdir", type=str, default="figures_cost")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("DV vs CV network cost model — both clock regimes")
    print(f"  c_min {_fmt_bps(args.cmin)}   CV {args.cv}   "
          f"floor {REACHABLE_BPS:.0f} bit/s")

    dep = Regime("deployed", args.dv_clock, args.cv_deployed_clock, args.cv)
    mat = Regime("matched", args.dv_clock, args.dv_clock, args.cv,
                 dv_table=dep.tables["dv"])          # DV clock is unchanged
    regimes = [dep, mat]
    for rg in regimes:
        d, c = rg.span(args.cmin)
        print(f"  {rg.label:>8}: longest hop at c_min — DV {d:.1f} km, "
              f"CV {c:.1f} km (ratio {d/c:.2f})")

    df1 = table1_provisioning(regimes, args.cmin, args.outdir)
    df2 = table2_breakeven(regimes, args.cmin_grid, args.outdir)
    print("\nTABLE 1 — provisioning")
    print(df1.to_string(index=False))
    print("\nTABLE 2 — break-even R* vs target")
    print(df2.to_string(index=False))

    asyms = fig_fom(regimes, args.cmin, args.outdir)
    fig_fd(regimes, args.cmin, args.site, args.rmin, args.rmax, args.outdir,
           band=args.band)
    df_site = fig_breakeven_vs_site(regimes, args.cmin, args.site_grid,
                                    args.outdir)
    df_site.to_csv(os.path.join(args.outdir, "cost_breakeven_vs_site.csv"),
                   index=False)
    if args.all_figures:
        fig_fd(regimes, args.cmin, 0.0, args.rmin, args.rmax, args.outdir,
               band=args.band)
    write_markdown(df1, df2, args.cmin, 0.0, args.outdir, regimes)

    print("\nBreak-even R* at s=%g:" % args.site)
    for rg in regimes:
        vals = [rg.model.break_even(n, args.cmin, args.site)
                for n in TOPOLOGY_FILES]
        print(f"  {rg.label:>8}: {min(vals):.2f}–{max(vals):.2f}   "
              f"(asymptote {asyms[rg.label][0]:.1f})")
    print("\nReading: for R above R*, the CV network is the cheaper build.")


if __name__ == "__main__":
    main()