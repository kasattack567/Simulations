"""
uncertainty_analysis.py — uncertainty quantification and global sensitivity
analysis for the DV/CV comparison.

WHAT THIS ANSWERS. The sensitivity sweeps (s1..s6) vary one parameter at a time
around a fixed baseline. That identifies which parameter matters but says
nothing about what the headline numbers are worth when every parameter is
uncertain at once, which is the situation an operator is actually in. This
module samples the full joint parameter space over the deployed bands and
propagates it through to the quantities the report reports:

    crossover distance          where the advantage passes from CV to DV
    DV usable reach             at the Eq.(RSKR) floor
    CV usable reach             at the same floor
    span ratio  DV/CV           the driver of the trusted-node budget
    relay counts per topology   the trusted-node budget itself
    relay ratio CV/DV           the headline of the work

and then asks three questions of the result: how wide are the distributions,
which inputs drive them, and why is the span ratio narrower than the crossover.

METHOD.
  Sampling      Latin hypercube over the deployed bands of Table 2, scrambled
                and seeded. LHS is used rather than plain Monte Carlo because
                it stratifies every marginal, so a given accuracy needs roughly
                an order of magnitude fewer engine evaluations -- and the DV
                engine costs ~64 ms per call, which is what bounds this study.
  Priors        Uniform over each band, except the dark-count rate, which is
                log-uniform because its band spans three decades and a uniform
                draw would place 90% of the mass in the top decade. Bands are
                the deployed ranges cited in Table 2. NOTE that alpha uses the
                deployed fibre band (0.18-0.21), NOT the wider hollow-core
                sweep axis, which is a what-if rather than a deployed range.
  Independence  Parameters are drawn independently. This is a real assumption:
                detector efficiency and electronic noise are plausibly
                correlated within a vendor's product line, and correlation
                would narrow the output distributions relative to what is
                reported here. The results are therefore conservative.
  Estimators    Percentiles with bootstrap confidence intervals, so the
                reported uncertainty carries its own uncertainty.
  Sensitivity   Standardised rank regression coefficients (SRRC) with the rank
                model's R^2 reported alongside. SRRC is used rather than raw
                correlation because it accounts for the other inputs, and rank
                rather than linear because the responses are monotone but not
                linear. Where R^2 is low the SRRC values should be distrusted
                and a variance-based method used instead.

RESUMABILITY. One sample costs ~1.9 s, almost all of it in the DV engine's
intensity optimisation. The design is generated up front for --max-samples and
evaluated incrementally, appending each completed row to the CSV, so the run
can be stopped and restarted and will continue where it left off. Analysis
reads the CSV, so it can be run at any point.

USAGE
    python uncertainty_analysis.py --samples 128           # sample, then analyse
    python uncertainty_analysis.py --samples 512           # continue to 512
    python uncertainty_analysis.py --analyse-only          # re-analyse, no sampling
    python uncertainty_analysis.py --samples 2048 --max-samples 2048
"""
import argparse
import os
import sys
import time

import numpy as np

# --------------------------------------------------------------------------
# path shim: this file lives beside sens_common.py, but topo_loader.py lives
# in a sibling folder (Networks/). Search this directory, the parent, and the
# parent's immediate subdirectories so the script runs from anywhere.
# --------------------------------------------------------------------------
def _add_module_paths(*modules):
    here = os.path.dirname(os.path.abspath(__file__))
    parent = os.path.dirname(here)
    candidates = [here, parent]
    if os.path.isdir(parent):
        candidates += [os.path.join(parent, d) for d in sorted(os.listdir(parent))
                       if os.path.isdir(os.path.join(parent, d))]
    for mod in modules:
        for d in candidates:
            if os.path.exists(os.path.join(d, mod + ".py")) and d not in sys.path:
                sys.path.insert(0, d)
                break

_add_module_paths("sens_common", "topo_loader")

# --------------------------------------------------------------------------
# priors: the scientific input to this study, kept visible rather than buried
# --------------------------------------------------------------------------
# name          -> (low, high, scale, applies_to)
PRIORS = {
    "dv_eta":  (0.65,   0.93,   "lin", "DV"),   # deployed SNSPD efficiency
    "dv_dark": (1.0,    1000.0, "log", "DV"),   # dark counts, cps, 3 decades
    "dv_qber": (0.005,  0.021,  "lin", "DV"),   # intrinsic error rate
    "dv_beta": (0.90,   0.96,   "lin", "DV"),   # reconciliation efficiency
    "cv_eta":  (0.60,   0.72,   "lin", "CV"),   # deployed coherent receiver
    "cv_vel":  (0.015,  0.11,   "lin", "CV"),   # electronic noise, SNU
    "cv_xib":  (0.0005, 0.002,  "lin", "CV"),   # Bob-side excess noise, SNU
    "cv_beta": (0.90,   0.96,   "lin", "CV"),   # reconciliation efficiency
    "alpha":   (0.18,   0.21,   "lin", "both"), # deployed fibre attenuation
}
OUTPUTS = ["crossover_km", "dv_reach_km", "cv_reach_km", "span_ratio"]
TOPOLOGIES = ["TATANID", "SAGO", "GERMANY50", "CESNET", "LAYER42", "RNPBRAZIL"]


# --------------------------------------------------------------------------
# forward model: one parameter draw -> the reported quantities
# --------------------------------------------------------------------------
def _bisect(f, thr, lo, hi, iters):
    """Largest x in [lo,hi] with f(x) > thr, or nan if f(lo) already below."""
    if f(lo) <= thr:
        return np.nan
    if f(hi) > thr:
        return hi
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if f(mid) > thr:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def evaluate(p, floor_bps=1e3, clock_hz=1e9, iters=14):
    """Run the two engines at one parameter draw. Returns a dict of outputs."""
    from sens_common import dv_rate, cv_rate, DARK_CPS_TO_PERGATE
    thr = floor_bps / clock_hz          # floor in bits per channel use

    DV = lambda L: dv_rate(L_km=L, eta=p["dv_eta"], alpha=p["alpha"],
                           dark_pergate=p["dv_dark"] * DARK_CPS_TO_PERGATE,
                           qber=p["dv_qber"], beta=p["dv_beta"])
    CV = lambda L: cv_rate(L_km=L, eta=p["cv_eta"], alpha=p["alpha"],
                           vel=p["cv_vel"], xi_b=p["cv_xib"],
                           beta=p["cv_beta"], detection="heterodyne")

    # crossover: the distance at which CV stops leading. Bracketed rather than
    # assumed, so a draw in which CV never leads (or always leads out to its
    # cutoff) is recorded as nan instead of returning a bracket endpoint.
    diff = lambda L: CV(L) - DV(L)
    lo, hi = 2.0, 120.0
    if diff(lo) <= 0 or diff(hi) > 0:
        crossover = np.nan
    else:
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            if diff(mid) > 0:
                lo = mid
            else:
                hi = mid
        crossover = 0.5 * (lo + hi)

    dv_reach = _bisect(DV, thr, 1.0, 400.0, iters)
    cv_reach = _bisect(CV, thr, 1.0, 300.0, iters)
    return dict(crossover_km=crossover, dv_reach_km=dv_reach,
                cv_reach_km=cv_reach,
                span_ratio=dv_reach / cv_reach if cv_reach else np.nan)


# --------------------------------------------------------------------------
# design
# --------------------------------------------------------------------------
def build_design(n, seed):
    """Scrambled Latin hypercube mapped onto the priors."""
    from scipy.stats import qmc
    names = list(PRIORS)
    u = qmc.LatinHypercube(d=len(names), scramble=True, seed=seed).random(n)
    cols = {}
    for j, k in enumerate(names):
        lo, hi, scale, _ = PRIORS[k]
        cols[k] = (10 ** (u[:, j] * (np.log10(hi) - np.log10(lo)) + np.log10(lo))
                   if scale == "log" else lo + u[:, j] * (hi - lo))
    return names, cols


def sample(csv_path, n_target, max_n, seed, floor_bps):
    names, cols = build_design(max_n, seed)
    header = names + OUTPUTS
    done = 0
    if os.path.exists(csv_path):
        with open(csv_path) as fh:
            done = max(sum(1 for _ in fh) - 1, 0)
        print(f"resuming: {done} samples already in {os.path.basename(csv_path)}")
    else:
        with open(csv_path, "w") as fh:
            fh.write(",".join(header) + "\n")
    if done >= n_target:
        print(f"already have {done} >= {n_target}; nothing to do")
        return
    t0 = time.time()
    for i in range(done, min(n_target, max_n)):
        p = {k: cols[k][i] for k in names}
        out = evaluate(p, floor_bps=floor_bps)
        with open(csv_path, "a") as fh:
            fh.write(",".join(f"{p[k]:.10g}" for k in names) + "," +
                     ",".join(f"{out[k]:.10g}" for k in OUTPUTS) + "\n")
        if (i + 1) % 5 == 0 or i + 1 == n_target:
            el = time.time() - t0
            rate = el / (i + 1 - done)
            print(f"  {i+1}/{n_target}  {rate:.2f} s/sample  "
                  f"eta {(n_target-i-1)*rate/60:.1f} min", flush=True)


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def boot_ci(x, stat, reps=4000, alpha=0.05, seed=0):
    """Percentile bootstrap CI for any statistic of a 1-D sample."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(reps, len(x)))
    vals = np.array([stat(x[i]) for i in idx])
    return np.percentile(vals, [100 * alpha / 2, 100 * (1 - alpha / 2)])


def srrc(X, y, names):
    """Standardised rank regression coefficients, with the rank model's R^2.

    Ranks are used because the responses are monotone in the inputs but not
    linear; standardising makes the coefficients comparable across inputs of
    different units. R^2 says how much of the rank variance the additive model
    explains: if it is low, interactions matter and these numbers understate.
    """
    from scipy.stats import rankdata
    R = np.column_stack([rankdata(X[:, j]) for j in range(X.shape[1])])
    r = rankdata(y)
    Rz = (R - R.mean(0)) / R.std(0)
    rz = (r - r.mean()) / r.std()
    beta, *_ = np.linalg.lstsq(Rz, rz, rcond=None)
    r2 = 1 - np.sum((rz - Rz @ beta) ** 2) / np.sum(rz ** 2)
    return dict(zip(names, beta)), r2


def convergence(x, step=16):
    """Running median and running 90% interval width, to show N is adequate."""
    out = []
    for n in range(step, len(x) + 1, step):
        s = x[:n]
        out.append((n, np.median(s),
                    np.percentile(s, 95) - np.percentile(s, 5)))
    return np.array(out)


def relay_budget(dv_span, cv_span):
    """Propagate sampled spans to trusted-node counts on the real topologies.

    Under the reach criterion a link of length L is cut into ceil(L/span) hops,
    needing ceil(L/span)-1 relays. No engine calls, so this is essentially free
    once the spans are known.
    """
    from topo_loader import load_topology
    res = {}
    for name in TOPOLOGIES:
        L = np.asarray(load_topology(name)["edge_len_km"], dtype=float)
        dv = np.array([np.ceil(L / s).sum() - len(L) for s in dv_span])
        cv = np.array([np.ceil(L / s).sum() - len(L) for s in cv_span])
        res[name] = (dv, cv)
    return res


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def analyse(csv_path, out_dir, seed=0):
    import csv as _csv
    with open(csv_path) as fh:
        rows = list(_csv.DictReader(fh))
    if not rows:
        sys.exit("no samples yet; run with --samples first")
    names = list(PRIORS)
    X = np.array([[float(r[k]) for k in names] for r in rows])
    Y = {k: np.array([float(r[k]) for r in rows]) for k in OUTPUTS}
    n = len(rows)
    L = []
    add = L.append

    add("=" * 78)
    add(f"UNCERTAINTY ANALYSIS — {n} Latin hypercube samples over the deployed bands")
    add("=" * 78)
    add("")
    add("Priors (independent; see module docstring for the independence caveat)")
    for k, (lo, hi, sc, who) in PRIORS.items():
        add(f"   {k:9s} {who:4s} {sc}-uniform  [{lo:g}, {hi:g}]")
    add("")

    # ---- 1. marginals -----------------------------------------------------
    add("-" * 78)
    add("1. MARGINAL DISTRIBUTIONS  (bootstrap 95% CI on each statistic)")
    add("-" * 78)
    add(f"{'quantity':16s} {'n':>4s} {'median':>9s} {'  95% CI':>16s} "
        f"{'p5':>8s} {'p95':>8s} {'CoV':>7s}")
    for k in OUTPUTS:
        v = Y[k][~np.isnan(Y[k])]
        med = np.median(v)
        ci = boot_ci(v, np.median, seed=seed)
        add(f"{k:16s} {len(v):4d} {med:9.2f} [{ci[0]:6.2f},{ci[1]:6.2f}] "
            f"{np.percentile(v,5):8.2f} {np.percentile(v,95):8.2f} "
            f"{v.std()/v.mean():7.3f}")
    nan_x = int(np.isnan(Y['crossover_km']).sum())
    if nan_x:
        add(f"\n   {nan_x} draw(s) had no crossover inside [2,120] km and are excluded")
        add("   from the crossover row only (CV never led, or still led at 120 km).")
    add("")

    # ---- 2. why the ratio is tighter -------------------------------------
    add("-" * 78)
    add("2. VARIANCE DECOMPOSITION — why the span ratio beats the crossover")
    add("-" * 78)
    ld, lc = np.log(Y["dv_reach_km"]), np.log(Y["cv_reach_km"])
    vd, vc, cov = ld.var(), lc.var(), np.cov(ld, lc)[0, 1]
    add(f"   Var(ln DV reach)            {vd:.5f}")
    add(f"   Var(ln CV reach)            {vc:.5f}")
    add(f"   2*Cov(ln DV, ln CV)         {2*cov:+.5f}   (corr {np.corrcoef(ld,lc)[0,1]:+.3f})")
    add(f"   Var(ln ratio) = sum         {vd+vc-2*cov:.5f}")
    add(f"   vs sum without cancellation {vd+vc:.5f}"
        f"   -> cancellation removes {100*(1-(vd+vc-2*cov)/(vd+vc)):.1f}%")
    cvx = Y["crossover_km"][~np.isnan(Y["crossover_km"])]
    add("")
    add(f"   coefficient of variation:  crossover {cvx.std()/cvx.mean():.3f}   "
        f"span ratio {Y['span_ratio'].std()/Y['span_ratio'].mean():.3f}   "
        f"-> ratio is {(cvx.std()/cvx.mean())/(Y['span_ratio'].std()/Y['span_ratio'].mean()):.1f}x tighter")
    add("")
    add("   A ratio cancels any input that moves both protocols the same way; a")
    add("   crossover is a difference of two curves and inherits both variances.")
    add("")

    # ---- 3. global sensitivity -------------------------------------------
    add("-" * 78)
    add("3. GLOBAL SENSITIVITY — standardised rank regression coefficients")
    add("-" * 78)
    for k in OUTPUTS:
        m = ~np.isnan(Y[k])
        coef, r2 = srrc(X[m], Y[k][m], names)
        add(f"\n   {k}   (rank model R^2 = {r2:.3f}"
            f"{'  — LOW, treat with caution' if r2 < 0.7 else ''})")
        for nm, b in sorted(coef.items(), key=lambda t: -abs(t[1])):
            if abs(b) < 0.05:
                continue
            bar = "#" * int(round(abs(b) * 40))
            add(f"      {nm:9s} {b:+.3f}  {bar}")
    add("")

    # ---- 4. convergence ---------------------------------------------------
    add("-" * 78)
    add("4. CONVERGENCE — is N adequate?")
    add("-" * 78)
    add(f"{'quantity':16s} {'N/4':>18s} {'N/2':>18s} {'N':>18s}")
    for k in OUTPUTS:
        v = Y[k][~np.isnan(Y[k])]
        c = convergence(v, step=max(len(v) // 4, 1))
        if len(c) >= 4:
            cells = [f"{c[i][1]:8.2f} +/-{c[i][2]/2:6.2f}" for i in (0, 1, 3)]
            add(f"{k:16s} " + " ".join(f"{x:>18s}" for x in cells))
    add("\n   Stable median and interval width across the last doubling indicates")
    add("   the sample is adequate; drift indicates it is not.")
    add("")

    # ---- 5. relay budget --------------------------------------------------
    add("-" * 78)
    add("5. PROPAGATION TO THE TRUSTED-NODE BUDGET (reach criterion)")
    add("-" * 78)
    try:
        rb = relay_budget(Y["dv_reach_km"], Y["cv_reach_km"])
        add(f"{'topology':11s} {'DV relays':>18s} {'CV relays':>18s} {'CV/DV ratio':>20s}")
        add(f"{'':11s} {'median [p5,p95]':>18s} {'median [p5,p95]':>18s} {'median [p5,p95]':>20s}")
        tot_d = np.zeros(n); tot_c = np.zeros(n)
        for t in TOPOLOGIES:
            d, c = rb[t]; tot_d += d; tot_c += c
            rr = np.where(d > 0, c / np.maximum(d, 1), np.nan)
            rs = ("undefined (DV needs 0)" if np.all(np.isnan(rr))
                  else f"{np.nanmedian(rr):5.2f} [{np.nanpercentile(rr,5):.2f},{np.nanpercentile(rr,95):.2f}]")
            add(f"{t:11s} {np.median(d):6.0f} [{np.percentile(d,5):.0f},{np.percentile(d,95):.0f}]"
                f"{'':4s} {np.median(c):6.0f} [{np.percentile(c,5):.0f},{np.percentile(c,95):.0f}]"
                f"{'':4s} {rs:>20s}")
        tr = tot_c / tot_d
        add("")
        add(f"{'TOTAL':11s} {np.median(tot_d):6.0f} [{np.percentile(tot_d,5):.0f},{np.percentile(tot_d,95):.0f}]"
            f"{'':4s} {np.median(tot_c):6.0f} [{np.percentile(tot_c,5):.0f},{np.percentile(tot_c,95):.0f}]"
            f"{'':4s} {np.median(tr):5.2f} [{np.percentile(tr,5):.2f},{np.percentile(tr,95):.2f}]")
        add("")
        add("   This is the headline of the work carried through the full joint")
        add("   parameter uncertainty rather than quoted at a single baseline.")
    except ImportError as e:                                  # noqa: BLE001
        add(f"   [SKIPPED: {e}]")
        add("   topo_loader.py was not found. This section needs the real edge")
        add("   lengths. Either run from a directory where Networks/ is a sibling,")
        add("   or: PYTHONPATH=/path/to/Networks python uncertainty_analysis.py --analyse-only")
    except Exception as e:                                    # noqa: BLE001
        add(f"   [skipped: {e}]")
    add("")
    add("=" * 78)

    txt = "\n".join(L)
    print(txt)
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "uncertainty_summary.txt"), "w") as fh:
        fh.write(txt + "\n")
    _figure(X, Y, names, out_dir, seed)
    print(f"\nWrote summary and figure to {out_dir}")


def _figure(X, Y, names, out_dir, seed):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.size": 10,
                         "savefig.dpi": 300, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))

    # (a) crossover vs (b) span ratio, the central comparison
    for a, k, col, lab in ((ax[0, 0], "crossover_km", "#c0392b", "Crossover (km)"),
                           (ax[0, 1], "span_ratio", "#1f4e9c", "Span ratio DV/CV")):
        v = Y[k][~np.isnan(Y[k])]
        a.hist(v, bins=28, color=col, alpha=.75, edgecolor="white", linewidth=.5)
        for q, ls in ((5, ":"), (50, "-"), (95, ":")):
            a.axvline(np.percentile(v, q), color="k", ls=ls, lw=1.2)
        a.set_xlabel(lab); a.set_ylabel("samples")
        a.set_title(f"{lab}: median {np.median(v):.2f}, "
                    f"5–95% [{np.percentile(v,5):.2f}, {np.percentile(v,95):.2f}]\n"
                    f"CoV = {v.std()/v.mean():.3f}", fontsize=10.5)
        a.grid(alpha=.25)

    # (c) tornado for the span ratio
    coef, r2 = srrc(X, Y["span_ratio"], names)
    items = sorted(coef.items(), key=lambda t: abs(t[1]))
    ax[1, 0].barh([i[0] for i in items], [i[1] for i in items],
                  color=["#c0392b" if v < 0 else "#1f4e9c" for _, v in items])
    ax[1, 0].axvline(0, color="k", lw=.8)
    ax[1, 0].set_xlabel("standardised rank regression coefficient")
    ax[1, 0].set_title(f"What drives the span ratio  ($R^2$ = {r2:.2f})", fontsize=10.5)
    ax[1, 0].grid(alpha=.25, axis="x")

    # (d) convergence of both headline quantities
    for k, col in (("crossover_km", "#c0392b"), ("span_ratio", "#1f4e9c")):
        v = Y[k][~np.isnan(Y[k])]
        c = convergence(v, step=max(len(v) // 20, 4))
        ax[1, 1].plot(c[:, 0], c[:, 2] / np.median(v), "-o", ms=3, color=col,
                      label=f"{k} (90% width / median)")
    ax[1, 1].set_xlabel("samples"); ax[1, 1].set_ylabel("relative 90% interval width")
    ax[1, 1].set_title("Convergence", fontsize=10.5)
    ax[1, 1].legend(fontsize=8.5); ax[1, 1].grid(alpha=.25)

    fig.suptitle("Joint parameter uncertainty over the deployed bands",
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, .96])
    fig.savefig(os.path.join(out_dir, "uncertainty_analysis.png"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=128,
                    help="evaluate up to this many samples (resumable)")
    ap.add_argument("--max-samples", type=int, default=4096,
                    help="size of the fixed LHS design; keep constant across runs")
    ap.add_argument("--seed", type=int, default=20260805)
    ap.add_argument("--floor-bps", type=float, default=1e3,
                    help="usability floor in bit/s at the matched clock")
    ap.add_argument("--output-dir", type=str, default="uncertainty")
    ap.add_argument("--analyse-only", action="store_true")
    a = ap.parse_args()

    os.makedirs(a.output_dir, exist_ok=True)
    csv_path = os.path.join(a.output_dir, "uncertainty_samples.csv")
    if not a.analyse_only:
        sample(csv_path, a.samples, a.max_samples, a.seed, a.floor_bps)
    analyse(csv_path, a.output_dir, seed=a.seed)


if __name__ == "__main__":
    main()