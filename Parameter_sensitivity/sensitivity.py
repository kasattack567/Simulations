"""
sensitivity_grid.py — the six parameter sweeps as ONE six-panel figure at
matched clock (DV and CV both at 1 GHz).

s1..s6 each show two panels (matched | deployed clock); twelve don't fit one
grid, so this recomputes the matched-clock reading from the same engines and
baseline in sens_common (so the figures cannot disagree).

Panel (a) uses RANGES['dv_eta'] (starts at 0.35), while s1_detector_eff.py
hardcodes 0.50-1.00; the DV band starts at 0.411, so its lower edge sits off
the standalone plot but is visible here.

Panels (a), (c), (f) share an x-axis, so crossings between curves are real.
(b), (d), (e) use twin x-axes with no common unit — crossings there are
axis-alignment artefacts. (f) is CV-only (no DV channel-noise floor).

Every plotted value is printed as the figure is built.

    python sensitivity_grid.py                 # -> figures/sensitivity_grid.png
    python sensitivity_grid.py --distance 50   # sweep at another distance
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sens_common import (dv_rate, cv_rate, band_label, band_scaled,
                         RANGES, BAND, BAND_DV, BAND_CV, BAND_SH,
                         DV_COLOR, CV_COLOR, DARK_CPS_TO_PERGATE,
                         CV_XI_B, L_KM)

HOM_COLOR = "#7d2d8c"
CLOCK_HZ = 1e9          # matched clock, both protocols

plt.rcParams.update({
    "font.family": "serif", "font.size": 10.5,
    "axes.titlesize": 11.5, "axes.labelsize": 10.5,
    "legend.fontsize": 8.2, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.linewidth": 0.9, "lines.linewidth": 1.9, "lines.markersize": 3.2,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})


def _bps(y):
    """bits/channel use -> bit/s at the matched clock, zeros masked.

    A zero cannot be drawn on a log axis, and the segment running into it reads
    as a real value; NaN breaks the line where the rate actually vanishes.
    """
    y = np.asarray(y, dtype=float) * CLOCK_HZ
    return np.where(y > 0, y, np.nan)


def _report(panel, x, curves, band, xname):
    """Print exactly what a panel draws: endpoints, span, and band-only change.

    `band` is the shaded deployed range for these curves; the change across it
    is the figure the report quotes, which is usually far smaller than the
    change across the whole swept axis.
    """
    print("\n  %s" % panel)
    print("    axis %s: %.4g to %.4g  (%d points)   deployed band %.4g to %.4g"
          % (xname, x.min(), x.max(), len(x), band[0], band[1]))
    for name, y in curves:
        y = np.asarray(y, dtype=float)
        ok = np.isfinite(y)
        if not ok.any():
            print("      %-16s all zero" % name)
            continue
        f, l = y[ok][0], y[ok][-1]
        line = "      %-16s %.4e -> %.4e   (%.3gx over axis" % (name, f, l, l / f)
        a, b = np.interp(band, x, y)
        if np.isfinite(a) and np.isfinite(b):
            line += ", %+.1f%% across band" % (100 * (b / a - 1))
        line += ")"
        if not ok.all():
            line += "   [zero at %s = %.4g]" % (xname, x[~ok][0])
        print(line)


def _shade(ax, band, color, label):
    lo, hi = band
    h = ax.axvspan(lo, hi, alpha=0.16, color=color, zorder=0, label=label)
    ax.axvline(lo, color=color, ls=":", lw=1.0, alpha=0.8, zorder=1)
    ax.axvline(hi, color=color, ls=":", lw=1.0, alpha=0.8, zorder=1)
    return h


def _finish(ax, xlabel, title, handles, logx=False, loc="best"):
    if logx:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title, fontsize=11.5)
    ax.grid(True, which="both", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(handles=handles, loc=loc, framealpha=0.92, edgecolor="0.7")


def _twin_style(ax, axt, cv_label, title):
    """Colour-code a twin-axis panel and lift its title clear of the top label."""
    axt.set_yscale("log")
    axt.set_xlabel(cv_label, color=CV_COLOR)
    axt.tick_params(axis="x", colors=CV_COLOR)
    ax.tick_params(axis="x", colors=DV_COLOR)
    ax.spines["bottom"].set_color(DV_COLOR)
    axt.spines["top"].set_color(CV_COLOR)
    ax.xaxis.label.set_color(DV_COLOR)
    axt.set_title(title, pad=34)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--distance", type=float, default=L_KM)
    p.add_argument("--output-dir", type=str, default="figures")
    p.add_argument("--save", type=str, default="sensitivity_grid.png")
    a = p.parse_args()
    L = a.distance

    print("=" * 78)
    print("SENSITIVITY GRID — matched clock, DV and CV both %.0e Hz" % CLOCK_HZ)
    print("L = %.3f km.  All rates in bit/s, exactly as plotted." % L)
    print("=" * 78)

    fig, axes = plt.subplots(3, 2, figsize=(13.5, 15.5))
    (axA, axB), (axC, axD), (axE, axF) = axes

    # ---------------- (a) receiver efficiency — shared axis ----------------
    # NOT detector efficiency. Both engines consume Bob's WHOLE-RECEIVER
    # transmission: TNO's efficiency_party and qosst-skr's eta are the same
    # physical quantity, so the axis is shared and crossings are meaningful.
    # The DV band is the detector SDE reduced by the receiver optical budget;
    # the CV band is already a whole-receiver figure in its sources.
    eta = RANGES["dv_eta"]
    dv = _bps([dv_rate(L_km=L, eta=e) for e in eta])
    het = _bps([cv_rate(L_km=L, eta=e, detection="heterodyne") for e in eta])
    hom = _bps([cv_rate(L_km=L, eta=e, detection="homodyne") for e in eta])
    h = [_shade(axA, BAND["dv_eta"], BAND_DV, band_label("dv_eta")),
         _shade(axA, BAND["cv_eta"], BAND_CV, band_label("cv_eta"))]
    h += [axA.plot(eta, dv, "-o", color=DV_COLOR, label="DV — decoy BB84")[0],
          axA.plot(eta, het, "-s", color=CV_COLOR, label="CV — heterodyne")[0],
          axA.plot(eta, hom, "--D", color=HOM_COLOR, label="CV — homodyne")[0]]
    _finish(axA, "Bob total receiver efficiency  $\\eta_{B}$",
            "(a) Receiver efficiency", h)
    _report("(a) RECEIVER EFFICIENCY — over the DV band", eta,
            [("DV", dv), ("CV het", het), ("CV hom", hom)],
            BAND["dv_eta"], "eta_B")
    _report("(a) RECEIVER EFFICIENCY — over the CV band", eta,
            [("CV het", het)], BAND["cv_eta"], "eta_B")
    dlo, dhi = BAND["dv_eta"]
    clo, chi = BAND["cv_eta"]
    lo, hi = max(dlo, clo), min(dhi, chi)
    if lo >= hi:
        rel = "the two bands do NOT overlap"
    elif clo >= dlo and chi <= dhi:
        rel = "the CV band lies entirely INSIDE the DV band"
    elif dlo >= clo and dhi <= chi:
        rel = "the DV band lies entirely INSIDE the CV band"
    else:
        rel = ("they overlap only partially; CV extends higher (%.3g vs %.3g)"
               % (chi, dhi) if chi > dhi else
               "they overlap only partially; DV extends higher (%.3g vs %.3g)"
               % (dhi, chi))
    print("      DV band %.3g-%.3g, CV band %.3g-%.3g: %s"
          % (dlo, dhi, clo, chi, rel))
    if lo < hi:
        print("      common region %.3g to %.3g" % (lo, hi))

    # ---------------- (b) reconciliation — TWIN axis ----------------
    f_ec, beta = RANGES["dv_f_ec"], RANGES["beta"]
    axBt = axB.twiny()
    dv = _bps([dv_rate(L_km=L, f_ec=x) for x in f_ec])
    het = _bps([cv_rate(L_km=L, beta=x, detection="heterodyne") for x in beta])
    hom = _bps([cv_rate(L_km=L, beta=x, detection="homodyne") for x in beta])
    h = [_shade(axB, BAND["dv_f_ec"], BAND_DV, band_label("dv_f_ec")),
         _shade(axBt, BAND["beta"], BAND_CV, band_label("beta"))]
    h += [axB.plot(f_ec, dv, "-o", color=DV_COLOR, label="DV — decoy BB84")[0],
          axBt.plot(beta, het, "-s", color=CV_COLOR, label="CV — heterodyne")[0],
          axBt.plot(beta, hom, "--D", color=HOM_COLOR,
                    label="CV — homodyne")[0]]
    _twin_style(axB, axBt, "CV reconciliation efficiency  $\\beta$",
                "(b) Reconciliation")
    _finish(axB, "DV error-correction inefficiency  $f$", "", h,
            loc="lower right")
    _report("(b) RECONCILIATION — DV axis (f_EC, larger is worse)", f_ec,
            [("DV", dv)], BAND["dv_f_ec"], "f_EC")
    _report("(b) RECONCILIATION — CV axis (beta, larger is better)", beta,
            [("CV het", het), ("CV hom", hom)], BAND["beta"], "beta")

    # ---------------- (c) fibre attenuation — shared axis ----------------
    al = RANGES["alpha"]
    dv = _bps([dv_rate(L_km=L, alpha=x) for x in al])
    het = _bps([cv_rate(L_km=L, alpha=x, detection="heterodyne") for x in al])
    hom = _bps([cv_rate(L_km=L, alpha=x, detection="homodyne") for x in al])
    h = [_shade(axC, BAND["alpha"], BAND_SH, band_label("alpha"))]
    h += [axC.plot(al, dv, "-o", color=DV_COLOR, label="DV — decoy BB84")[0],
          axC.plot(al, het, "-s", color=CV_COLOR, label="CV — heterodyne")[0],
          axC.plot(al, hom, "--D", color=HOM_COLOR, label="CV — homodyne")[0]]
    _finish(axC, "Fibre attenuation  $\\alpha$  (dB/km)",
            "(c) Fibre attenuation", h)
    _report("(c) FIBRE ATTENUATION", al,
            [("DV", dv), ("CV het", het), ("CV hom", hom)],
            BAND["alpha"], "alpha")
    print("      CV/DV ratio   " + "   ".join(
        "%.3g dB/km: %.2f" % (x, cv_rate(L_km=L, alpha=x,
                                         detection="heterodyne")
                              / dv_rate(L_km=L, alpha=x))
        for x in (0.25, 0.20, 0.091, 0.01)))

    # ---------------- (d) detector noise — TWIN axis ----------------
    dk, vel = RANGES["dv_dark"], RANGES["cv_vel"]
    axDt = axD.twiny()
    dv = _bps([dv_rate(L_km=L, dark_pergate=c * DARK_CPS_TO_PERGATE)
               for c in dk])
    het = _bps([cv_rate(L_km=L, vel=v, detection="heterodyne") for v in vel])
    hom = _bps([cv_rate(L_km=L, vel=v, detection="homodyne") for v in vel])
    h = [_shade(axD, BAND["dv_dark"], BAND_DV, band_label("dv_dark")),
         _shade(axDt, BAND["cv_vel"], BAND_CV, band_label("cv_vel"))]
    h += [axD.plot(dk, dv, "-o", color=DV_COLOR, label="DV — decoy BB84")[0],
          axDt.plot(vel, het, "-s", color=CV_COLOR, label="CV — heterodyne")[0],
          axDt.plot(vel, hom, "--D", color=HOM_COLOR,
                    label="CV — homodyne")[0]]
    axD.set_xscale("log")
    _twin_style(axD, axDt, "CV electronic noise  $v_{el}$  (SNU)",
                "(d) Detector noise")
    _finish(axD, "DV dark-count rate  (cps)", "", h, logx=True,
            loc="lower left")
    _report("(d) DETECTOR NOISE — DV axis", dk, [("DV", dv)],
            BAND["dv_dark"], "dark cps")
    _report("(d) DETECTOR NOISE — CV axis", vel,
            [("CV het", het), ("CV hom", hom)], BAND["cv_vel"], "v_el")

    # ---------------- (e) channel noise — TWIN axis ----------------
    qb, xb = RANGES["dv_qber"], RANGES["cv_xi_b"]
    axEt = axE.twiny()
    dv = _bps([dv_rate(L_km=L, qber=q) for q in qb])
    het = _bps([cv_rate(L_km=L, xi_b=x, detection="heterodyne") for x in xb])
    hom = _bps([cv_rate(L_km=L, xi_b=x, detection="homodyne") for x in xb])
    h = [_shade(axE, band_scaled("dv_qber"), BAND_DV, band_label("dv_qber")),
         _shade(axEt, band_scaled("cv_xi_b"), BAND_CV, band_label("cv_xi_b"))]
    h += [axE.plot(qb * 100, dv, "-o", color=DV_COLOR,
                   label="DV — decoy BB84")[0],
          axEt.plot(xb * 1e3, het, "-s", color=CV_COLOR,
                    label="CV — heterodyne")[0],
          axEt.plot(xb * 1e3, hom, "--D", color=HOM_COLOR,
                    label="CV — homodyne")[0]]
    _twin_style(axE, axEt,
                "CV Bob-side excess noise  $\\epsilon_b$  ($10^{-3}$ SNU)",
                "(e) Channel noise")
    _finish(axE, "DV QBER  (%)", "", h, loc="lower left")
    _report("(e) CHANNEL NOISE — DV axis", qb * 100, [("DV", dv)],
            band_scaled("dv_qber"), "QBER %")
    _report("(e) CHANNEL NOISE — CV axis", xb * 1e3,
            [("CV het", het), ("CV hom", hom)], band_scaled("cv_xi_b"),
            "eps_b (1e-3)")

    # ---------------- (f) CV excess-noise floor — CV only ----------------
    xa = RANGES["cv_xi_al"]
    het = _bps([cv_rate(L_km=L, xi_al=x, detection="heterodyne") for x in xa])
    hom = _bps([cv_rate(L_km=L, xi_al=x, detection="homodyne") for x in xa])
    h = [_shade(axF, BAND["cv_xi_al"], BAND_SH, band_label("cv_xi_al"))]
    h += [axF.plot(xa, het, "-s", color=CV_COLOR, label="CV — heterodyne")[0],
          axF.plot(xa, hom, "--D", color=HOM_COLOR, label="CV — homodyne")[0]]
    _finish(axF, "CV channel noise floor  $\\epsilon_a + \\epsilon_l$  (SNU)",
            "(f) CV excess-noise floor  ($\\epsilon_b = %.0e$ SNU)" % CV_XI_B,
            h, loc="upper right")
    _report("(f) CV EXCESS-NOISE FLOOR", xa,
            [("CV het", het), ("CV hom", hom)], BAND["cv_xi_al"],
            "eps_a+eps_l")

    for ax in (axA, axC, axE):
        ax.set_ylabel("Secret key rate  (bit / s)")

    fig.suptitle("Parameter sensitivity of DV and CV QKD at the crossover "
                 "($L = %.0f$ km, matched clock: DV and CV both 1 GHz)" % L,
                 fontsize=13.5, weight="bold", y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.subplots_adjust(hspace=0.50, wspace=0.22)

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          a.output_dir)
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, a.save)
    fig.savefig(path)
    plt.close(fig)

    # ---------------- closing anchors ----------------
    def _zero(fn, lo, hi):
        for _ in range(60):
            m = 0.5 * (lo + hi)
            if fn(m) > 0:
                lo = m
            else:
                hi = m
        return lo

    print("\n" + "=" * 78)
    print("BASELINE ANCHORS at this distance (bit/s, matched clock)")
    print("=" * 78)
    print("    DV                 %.4e" % (dv_rate(L_km=L) * CLOCK_HZ))
    print("    CV heterodyne      %.4e"
          % (cv_rate(L_km=L, detection="heterodyne") * CLOCK_HZ))
    print("    CV homodyne        %.4e"
          % (cv_rate(L_km=L, detection="homodyne") * CLOCK_HZ))
    print("    hom/het            %+.2f %%"
          % (100 * (cv_rate(L_km=L, detection="homodyne")
                    / cv_rate(L_km=L, detection="heterodyne") - 1)))
    print("    CV/DV              %.3f" % (cv_rate(L_km=L,
                                                   detection="heterodyne")
                                           / dv_rate(L_km=L)))
    print("    DV zero at QBER    %.2f %%"
          % (100 * _zero(lambda q: dv_rate(L_km=L, qber=q), 0.02, 0.15)))
    print("    CV zero at eps_b   %.3e SNU"
          % _zero(lambda x: cv_rate(L_km=L, xi_b=x, detection="heterodyne"),
                  1e-4, 5e-2))
    print("=" * 78)
    print("Saved: %s" % path)


if __name__ == "__main__":
    main()