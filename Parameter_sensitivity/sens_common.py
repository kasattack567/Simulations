"""
Shared engine wrappers, locked baseline, and plot helpers for the DV vs CV
parameter sensitivity analysis (fixed distance, near the crossover).

Imported by the per-parameter scripts s1..s5; not run directly. Keeping the
baseline and engine calls in ONE place guarantees they stay identical to
cv_hetro.py / param.md (no per-script drift).

METHOD NOTE: qosst-skr and TNO are DETERMINISTIC analytic formulas — one input,
one output. There is no run-to-run variance, so these plots carry NO error bars
(unlike a Monte-Carlo simulator such as NetSquid, where the spread is real).

DV RECONCILIATION: TNO's asymptotic BB84 hardcodes error correction at the
Shannon limit (beta=1). DV_BETA<1 is applied post-hoc by scaling the EC cost,
exactly as in cv_hetro.py.

CV EXCESS NOISE — Wang et al., Opt. Express 27, 13372 (2019), Sec. 5, Eq. 12:
    xi_r(T) = (eps_a + eps_l) + eps_b / (eta * T)
referred to the channel input. eps_b (Bob-side measurement noise) arises AFTER
the channel and so does not attenuate; referred to the input it is amplified by
1/(eta*T). This replaces the earlier single-parameter xi_bob/(T*eta) form, which
had no constant channel floor and an eps_b ~30x too large (it implied 2.5 SNU
input-referred at 100 km, against ~0.098 SNU measured in Wang Table 2).

CONSEQUENCE FOR THESE SWEEPS: eta and alpha now enter TWICE for CV — once through
detection, once through xi_r. That coupling is physical and intended, but it means
the CV curves in s1 and s3 are steeper than under the old constant-xi model.

UNITS — bit/s, TWO PANELS PER FIGURE. The engines return bits per channel use;
every figure now multiplies by a repetition rate and is plotted in bit/s as a
side-by-side pair sharing one y-axis:
    (a) matched clock   — DV 1 GHz, CV 1 GHz  (protocol physics, like-for-like)
    (b) deployed clock  — DV 1 GHz, CV 100 MHz (what the hardware delivers)
The DV curve is identical in both panels; CV drops by exactly one decade in (b).
Only the plot helpers changed — cv_rate/dv_rate still return bits/channel use, so
nothing downstream of them (crossover distances, network suite, cost model) is
affected unless it calls these plotters. See CLOCK RATES below.
"""
import os
import textwrap
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

from qosst_skr.gaussian_trusted_heterodyne_asymptotic import (
    GaussianTrustedHeterodyneAsymptotic)
from qosst_skr.gaussian_trusted_homodyne_asymptotic import (
    GaussianTrustedHomodyneAsymptotic)
from tno.quantum.communication.qkd_key_rate.quantum import standard_detector
from tno.quantum.communication.qkd_key_rate.quantum.bb84 import (
    BB84FullyAsymptoticKeyRateEstimate, compute_gain_and_error_rate)


# ============================================================
# LOCKED BASELINE (must match cv_hetro.py / param.md section 0)
# ============================================================
# Sensitivity test distance. Chosen to sit near the DV/CV crossover so that both
# protocols are competitive and the sweeps are informative. Set to the exact
# DV/CV crossover at the baseline beta = 0.95, located by bisection (49.891 km).
# NB the crossover is beta-dependent: 27.2 km at beta = 0.90, 56.8 km at 0.96,
# so this is the crossover at baseline, not a universal one. The report rounds
# it to 50 km. Override per-run with --distance.
L_KM = 49.891
ALPHA_DB_KM = 0.20

DV_EFFICIENCY = 0.65
DV_DARK_COUNT = 1e-7       # per gate (100 cps @ ~1 ns / 1 GHz)
DV_QBER       = 0.005
DV_BETA       = 0.95       # matched to CV; applied post-hoc (engine is ideal)

CV_ETA    = 0.60
CV_VEL    = 0.10
CV_BETA   = 0.95

# CV excess noise, Wang Eq. 12 (channel-input referred). Calibrated by Wang
# against their measured prototype (Table 2 / Fig. 8).
CV_XI_AL  = 0.005     # eps_a + eps_l: Alice + fibre channel (SNU)
CV_XI_B   = 0.0005    # eps_b: Bob-side measurement noise (SNU)

# Va capped at 10 SNU per param.md (deployed modulator headroom). Does not affect
# any result beyond ~20 km; below that the uncapped optimum was unphysical.
VA_LO, VA_HI = 1e-2, 10.0

DARK_CPS_TO_PERGATE = 1e-9  # 100 cps @ 1 GHz, ~1 ns window -> 1e-7 per gate


# ============================================================
# CLOCK RATES — bits/channel-use -> bits/second
# ============================================================
# Both engines return an ASYMPTOTIC rate per channel use (DV: per pulse;
# CV: per transmitted coherent state / symbol). Converting to bit/s requires a
# repetition rate, and DV and CV do not run at the same one in deployed hardware:
# modern SNSPD-based DV systems clock at ~1 GHz, while CV-QKD symbol rates are
# ~100 MHz (limited by the DAC/ADC, the DSP chain and the shot-noise-limited
# bandwidth of the balanced receiver, not by the protocol).
#
# Reporting only one of these is a choice with a thumb on the scale, so every
# sensitivity figure is drawn as TWO PANELS:
#   left  — matched clock (1 GHz both): isolates the PROTOCOL physics. Any
#           DV/CV separation here is intrinsic to the protocols, not to the
#           electronics, so this is the like-for-like comparison.
#   right — deployed clock (DV 1 GHz, CV 100 MHz): what the two technologies
#           actually deliver today. CV shifts down by exactly one decade; the DV
#           curve is identical to the left panel by construction.
#
# The two panels SHARE a y-axis so the decade of separation is read directly off
# the figure rather than hidden by independent autoscaling.
DV_CLOCK_HZ  = 1e9    # SNSPD-based DV, ~1 GHz gating/free-running
CV_CLOCK_HZ  = 1e8    # deployed CV symbol rate, ~100 MHz
CV_CLOCK_MATCHED_HZ = DV_CLOCK_HZ   # hypothetical parity, left panel


def _fmt_hz(f):
    return f"{f/1e9:g} GHz" if f >= 1e9 else f"{f/1e6:g} MHz"


# Panel order is left -> right. Each entry gives the (dv, cv) clock in Hz.
PANELS = (
    dict(key='matched',
         dv=DV_CLOCK_HZ, cv=CV_CLOCK_MATCHED_HZ,
         title=f"(a) Matched clock — DV {_fmt_hz(DV_CLOCK_HZ)}, "
               f"CV {_fmt_hz(CV_CLOCK_MATCHED_HZ)}"),
    dict(key='deployed',
         dv=DV_CLOCK_HZ, cv=CV_CLOCK_HZ,
         title=f"(b) Deployed clock — DV {_fmt_hz(DV_CLOCK_HZ)}, "
               f"CV {_fmt_hz(CV_CLOCK_HZ)}"),
)

YLABEL = 'Secret key rate  (bit / s)'
PANEL_FIGSIZE = (13.5, 5.8)
Y_DECADES = 9          # clamp the log y-axis to this many decades below the peak


# ============================================================
# CV ENGINE (identical Va optimiser to cv_hetro.py)
# ============================================================
def T_of_L(L_km, alpha=ALPHA_DB_KM):
    return 10**(-alpha * L_km / 10.0)


def cv_xi_input(T, eta=None, xi_al=None, xi_b=None):
    """Channel-input-referred CV excess noise, Wang et al. 2019 Eq. 12."""
    eta = CV_ETA if eta is None else eta
    xi_al = CV_XI_AL if xi_al is None else xi_al
    xi_b = CV_XI_B if xi_b is None else xi_b
    return xi_al + xi_b / (eta * T)


def cv_rate(L_km=L_KM, eta=None, vel=None, xi_al=None, xi_b=None, beta=None,
            alpha=None, detection="heterodyne"):
    """CV bits/symbol, Va optimised. detection='heterodyne' or 'homodyne'.

    Excess noise follows Wang Eq. 12; sweep it via xi_b (Bob-side, carries the
    distance scaling) and/or xi_al (constant Alice+fibre floor).
    """
    eta = CV_ETA if eta is None else eta
    vel = CV_VEL if vel is None else vel
    beta = CV_BETA if beta is None else beta
    alpha = ALPHA_DB_KM if alpha is None else alpha
    T = T_of_L(L_km, alpha)
    if T < 1e-12:
        return 0.0
    xi = cv_xi_input(T, eta=eta, xi_al=xi_al, xi_b=xi_b)
    skr_fn = (GaussianTrustedHomodyneAsymptotic.skr if detection == "homodyne"
              else GaussianTrustedHeterodyneAsymptotic.skr)

    def neg(va):
        try:
            r = skr_fn(Va=va, T=T, xi=xi, eta=eta, Vel=vel, beta=beta)
            return -r if (r is not None and np.isfinite(r)) else 0.0
        except Exception:
            return 0.0

    grid = np.logspace(np.log10(VA_LO), np.log10(VA_HI), 120)
    vals = np.array([neg(v) for v in grid])
    i = int(np.argmin(vals))
    best = vals[i]
    lo, hi = grid[max(i-1, 0)], grid[min(i+1, len(grid)-1)]
    if hi > lo:
        try:
            res = minimize_scalar(neg, bounds=(lo, hi), method="bounded",
                                  options={"xatol": 1e-4})
            best = min(best, res.fun)
        except Exception:
            pass
    return max(-best, 0.0)


# ============================================================
# DV ENGINE (TNO) + post-hoc reconciliation
# ============================================================
def _h(p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -p*np.log2(p) - (1-p)*np.log2(1-p)


def dv_rate(L_km=L_KM, eta=None, dark_pergate=None, qber=None, beta=None, alpha=None):
    """DV decoy BB84 bits/pulse, asymptotic, with post-hoc reconciliation beta.

    NB: TNO's asymptotic BB84 derives the channel error rate from dark counts and
    POLARIZATION DRIFT, not from the detector's `error_detector` field (which is
    unused in this rate path). Intrinsic/misalignment QBER is therefore encoded as
    a drift angle: polarization_drift = arcsin(sqrt(QBER)) reproduces the target
    QBER to <1e-3 across 0-11% (verified against compute_gain_and_error_rate).
    """
    eta = DV_EFFICIENCY if eta is None else eta
    dark_pergate = DV_DARK_COUNT if dark_pergate is None else dark_pergate
    qber = DV_QBER if qber is None else qber
    beta = DV_BETA if beta is None else beta
    alpha = ALPHA_DB_KM if alpha is None else alpha
    att = alpha * L_km
    drift = float(np.arcsin(np.sqrt(np.clip(qber, 0.0, 0.5))))  # QBER -> drift angle
    det = standard_detector.customise(
        efficiency_party=eta, dark_count_rate=dark_pergate,
        error_detector=qber, polarization_drift=drift)
    try:
        mu_opt, r_ideal = BB84FullyAsymptoticKeyRateEstimate(
            detector=det).optimize_rate(attenuation=att)
    except Exception:
        return 0.0
    r_ideal = max(r_ideal, 0.0)
    if beta >= 1.0 or r_ideal <= 0.0:
        return r_ideal
    try:
        mu = float(np.atleast_1d(mu_opt["mu"])[0])   # optimize_rate returns {"mu": array}
        gain, err = compute_gain_and_error_rate(det, mu, att)
        ec_extra = (1.0/beta - 1.0) * float(np.atleast_1d(gain)[0]) * _h(float(np.atleast_1d(err)[0]))
        return max(r_ideal - ec_extra, 0.0)
    except Exception:
        return r_ideal


# ============================================================
# PLOT STYLE — dissertation-grade defaults
# ============================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.linewidth": 0.9,
    "lines.linewidth": 2.2,
    "lines.markersize": 4,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

DV_COLOR = "#1f4e9c"   # deep blue
CV_COLOR = "#c0392b"   # deep red
BAND_DV  = "#5b8def"   # light blue band
BAND_CV  = "#e8897e"   # light red band
BAND_SH  = "#7cc47f"   # green band (shared parameter)


def _save_or_show(fig, output_dir, fname):
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, fname)
        fig.savefig(path)
        plt.close(fig)
        print(f"Saved: {path}")
    else:
        plt.show()


def _bps(y, clock_hz):
    """bits/channel-use -> bit/s, with non-positive values masked.

    Zeros are masked rather than plotted because a log axis silently drops them
    and the line segment leading into a dropped point is misleading; NaN breaks
    the line cleanly at the point where the rate goes to zero.
    """
    y = np.asarray(y, dtype=float) * float(clock_hz)
    return np.where(y > 0, y, np.nan)


def _autoscale_y(axes, series, decades=Y_DECADES):
    """Clamp the shared log y-axis to `decades` below the global peak.

    Without this a single curve diving toward zero stretches the axis over 20+
    decades and flattens every other curve into a horizontal line.
    """
    vals = np.concatenate([np.asarray(s, dtype=float).ravel() for s in series])
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        return
    top = 10 ** np.ceil(np.log10(vals.max()))
    bot = max(vals.min(), top / 10 ** decades)
    bot = 10 ** np.floor(np.log10(bot))
    for ax in np.atleast_1d(axes).ravel():
        ax.set_ylim(bot, top)


def _panel_axes():
    """Two side-by-side panels sharing a y-axis (matched clock | deployed clock)."""
    return plt.subplots(1, 2, figsize=PANEL_FIGSIZE, sharey=True)


NOTE_WRAP = 135        # caption characters per line before wrapping


def _finish_panels(fig, axes, title, note=None, note_y=None):
    """Label, title and lay out the two panels, then write the caption.

    THE CAPTION MUST BE A FIGURE ARTIST, NOT AN AXES ARTIST. It was previously
    drawn with axes[0].text(1.0, -0.16, ..., transform=axes[0].transAxes), which
    makes a ~200-character string a CHILD of the left axes. tight_layout sizes
    each axes to its tight bounding box INCLUDING its children, so the left axes
    was shrunk until the overhanging caption fitted its half of the figure — the
    panels collapsed to narrow strips with a large gap between them, and
    savefig(bbox='tight') then grew the canvas to the caption's width (13.5 in
    figure -> 14.6 in saved, 18.6 in for the twin-axis figures, whose caption is
    half as long again). fig.text() is laid out after tight_layout has run and
    cannot deform the axes, so the panels keep their allotted width.

    `note_y` is accepted and ignored: the caption is now placed automatically
    below the axes and vertical space is reserved for it in the layout rect.
    """
    axes[0].set_ylabel(YLABEL)
    fig.suptitle(title, y=0.985, fontsize=13.5)
    lines = textwrap.wrap(note, width=NOTE_WRAP) if note else []
    # Reserve one caption line's worth of figure height, plus one per extra line.
    bottom = 0.055 + 0.030 * max(len(lines) - 1, 0)
    fig.tight_layout(rect=(0, bottom, 1, 0.95))
    if lines:
        fig.text(0.5, 0.012, "\n".join(lines), ha='center', va='bottom',
                 fontsize=8.5, style='italic', color='0.35')


CLOCK_NOTE = ("Both panels use the same simulated rates per channel use; only the "
              "repetition rate differs, so the DV curve is identical in (a) and (b) "
              "and CV is shifted down by one decade in (b).")

TWIN_NOTE = ("Note: DV and CV use separate x-axes — curve crossings are not "
             "physically meaningful.")


# ============================================================
# SHARED-AXIS PLOT (parameter identical for both protocols)
# ============================================================
def shared_plot(x, dv, cv, xlabel, title, fname,
                band=None, band_label=None, bands=None, logx=False, output_dir=None,
                labels=('DV — decoy BB84', 'CV — GG02 heterodyne'),
                colors=(DV_COLOR, CV_COLOR), markers=('-o', '-s'),
                kinds=('dv', 'cv'),
                legend_loc='best'):
    """Shared x-axis (parameter identical for both protocols), two clock panels.
    Pass a single band via (band, band_label), or multiple via
    bands=[((lo,hi), color, label), ...]. Crossings here ARE meaningful.

    labels/colors/markers are overridable so this helper can also plot two curves
    that are NOT one-DV-one-CV (e.g. CV heterodyne vs CV homodyne in s6). Leaving
    them at the defaults while passing non-DV data mislabels the figure.

    `kinds` says which clock each curve is scaled by — ('dv','cv') by default.
    Two CV curves must pass kinds=('cv','cv'), otherwise the first would be
    scaled at the DV clock and the deployed panel would be wrong.
    """
    band_list = []
    if bands is not None:
        band_list = bands
    elif band is not None:
        band_list = [(band, BAND_SH, band_label or "Deployed regime")]

    fig, axes = _panel_axes()
    all_series = []
    for ax, panel in zip(axes, PANELS):
        handles = []
        for (lo, hi), color, label in band_list:
            sp = ax.axvspan(lo, hi, alpha=0.16, color=color, zorder=0, label=label)
            ax.axvline(lo, color=color, ls=':', lw=1.0, alpha=0.8, zorder=1)
            ax.axvline(hi, color=color, ls=':', lw=1.0, alpha=0.8, zorder=1)
            handles.append(sp)
        y0 = _bps(dv, panel[kinds[0]])
        y1 = _bps(cv, panel[kinds[1]])
        all_series += [y0, y1]
        h0, = ax.plot(x, y0, markers[0], color=colors[0], label=labels[0])
        h1, = ax.plot(x, y1, markers[1], color=colors[1], label=labels[1])
        handles += [h0, h1]
        if logx:
            ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(xlabel)
        ax.grid(True, which='both', alpha=0.25, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_title(panel['title'], pad=8, fontsize=11.5)
        if ax is axes[0]:
            ax.legend(handles=handles, loc=legend_loc, framealpha=0.92,
                      edgecolor='0.7')
    _autoscale_y(axes, all_series)
    _finish_panels(fig, axes, title, note=CLOCK_NOTE)
    _save_or_show(fig, output_dir, fname)


# ============================================================
# TWIN-AXIS PLOT (DV / CV parameters analogous, not identical)
# Crossing points between the two curves are NOT physically meaningful
# (different x-axes); this is stated in the caption text.
# ============================================================
def twin_plot(xdv, dv, xcv, cv, xl_dv, xl_cv, title, fname,
              band_dv=None, band_cv=None, band_dv_label=None, band_cv_label=None,
              logx_dv=False, logx_cv=False, output_dir=None):
    fig, axes = _panel_axes()
    all_series = []
    for ax, panel in zip(axes, PANELS):
        axt = ax.twiny()
        handles = []

        if band_dv is not None:
            s1 = ax.axvspan(band_dv[0], band_dv[1], alpha=0.15, color=BAND_DV,
                            zorder=0, label=band_dv_label or "DV deployed regime")
            ax.axvline(band_dv[0], color=DV_COLOR, ls=':', lw=1.0, alpha=0.7, zorder=1)
            ax.axvline(band_dv[1], color=DV_COLOR, ls=':', lw=1.0, alpha=0.7, zorder=1)
            handles.append(s1)
        if band_cv is not None:
            s2 = axt.axvspan(band_cv[0], band_cv[1], alpha=0.15, color=BAND_CV,
                             zorder=0, label=band_cv_label or "CV deployed regime")
            axt.axvline(band_cv[0], color=CV_COLOR, ls=':', lw=1.0, alpha=0.7, zorder=1)
            axt.axvline(band_cv[1], color=CV_COLOR, ls=':', lw=1.0, alpha=0.7, zorder=1)
            handles.append(s2)

        y_dv = _bps(dv, panel['dv'])
        y_cv = _bps(cv, panel['cv'])
        all_series += [y_dv, y_cv]
        h_dv, = ax.plot(xdv, y_dv, '-o', color=DV_COLOR, label='DV — decoy BB84')
        h_cv, = axt.plot(xcv, y_cv, '-s', color=CV_COLOR,
                         label='CV — GG02 heterodyne')
        handles += [h_dv, h_cv]

        if logx_dv:
            ax.set_xscale('log')
        if logx_cv:
            axt.set_xscale('log')
        ax.set_yscale('log')
        axt.set_yscale('log')

        ax.set_xlabel(xl_dv, color=DV_COLOR)
        axt.set_xlabel(xl_cv, color=CV_COLOR)
        ax.tick_params(axis='x', colors=DV_COLOR)
        axt.tick_params(axis='x', colors=CV_COLOR)
        ax.spines['bottom'].set_color(DV_COLOR)
        axt.spines['top'].set_color(CV_COLOR)

        ax.grid(True, which='both', alpha=0.25, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_title(panel['title'], pad=32, fontsize=11.5)  # pad: top axis label
        if ax is axes[0]:
            ax.legend(handles=handles, loc='best', framealpha=0.92, edgecolor='0.7')
    # Caveat: the two curves use DIFFERENT x-axes (DV bottom, CV top), so any
    # apparent intersection is an artefact of axis alignment, not a physical
    # equivalence. State this on the figure so it cannot be misread.
    # (twiny() shares the parent's y-axis, so clamping `axes` covers the twins.)
    _autoscale_y(axes, all_series)
    _finish_panels(fig, axes, title, note=TWIN_NOTE + " " + CLOCK_NOTE,
                   note_y=-0.155)
    _save_or_show(fig, output_dir, fname)


def shared_plot3(x, dv, cv_het, cv_hom, xlabel, title, fname,
                 bands=None, logx=False, output_dir=None, legend_loc='best'):
    """Shared x-axis with three curves: DV, CV heterodyne, CV homodyne, drawn as
    two clock panels (matched 1 GHz | deployed DV 1 GHz vs CV 100 MHz).
    Used where the parameter is identical for all (beta, alpha, detector eff)."""
    fig, axes = _panel_axes()
    all_series = []
    for ax, panel in zip(axes, PANELS):
        handles = []
        for (lo, hi), color, label in (bands or []):
            sp = ax.axvspan(lo, hi, alpha=0.16, color=color, zorder=0, label=label)
            ax.axvline(lo, color=color, ls=':', lw=1.0, alpha=0.8, zorder=1)
            ax.axvline(hi, color=color, ls=':', lw=1.0, alpha=0.8, zorder=1)
            handles.append(sp)
        y_dv  = _bps(dv,     panel['dv'])
        y_het = _bps(cv_het, panel['cv'])
        y_hom = _bps(cv_hom, panel['cv'])
        all_series += [y_dv, y_het, y_hom]
        h_dv,  = ax.plot(x, y_dv,  '-o', color=DV_COLOR, label='DV — decoy BB84')
        h_het, = ax.plot(x, y_het, '-s', color=CV_COLOR,
                         label='CV — GG02 heterodyne')
        h_hom, = ax.plot(x, y_hom, '--D', color='#7d2d8c', ms=3.5,
                         label='CV — GG02 homodyne')
        handles += [h_dv, h_het, h_hom]
        if logx:
            ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(xlabel)
        ax.grid(True, which='both', alpha=0.25, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_title(panel['title'], pad=8, fontsize=11.5)
        if ax is axes[0]:
            ax.legend(handles=handles, loc=legend_loc, framealpha=0.92,
                      edgecolor='0.7')
    _autoscale_y(axes, all_series)
    _finish_panels(fig, axes, title, note=CLOCK_NOTE)
    _save_or_show(fig, output_dir, fname)


def twin_plot3(xdv, dv, xcv, cv_het, cv_hom, xl_dv, xl_cv, title, fname,
               band_dv=None, band_cv=None, band_dv_label=None, band_cv_label=None,
               logx_dv=False, output_dir=None, legend_loc='best'):
    """Twin x-axes: DV (bottom) vs CV het+hom (top), drawn as two clock panels
    (matched 1 GHz | deployed DV 1 GHz vs CV 100 MHz). Crossings not meaningful."""
    fig, axes = _panel_axes()
    all_series = []
    for ax, panel in zip(axes, PANELS):
        axt = ax.twiny()
        handles = []
        if band_dv is not None:
            s1 = ax.axvspan(band_dv[0], band_dv[1], alpha=0.15, color=BAND_DV,
                            zorder=0, label=band_dv_label or "DV deployed regime")
            ax.axvline(band_dv[0], color=DV_COLOR, ls=':', lw=1.0, alpha=0.7, zorder=1)
            ax.axvline(band_dv[1], color=DV_COLOR, ls=':', lw=1.0, alpha=0.7, zorder=1)
            handles.append(s1)
        if band_cv is not None:
            s2 = axt.axvspan(band_cv[0], band_cv[1], alpha=0.15, color=BAND_CV,
                             zorder=0, label=band_cv_label or "CV deployed regime")
            axt.axvline(band_cv[0], color=CV_COLOR, ls=':', lw=1.0, alpha=0.7, zorder=1)
            axt.axvline(band_cv[1], color=CV_COLOR, ls=':', lw=1.0, alpha=0.7, zorder=1)
            handles.append(s2)
        y_dv  = _bps(dv,     panel['dv'])
        y_het = _bps(cv_het, panel['cv'])
        y_hom = _bps(cv_hom, panel['cv'])
        all_series += [y_dv, y_het, y_hom]
        h_dv,  = ax.plot(xdv, y_dv, '-o', color=DV_COLOR, label='DV — decoy BB84')
        h_het, = axt.plot(xcv, y_het, '-s', color=CV_COLOR, label='CV — heterodyne')
        h_hom, = axt.plot(xcv, y_hom, '--D', color='#7d2d8c', ms=3.5,
                          label='CV — homodyne')
        handles += [h_dv, h_het, h_hom]
        if logx_dv:
            ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(xl_dv, color=DV_COLOR); axt.set_xlabel(xl_cv, color=CV_COLOR)
        ax.tick_params(axis='x', colors=DV_COLOR)
        axt.tick_params(axis='x', colors=CV_COLOR)
        ax.spines['bottom'].set_color(DV_COLOR); axt.spines['top'].set_color(CV_COLOR)
        ax.grid(True, which='both', alpha=0.25, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_title(panel['title'], pad=32, fontsize=11.5)
        if ax is axes[0]:
            ax.legend(handles=handles, loc=legend_loc, framealpha=0.92,
                      edgecolor='0.7')
    _autoscale_y(axes, all_series)
    _finish_panels(fig, axes, title, note=TWIN_NOTE + " " + CLOCK_NOTE,
                   note_y=-0.155)
    _save_or_show(fig, output_dir, fname)


def arg_common():
    """Shared CLI: --output-dir PATH (save; omit to show) and --distance KM."""
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=str, default=None,
                   help="Directory to save the figure (omit to show interactively)")
    p.add_argument("--distance", type=float, default=L_KM,
                   help=f"Test distance in km (default {L_KM:.0f}, near the crossover)")
    return p.parse_args()


def arg_output_dir():
    """Back-compat shim for callers that only need the output directory."""
    return arg_common().output_dir


# ============================================================
# SWEEP RANGES + DEPLOYED BANDS
# Ranges: param.md deployed ranges, two modestly extended so the curve's full
# response is visible (beta to 0.80; dv_qber to the ~11% BB84 cutoff). Dark
# count kept realistic (1-1000 cps) — flat there is correct physics, corroborated
# by an independent NetSquid simulation, so the axis is not stretched.
#
# DEPLOYED/REALISTIC BANDS — each end sourced to a specific paper (June 2026).
# Bands reflect "best-demonstrated deployed/field" hardware for each protocol.
# Outliers excluded where noted (different detector tier / amplifier tricks).
#
#  DV eta   0.65-0.93  lo: Tang 2016 Hefei deployed metro SNSPD 0.64-0.66
#                      hi: Swedish 303km field trial 0.92-0.93 (arXiv:2606.06107)
#                      (NbN intercity 0.30, arXiv:1708.00434, excluded as low outlier)
#  CV eta   0.60-0.72  lo: Lodewyck 2007 (0.606); "most CV uses 0.6" (arXiv:1006.4216)
#                      hi: Si-photonic TBHD 0.72 (arXiv:2305.03419)
#                      (PSA-enhanced 0.965 and integrated 0.37 excluded: diff. tier)
#  beta     0.90-0.96  lo: Lodewyck/SECOQC field beta=0.9
#                      hi: rate-adaptive MET-LDPC 0.964 (arXiv:1703.04916)
#                      (older homodyne beta=0.8, arXiv:1006.4216, noted as floor)
#  alpha    0.18-0.21  lo: better-grade deployed SMF ~0.18 (ITU-T G.652.D)
#                      hi: deployed field fibre 0.19-0.21 (Tang 2016 et al.)
#  DV dark  1-1000 cps lo: Swedish 303km field SNSPD <=1 cps (arXiv:2606.06107)
#                      hi: Boston metro WSi/NbN ~1000 cps (arXiv:1708.00434)
#  CV v_el  0.015-0.11 lo: Jouguet 2013 deployed 80 km system, v_el=0.015
#                      hi: deployed field BPD 0.11 SNU (npj QI 2025, s41534-025-01060-7)
#  CV eps_b 5e-4-2e-3 lo: Wang 2019 prototype calibration; also reproduces the
#                          LuxQuanta NOVA LQ Gen-2 spec (100 km / 20 dB)
#                      hi: pessimistic deployed (reach ~69 km)
#                      NB total input-referred xi_r follows Eq. 12 and is
#                          distance-dependent; measured values run 0.005 SNU at
#                          0 km to ~0.098 SNU at 100 km (Wang Table 2)
#  DV QBER  0.5-2.1%   lo: Swedish field SNSPD link 0.5% (arXiv:2606.06107)
#                      hi: same trial InGaAs detector 2.1% (arXiv:2606.06107)
#
# NB asymmetry is real, not an artefact: DV eta band (0.65-0.93) sits above CV eta
# (0.60-0.72) because deployed SNSPDs genuinely outperform deployed homodyne. This
# supports the matched-maturity-tier fairness argument (see methodology.md).
# ============================================================
N = 40
RANGES = dict(
    dv_eta  = np.linspace(0.50, 1.00, N),    # axis 0.5-1.0; deployed band sits inside
    cv_eta  = np.linspace(0.50, 1.00, N),    # axis 0.5-1.0; deployed band sits inside
    beta    = np.linspace(0.80, 1.00, N),    # deployed 0.90-0.96 is a slice
    alpha   = np.linspace(0.01, 0.25, N),    # HCF-motivated: 0.01 (HCF projected) to 0.25
    dv_dark = np.logspace(0, 5, N),          # 1-1e5 cps; deployed 1-1000 is a slice
    cv_vel  = np.linspace(0.00, 0.50, N),    # deployed 0.05-0.11 is a slice
    dv_qber = np.linspace(0.001, 0.11, N),   # to ~11% cutoff; deployed 0.5-2.1% is a slice
    # Wang Eq. 12 decomposition. eps_b is the term that carries the distance
    # scaling (amplified by 1/(eta*T) at the channel input), so it is the
    # informative sensitivity axis; eps_a+eps_l is a constant floor.
    cv_xi_b  = np.linspace(1e-4, 5e-3, N),   # Bob-side (SNU); deployed 5e-4-2e-3
    cv_xi_al = np.linspace(0.0, 0.02, N),    # Alice+fibre floor (SNU)
)
BAND = dict(
    dv_eta=(0.65, 0.93), cv_eta=(0.60, 0.72), beta=(0.90, 0.96),
    alpha=(0.091, 0.20), dv_dark=(1, 1000), cv_vel=(0.015, 0.11),
    dv_qber=(0.005, 0.021),
    # eps_b lo: Wang 2019 calibration to their prototype (5e-4), which also
    # reproduces the 100 km / 20 dB reach of the LuxQuanta NOVA LQ Gen-2
    # commercial CV-QKD system. hi: 2e-3, a more pessimistic deployed figure
    # (reach ~69 km). Reach at the lo/hi ends: 94 km / 69 km.
    cv_xi_b=(5e-4, 2e-3),
    cv_xi_al=(0.003, 0.008),
)