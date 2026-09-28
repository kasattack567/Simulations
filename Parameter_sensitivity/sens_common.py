"""
Shared engine wrappers, locked baseline, and plot helpers for the DV vs CV
parameter sensitivity analysis (fixed distance, near the crossover).

Imported by s1..s5; not run directly. Centralising the baseline and engine
calls keeps them identical to cv_hetro.py / param.md.

qosst-skr and TNO are deterministic analytic formulas, so no error bars.

CV excess noise follows Wang et al. 2019 Eq. 12 (channel-input referred):
    xi_r(T) = (eps_a + eps_l) + eps_b / (eta * T)
DV f_EC > 1 is applied post-hoc, as in cv_hetro.py (TNO is at Shannon f_EC=1).

Figures plot bit/s as two panels sharing a y-axis: (a) matched clock at 1 GHz
for both, (b) deployed clocks DV 1 GHz / CV 100 MHz. cv_rate/dv_rate still
return bits per channel use; only the plot helpers apply the clock scaling.
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

L_KM = 59.015
ALPHA_DB_KM = 0.20

DV_DETECTOR_ETA = 0.92    
DV_RX_LOSS_DB   = 2.0      
DV_EFFICIENCY = DV_DETECTOR_ETA * 10**(-DV_RX_LOSS_DB/10.0)   
DV_DARK_COUNT = 3e-8      
DV_QBER       = 0.005      
DV_F_EC       = 1.036      

CV_ETA    = 0.60
CV_VEL    = 0.10
CV_BETA   = 0.95           
CV_XI_AL  = 0.005     
CV_XI_B   = 0.0005    
VA_LO, VA_HI = 1e-2, 10.0
DARK_CPS_TO_PERGATE = 1e-9  

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


def dv_rate(L_km=L_KM, eta=None, dark_pergate=None, qber=None, f_ec=None, alpha=None):
    """DV decoy BB84 bits/pulse, asymptotic, with post-hoc error correction f_EC.

    f_ec is the DV-convention error-correction efficiency (leak = f_EC * Q * h(E)),
    NOT a CV-convention reconciliation efficiency beta. A caller sweeping a shared
    beta axis against CV must pass f_ec = 1/beta and label the axis accordingly.

    NB: TNO's asymptotic BB84 derives the channel error rate from dark counts and
    POLARIZATION DRIFT, not from the detector's `error_detector` field (which is
    unused in this rate path). Intrinsic/misalignment QBER is therefore encoded as
    a drift angle: polarization_drift = arcsin(sqrt(QBER)) reproduces the target
    QBER to <1e-3 across 0-11% (verified against compute_gain_and_error_rate).
    """
    eta = DV_EFFICIENCY if eta is None else eta
    dark_pergate = DV_DARK_COUNT if dark_pergate is None else dark_pergate
    qber = DV_QBER if qber is None else qber
    f_ec = DV_F_EC if f_ec is None else f_ec
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
    if f_ec <= 1.0 or r_ideal <= 0.0:
        return r_ideal
    try:
        mu = float(np.atleast_1d(mu_opt["mu"])[0])   # optimize_rate returns {"mu": array}
        gain, err = compute_gain_and_error_rate(det, mu, att)
        ec_extra = (f_ec - 1.0) * float(np.atleast_1d(gain)[0]) * _h(float(np.atleast_1d(err)[0]))
        return max(r_ideal - ec_extra, 0.0)
    except Exception:
        return r_ideal


# ============================================================
# PLOT STYLE
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
# SWEEP RANGES + DEPLOYED BANDS (param.md; June 2026 refs)
# Ranges cover deployed hardware; beta and DV QBER extended to the curve's
# full response. Bands are "best-demonstrated deployed/field" per protocol;
# outliers (different detector tier, amplifier tricks) excluded.
#
#  DV eta   0.65-0.93   Tang 2016 Hefei / Swedish 303km field (arXiv:2606.06107)
#  CV eta   0.60-0.72   Lodewyck 2007 / Si-photonic TBHD (arXiv:2305.03419)
#  beta     0.90-0.96   Lodewyck-SECOQC / MET-LDPC (arXiv:1703.04916)
#  alpha    0.18-0.21   ITU-T G.652.D / deployed field fibre
#  DV dark  1-1000 cps  Swedish 303km / Boston metro (arXiv:1708.00434)
#  CV v_el  0.015-0.11  Jouguet 2013 / deployed BPD (npj QI 2025)
#  CV eps_b 5e-4-2e-3   Wang 2019 (also matches LuxQuanta NOVA LQ Gen-2 spec)
#  DV QBER  0.5-2.1 %   Swedish trial SNSPD / InGaAs (arXiv:2606.06107)
#
# xi_r follows Wang Eq. 12 and is distance-dependent (~0.005 SNU at 0 km to
# ~0.098 at 100 km, Wang Table 2). DV eta band sits above CV eta because
# ============================================================

N = 40
RANGES = dict(
    dv_eta  = np.linspace(0.35, 1.00, N),    # Bob TOTAL efficiency (eta_B), not detector SDE
    cv_eta  = np.linspace(0.50, 1.00, N),    # axis 0.5-1.0; deployed band sits inside
    beta    = np.linspace(0.80, 1.00, N),    # CV reconciliation; deployed 0.90-0.96 is a slice
    dv_f_ec = np.linspace(1.00, 1.40, N),    # DV error correction; deployed 1.036-1.20 is a slice
    alpha   = np.linspace(0.01, 0.25, N),    # HCF-motivated: 0.01 (HCF projected) to 0.25
    dv_dark = np.logspace(0, 5, N),          # 1-1e5 cps; Clason band 1-70 is a slice
    cv_vel  = np.linspace(0.00, 0.50, N),    # deployed 0.05-0.11 is a slice
    dv_qber = np.linspace(0.001, 0.11, N),   # to ~11% cutoff; deployed 0.5-2.1% is a slice
    # Wang Eq. 12 decomposition. eps_b is the term that carries the distance
    # scaling (amplified by 1/(eta*T) at the channel input), so it is the
    # informative sensitivity axis; eps_a+eps_l is a constant floor.
    cv_xi_b  = np.linspace(1e-4, 5e-3, N),   # Bob-side (SNU); deployed 5e-4-2e-3
    cv_xi_al = np.linspace(0.0, 0.02, N),    # Alice+fibre floor (SNU)
)
BAND = dict(
    # --- DV ---
    # Bob TOTAL eff = SDE 0.92 (Clason 2026) x receiver optics 1.0-3.5 dB
    # (Kelsey 2026 / Chen 2009). Baseline 2.0 dB -> 0.580.
    dv_eta=(0.411, 0.731),
    # Clason 2026 (<=1 cps X-basis, ~30 key, ~10 ambient) to UK-Ireland 70 cps.
    dv_dark=(1, 70),
    # Clason 2026: 0.5% (110 km SNSPD) to 1.8% (143 km MCF).
    dv_qber=(0.005, 0.018),

    # --- CV ---
    # Measured receiver eff: Pi 2023, Fossier 2009, Zhang 2019/2020, Hajomer 2024.
    cv_eta=(0.56, 0.68),
    # Fossier 2009 to Zhang 2020 in-system.
    cv_vel=(0.01, 0.27),
    # Inverted from Wang 2019 Table 2 at 0/20/50/80/100 km (eps_al = 0.005).
    cv_xi_b=(1e-4, 2e-3),
    # Wang 0 km point, +/-1 SD.
    cv_xi_al=(0.002, 0.007),

    # --- reconciliation (two distinct parameters) ---
    # CV beta, multiplies I_AB. Hajomer 2024 to Wang 2017 MET-LDPC. Zhang's 0.98
    # (202 km lab) and Milicevic 0.99 (sim) excluded. Engine ignores FER, so
    # high beta w/o its FER isn't face-value usable — argues for the low end.
    beta=(0.90, 0.96),
    # DV f_EC, scales leakage. Related to beta by Mueller 2024 Eq. 6, NOT 1/beta.
    # Mueller 2024 Cascade 1.036 / LDPC 1.166; Xu et al. RMP 2020 1.1-1.2.
    # Band shifts crossover ~1 km, reach 303.4 -> 298.2 km — nearly inert vs beta.
    dv_f_ec=(1.036, 1.20),

    # --- shared ---
    # Not a deployed spread (silica ~0.2 dB/km, flat for 40 yrs). Motivated by
    # hollow-core: 0.091 dB/km demonstrated (Petrovich et al., Nat. Photon. 19,
    # 1203, 2025), ~0.01 projected. Forward-looking; CV is loss-limited.
    alpha=(0.091, 0.20),
)


# ============================================================
# BAND LEGEND LABELS — generated from BAND, never hardcoded
# ============================================================
# Every band label used to be a hand-typed string duplicating the numbers in
# BAND. They drifted apart the first time a band moved. band_label() formats the
# label straight from BAND so the legend and the shaded region cannot disagree.
_BAND_FMT = {
    "dv_eta":   ("DV Bob total efficiency",            "{:.2f}",  ""),
    "cv_eta":   ("CV receiver efficiency",             "{:.2f}",  ""),
    "beta":     (r"CV deployed $\beta$",               "{:.2f}",  ""),
    "dv_f_ec":  (r"DV deployed $f_{EC}$",              "{:.3f}",  ""),
    "alpha":    ("HCF demonstrated to silica",         "{:.3g}",  " dB/km"),
    "dv_dark":  ("DV deployed",                        "{:.0f}",  " cps"),
    "cv_vel":   ("CV deployed",                        "{:.3g}",  " SNU"),
    "dv_qber":  ("DV deployed",                        "{:.1f}",  "%"),
    "cv_xi_b":  (r"Wang Table 2 implied $\epsilon_b$", "{:.1f}",  r" $\times10^{-3}$ SNU"),
    "cv_xi_al": (r"Wang 0 km $\pm1\sigma$",            "{:.3g}",  " SNU"),
}
_BAND_SCALE = {"dv_qber": 100.0, "cv_xi_b": 1e3}


def band_label(key, prefix=None):
    """Legend label for BAND[key], with the numbers taken from BAND itself."""
    lo, hi = BAND[key]
    name, fmt, unit = _BAND_FMT[key]
    k = _BAND_SCALE.get(key, 1.0)
    return f"{prefix or name} ({fmt.format(lo*k)}-{fmt.format(hi*k)}{unit})"


def band_scaled(key):
    """BAND[key] in the units the corresponding figure plots in."""
    k = _BAND_SCALE.get(key, 1.0)
    return (BAND[key][0]*k, BAND[key][1]*k)