"""
Point-to-point DV vs CV secret key rate comparison (thesis Figures 2-4).

DV: TNO asymptotic decoy BB84. CV: qosst-skr GG02 heterodyne and homodyne,
with Wang et al. 2019 Eq. 12 excess noise. Parameters and sourcing in param.md.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar  # bounded Va refine

from qosst_skr.gaussian_trusted_heterodyne_asymptotic import (
    GaussianTrustedHeterodyneAsymptotic)
from qosst_skr.gaussian_trusted_homodyne_asymptotic import (
    GaussianTrustedHomodyneAsymptotic)
from tno.quantum.communication.qkd_key_rate.quantum import standard_detector
from tno.quantum.communication.qkd_key_rate.quantum.bb84 import (
    BB84FullyAsymptoticKeyRateEstimate, compute_gain_and_error_rate)


# ============================================================
# EVALUATION PARAMETERS — deployed baseline (see param.md)
# ============================================================
ALPHA_DB_KM = 0.20         # SMF-28 at 1550 nm

# DV (TNO decoy BB84) — Clason et al. 2026, ID281 SNSPD.
# TNO's `efficiency_party` is Bob's TOTAL receiver efficiency, so a receiver
# optical budget is applied before Clason's detector SDE enters the engine.
DV_DETECTOR_ETA   = 0.92   # ID281 SNSPD SDE (Clason 2026)
DV_RX_LOSS_DB     = 2.0    # lumped receiver insertion loss (cf. Kelsey 2026)
DV_EFFICIENCY     = DV_DETECTOR_ETA * 10**(-DV_RX_LOSS_DB/10.0)   # = 0.580
DV_DARK_COUNT     = 3e-8   # per gate, per detector (Clason 2026)
DV_QBER           = 0.005  # intrinsic QBER (Clason 2026)
DV_F_EC           = 1.036  # Cascade f_EC (Mueller 2024); applied post-hoc
                           # in dv_skr since TNO is at the Shannon limit

# CV (qosst-skr GG02) — deployed homodyne
CV_ETA            = 0.60   # homodyne efficiency (Lodewyck 2007, Zhang 2020)
CV_VEL            = 0.10   # electronic noise, SNU
CV_BETA           = 0.95   # Zhang et al. 2020 Table I

# CV excess noise, channel-input referred (Wang et al. 2019, Eq. 12):
#   xi_r(T) = CV_XI_AL + CV_XI_B / (CV_ETA * T)
CV_XI_AL          = 0.005   # Alice + fibre (SNU)
CV_XI_B           = 0.0005  # Bob-side measurement (SNU)

# Idealised CV reference (Wang Fig. 8 fixed-noise line), for Figure 2.
CV_XI_IDEAL       = 0.01    # SNU, constant
INCLUDE_IDEALISED = True

# Clock rates for Figure 4
DV_CLOCK_HZ       = 1e9
CV_CLOCK_HZ       = 100e6

DISTANCES_KM = np.arange(0, 321, 0.5)
RATE_FLOOR   = 1e-10        # bits/channel-use plotting floor
VA_LO, VA_HI = 1e-2, 10.0   # Va search bounds (SNU)

# ============================================================
# HELPERS
# ============================================================
def T_of_L(L_km):
    return 10**(-ALPHA_DB_KM * L_km / 10.0)


def cv_xi_input(T):
    """Channel-input-referred CV excess noise, Wang et al. 2019 Eq. 12."""
    return CV_XI_AL + CV_XI_B / (CV_ETA * T)


def _h(p):
    """Binary entropy in bits."""
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -p*np.log2(p) - (1-p)*np.log2(1-p)


# ============================================================
# Va OPTIMISATION (shared by heterodyne and homodyne)
# ============================================================
def _optimise_va(skr_fn, T, xi):
    """
    Maximise secret key rate over modulation variance Va.

    Coarse log-spaced scan to locate the basin, then a bounded scalar refine
    inside it. Log spacing because the GG02 optimum sits at small Va and the
    rate-vs-Va curve is sharply peaked there; a uniform linear grid straddles
    the peak and under-reports the rate. Returns (best_rate, best_va, hit_edge).
    """
    def neg_rate(va):
        try:
            r = skr_fn(Va=va, T=T, xi=xi, eta=CV_ETA, Vel=CV_VEL, beta=CV_BETA)
            return -r if (r is not None and np.isfinite(r)) else 0.0
        except Exception:
            return 0.0

    grid = np.logspace(np.log10(VA_LO), np.log10(VA_HI), 120)
    vals = np.array([neg_rate(va) for va in grid])
    i = int(np.argmin(vals))
    best_va, best_neg = grid[i], vals[i]

    lo = grid[max(i - 1, 0)]
    hi = grid[min(i + 1, len(grid) - 1)]
    if hi > lo:
        try:
            res = minimize_scalar(neg_rate, bounds=(lo, hi), method="bounded",
                                  options={"xatol": 1e-4})
            if res.fun < best_neg:
                best_va, best_neg = res.x, res.fun
        except Exception:
            pass

    hit_edge = (i == 0) or (i == len(grid) - 1)
    return max(-best_neg, 0.0), best_va, hit_edge


# ============================================================
# CV ENGINES (qosst-skr), each optimised over Va
# ============================================================
def cv_het_skr(L_km):
    """CV heterodyne, bits/symbol, qosst-skr, Wang Eq. 12 excess noise."""
    T = T_of_L(L_km)
    if T < 1e-12:
        return 0.0
    rate, va, hit_edge = _optimise_va(
        GaussianTrustedHeterodyneAsymptotic.skr, T, cv_xi_input(T))
    # Only warn on lower-bound hits (a real problem); upper bound is intentional.
    if hit_edge and rate > 0 and va < VA_LO * 1.5:
        print(f"  [warn] het Va optimum at LOWER bound (Va={va:.3g}) @ {L_km} km")
    return rate


def cv_hom_skr(L_km):
    """CV homodyne, bits/symbol, qosst-skr, Wang Eq. 12 excess noise."""
    T = T_of_L(L_km)
    if T < 1e-12:
        return 0.0
    rate, va, hit_edge = _optimise_va(
        GaussianTrustedHomodyneAsymptotic.skr, T, cv_xi_input(T))
    if hit_edge and rate > 0 and va < VA_LO * 1.5:
        print(f"  [warn] hom Va optimum at LOWER bound (Va={va:.3g}) @ {L_km} km")
    return rate


def cv_hom_ideal_skr(L_km):
    """Idealised CV homodyne reference: constant xi (Wang Fig. 8 fixed line)."""
    T = T_of_L(L_km)
    if T < 1e-12:
        return 0.0
    rate, _, _ = _optimise_va(
        GaussianTrustedHomodyneAsymptotic.skr, T, CV_XI_IDEAL)
    return rate


# ============================================================
# DV ENGINE (TNO) + post-hoc reconciliation DV_BETA
# ============================================================
def dv_skr(L_km):
    """DV decoy BB84 bits/pulse, asymptotic, with error correction at DV_F_EC.

    TNO's engine hardcodes error correction at the Shannon limit (f_EC = 1);
    a realistic f_EC > 1 is applied post-hoc by scaling the EC cost:
        rate = rate_ideal - (f_EC - 1) * gain * h(QBER).
    This is the standard f-notation cost in the decoy-BB84 rate.

    NB: TNO's asymptotic BB84 derives the channel error rate from dark counts and
    POLARIZATION DRIFT, not from `error_detector` (unused in this rate path). The
    intrinsic QBER is encoded as a drift angle: polarization_drift =
    arcsin(sqrt(QBER)) reproduces the target QBER to <1e-3 (verified against
    compute_gain_and_error_rate). Without this, DV runs at ~zero misalignment
    error and is optimistic.
    """
    drift = float(np.arcsin(np.sqrt(np.clip(DV_QBER, 0.0, 0.5))))  # QBER -> drift
    det = standard_detector.customise(
        efficiency_party=DV_EFFICIENCY,
        dark_count_rate=DV_DARK_COUNT,
        error_detector=DV_QBER,
        polarization_drift=drift
    )
    att = ALPHA_DB_KM * L_km
    try:
        mu_opt, r_ideal = BB84FullyAsymptoticKeyRateEstimate(
            detector=det).optimize_rate(attenuation=att)
    except Exception:
        return 0.0
    r_ideal = max(r_ideal, 0.0)
    if DV_F_EC <= 1.0 or r_ideal <= 0.0:
        return r_ideal
    try:
        mu = float(np.atleast_1d(mu_opt["mu"])[0])   # optimize_rate returns {"mu": array}
        gain, err = compute_gain_and_error_rate(det, mu, att)
        ec_extra = (DV_F_EC - 1.0) * float(np.atleast_1d(gain)[0]) * _h(float(np.atleast_1d(err)[0]))
        return max(r_ideal - ec_extra, 0.0)
    except Exception:
        return r_ideal


# ============================================================
# COMPUTE SWEEPS
# ============================================================
print("Computing DV sweep...")
dv_rates = np.array([dv_skr(L) for L in DISTANCES_KM])
print("Computing CV heterodyne sweep (Wang Eq. 12 noise)...")
cv_het_rates = np.array([cv_het_skr(L) for L in DISTANCES_KM])
print("Computing CV homodyne sweep (Wang Eq. 12 noise)...")
cv_hom_rates = np.array([cv_hom_skr(L) for L in DISTANCES_KM])
if INCLUDE_IDEALISED:
    print(f"Computing idealised CV reference (const xi = {CV_XI_IDEAL})...")
    cv_ideal_rates = np.array([cv_hom_ideal_skr(L) for L in DISTANCES_KM])
else:
    cv_ideal_rates = np.zeros_like(dv_rates)


# ---- Validation: excess-noise model vs Wang Table 2 (measured prototype) ----
WANG_TABLE2 = {0: 0.00517, 20: 0.01337, 50: 0.01582, 80: 0.07212, 100: 0.09791}
print("\nExcess-noise model vs Wang et al. Table 2 (measured, input-referred):")
print(f"{'L(km)':>6} {'model xi_r':>11} {'measured':>9}")
for L, meas in WANG_TABLE2.items():
    print(f"{L:6d} {cv_xi_input(T_of_L(L)):11.5f} {meas:9.5f}")

# ---- Crossover: where DV overtakes the better CV scheme (bits/channel use) ----
cv_best = np.maximum(cv_het_rates, cv_hom_rates)
both_live = (dv_rates > 1e-7) & (cv_best > 1e-7)
diff = np.sign(dv_rates - cv_best)
crossover_km = None
for i in range(1, len(DISTANCES_KM)):
    if both_live[i] and both_live[i - 1] and diff[i] != diff[i - 1] and diff[i] != 0:
        d0, d1 = (dv_rates - cv_best)[i - 1], (dv_rates - cv_best)[i]
        x0, x1 = DISTANCES_KM[i - 1], DISTANCES_KM[i]
        crossover_km = x0 + (x1 - x0) * (0 - d0) / (d1 - d0)
        break
if crossover_km is not None:
    print(f"\nDV/CV crossover (bits per channel use): {crossover_km:.1f} km")
else:
    side = '>' if dv_rates[1] > cv_best[1] else '<'
    print(f"\nNo DV/CV crossover found where both protocols are live "
          f"(DV {side} CV throughout)")

# ---- Reach of each protocol under this parameter set ----
reach_rows = [("DV", dv_rates), ("CV-het", cv_het_rates), ("CV-hom", cv_hom_rates)]
if INCLUDE_IDEALISED:
    reach_rows.append((f"CV idealised (xi={CV_XI_IDEAL})", cv_ideal_rates))
for name, rates in reach_rows:
    live = DISTANCES_KM[rates > RATE_FLOOR]
    reach = live.max() if live.size else 0
    note = "  (runs past axis limit)" if live.size and reach == DISTANCES_KM.max() else ""
    print(f"{name} maximum distance with key: {reach} km{note}")


# ============================================================
# FIGURE 2: bits per channel use (protocol comparison)
# ============================================================
fig2, ax2 = plt.subplots(figsize=(9, 6))

mask_dv = dv_rates > RATE_FLOOR
ax2.semilogy(DISTANCES_KM[mask_dv], dv_rates[mask_dv], 'b-', linewidth=2.5,
             label='DV: Decoy BB84 (TNO)')

mask_het = cv_het_rates > RATE_FLOOR
ax2.semilogy(DISTANCES_KM[mask_het], cv_het_rates[mask_het], 'r-', linewidth=2.5,
             label='CV: GG02 Heterodyne (qosst-skr)')

mask_hom = cv_hom_rates > RATE_FLOOR
ax2.semilogy(DISTANCES_KM[mask_hom], cv_hom_rates[mask_hom], 'm--', linewidth=2.5,
             label='CV: GG02 Homodyne (qosst-skr)')

if INCLUDE_IDEALISED:
    mask_id = cv_ideal_rates > RATE_FLOOR
    ax2.semilogy(DISTANCES_KM[mask_id], cv_ideal_rates[mask_id], color='0.45',
                 linewidth=1.5, linestyle=':', zorder=1,
                 label=rf'CV idealised: const $\xi$={CV_XI_IDEAL} (hom)')

ax2.set_xlabel('Distance (km)', fontsize=12)
ax2.set_ylabel('Secret key rate (bits per channel use)', fontsize=12)
ax2.set_title('DV vs CV — Protocol Comparison (asymptotic)', fontsize=12)
ax2.set_xlim(0, 320)
ax2.set_ylim(RATE_FLOOR, 10)
ax2.legend(fontsize=10, loc='lower left')
ax2.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('figure2_bits_per_channel_use.png', dpi=150, bbox_inches='tight')
print("Saved: figure2_bits_per_channel_use.png")

# ============================================================
# FIGURE 3: bits per second at EQUAL clock rates (1 GHz)
# ============================================================
fig3, ax3 = plt.subplots(figsize=(9, 6))

EQUAL_CLOCK_HZ = 1e9

dv_equal_bps     = dv_rates * EQUAL_CLOCK_HZ
cv_het_equal_bps = cv_het_rates * EQUAL_CLOCK_HZ
cv_hom_equal_bps = cv_hom_rates * EQUAL_CLOCK_HZ

mask_dv3 = dv_equal_bps > 1
ax3.semilogy(
    DISTANCES_KM[mask_dv3],
    dv_equal_bps[mask_dv3],
    'b-',
    linewidth=2.5,
    label='DV: Decoy BB84 @ 1 GHz'
)

mask_het3 = cv_het_equal_bps > 1
ax3.semilogy(
    DISTANCES_KM[mask_het3],
    cv_het_equal_bps[mask_het3],
    'r-',
    linewidth=2.5,
    label='CV: GG02 Heterodyne @ 1 GHz'
)

mask_hom3 = cv_hom_equal_bps > 1
ax3.semilogy(
    DISTANCES_KM[mask_hom3],
    cv_hom_equal_bps[mask_hom3],
    'm--',
    linewidth=2.5,
    label='CV: GG02 Homodyne @ 1 GHz'
)

ax3.set_xlabel('Distance (km)', fontsize=12)
ax3.set_ylabel('Secret key rate (bits/s)', fontsize=12)
ax3.set_title('DV vs CV — Equal Clock Rate (1 GHz)', fontsize=12)

ax3.set_xlim(0, 320)
ax3.set_ylim(1e0, 1e10)

ax3.legend(fontsize=10, loc='lower left')
ax3.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('figure3_equal_clock_1GHz.png', dpi=150, bbox_inches='tight')
print("Saved: figure3_equal_clock_1GHz.png")


# ============================================================
# FIGURE 4: bits per second at realistic clocks (deployment view)
# ============================================================
fig4, ax4 = plt.subplots(figsize=(9, 6))

dv_bps     = dv_rates     * DV_CLOCK_HZ
cv_het_bps = cv_het_rates * CV_CLOCK_HZ
cv_hom_bps = cv_hom_rates * CV_CLOCK_HZ

mask_dv4 = dv_bps > 1
ax4.semilogy(DISTANCES_KM[mask_dv4], dv_bps[mask_dv4], 'b-', linewidth=2.5,
             label=f'DV: Decoy BB84 @ {DV_CLOCK_HZ/1e9:.1f} GHz')

mask_het4 = cv_het_bps > 1
ax4.semilogy(DISTANCES_KM[mask_het4], cv_het_bps[mask_het4], 'r-', linewidth=2.5,
             label=f'CV: Heterodyne @ {CV_CLOCK_HZ/1e6:.0f} MHz')

mask_hom4 = cv_hom_bps > 1
ax4.semilogy(DISTANCES_KM[mask_hom4], cv_hom_bps[mask_hom4], 'm--', linewidth=2.5,
             label=f'CV: Homodyne @ {CV_CLOCK_HZ/1e6:.0f} MHz')

ax4.set_xlabel('Distance (km)', fontsize=12)
ax4.set_ylabel('Secret key rate (bits/s)', fontsize=12)
ax4.set_title('DV vs CV — Deployment View (realistic per-protocol clocks)', fontsize=12)
ax4.set_xlim(0, 320)
ax4.set_ylim(1e0, 1e9)
ax4.legend(fontsize=10, loc='lower left')
ax4.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('figure4_bits_per_second.png', dpi=150, bbox_inches='tight')
print("Saved: figure4_bits_per_second.png")


# ============================================================
# PRINT TABULAR RESULTS
# ============================================================
print(f"\n{'L(km)':>6} {'xi_r':>9} {'DV b/pulse':>11} {'CV-het b/sym':>13} "
      f"{'CV-hom b/sym':>13} {'DV bps':>11} {'CV-het bps':>12} {'CV-hom bps':>12}")
for i, L in enumerate(DISTANCES_KM):
    if L % 5 == 0:
        print(f"{L:6.0f} {cv_xi_input(T_of_L(L)):>9.4g} {dv_rates[i]:11.3e} "
              f"{cv_het_rates[i]:13.3e} {cv_hom_rates[i]:13.3e} "
              f"{dv_rates[i]*DV_CLOCK_HZ:11.3e} {cv_het_rates[i]*CV_CLOCK_HZ:12.3e} "
              f"{cv_hom_rates[i]*CV_CLOCK_HZ:12.3e}")