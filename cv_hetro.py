"""
Final point-to-point DV vs CV comparison — thesis figures.

CV computed via qosst-skr (consistent framework for both detection schemes):
  - GaussianTrustedHeterodyneAsymptotic
  - GaussianTrustedHomodyneAsymptotic
DV via TNO asymptotic decoy BB84.

EVALUATION PARAMETERS: "best demonstrated in metropolitan deployments /
current commercial-grade hardware" (see param.md for full sourcing).
DV uses cooled-SNSPD metropolitan numbers (Tang 2016 deployed Hefei network);
CV uses deployed homodyne numbers (Lodewyck 2007 / Zhang 2020). This puts the
two protocols on comparable detector footing (eta 0.65 vs 0.60).

DV reconciliation: TNO's BB84FullyAsymptoticKeyRateEstimate hardcodes the
error-correction term at the Shannon limit (beta=1, ideal) with no argument to
change it. DV_BETA = 0.95 is therefore applied POST-HOC in dv_skr (matching CV's
tier and the locked param.md set), by scaling the EC cost term. See dv_skr.

Figure 2: bits per channel use (protocol comparison), DV + CV-het + CV-hom + bounds
Figure 4: bits/s at realistic clocks (deployment view)
Figure 3: DV extended range sanity check (0-400 km)
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
# EVALUATION PARAMETERS — metropolitan deployed baseline (see param.md)
# ============================================================
ALPHA_DB_KM = 0.20         # SMF-28 at 1550 nm; deployed fibre 0.19-0.21 dB/km

# DV (TNO decoy BB84) — deployed metro SNSPD receiver (Tang 2016 Hefei network)
DV_EFFICIENCY     = 0.65   # SNSPD system efficiency (Tang 2016: 0.64-0.66)
DV_DARK_COUNT     = 1e-7   # per gate: 100 cps at ~1 ns window / 1 GHz clock
DV_QBER           = 0.005  # well-aligned source; field metro QBER <1%
DV_BETA           = 0.95   # reconciliation, matched to CV tier. NB: TNO engine is
                           # ideal (beta=1); 0.95 applied post-hoc in dv_skr.

# CV (qosst-skr GG02) — deployed homodyne (Lodewyck 2007 / Zhang 2020)
CV_ETA            = 0.60   # homodyne efficiency (Lodewyck 0.606, Zhang 0.613)
CV_VEL            = 0.10   # electronic noise, deployment-grade homodyne (SNU)
CV_XI_BOB         = 0.01  # Bob-side excess noise (SNU); factor-of-2 corrected below
CV_BETA           = 0.95   # reconciliation efficiency, deployed CV standard

# Clock rates for Figure 4 (realistic deployment)
DV_CLOCK_HZ       = 1e9     # 1 GHz, typical modern DV system
CV_CLOCK_HZ       = 100e6   # 100 MHz, representative deployed CV symbol rate

# Distance range
DISTANCES_KM = np.arange(0, 81, 1)

# Va search bounds (modulation variance, SNU). Kept unconstrained: the GG02
# asymptotic optimum exceeds ~10 SNU only at sub-km distance (not physically
# achievable on real modulators, but outside the metro range of interest). The
# crossover and all points >=5 km are unaffected; 0-1 km CV is upper-idealised.
VA_LO, VA_HI = 1e-2, 100.0


# ============================================================
# HELPERS
# ============================================================
def T_of_L(L_km):
    return 10**(-ALPHA_DB_KM * L_km / 10.0)


def cv_xi_input(T):
    """Convert Bob-referred excess noise to channel-input excess noise for qosst-skr.

    Convention taken directly from QOSST's own documentation example, which
    initialises the channel as:
        GaussianChannel(transmittance, xi_bob / (transmittance * detector.eta))
    i.e. xi_input = xi_bob / (T * eta)   — NO factor of 2.

    NOTE: an earlier version used 2*xi_bob/(T*eta); the factor of 2 is NOT in
    QOSST's reference example and made the excess noise 2x too high, shortening
    the CV range. This form matches the QOSST example and is calibrated against
    it (see check_qosst_example.py).
    """
    return CV_XI_BOB / (T * CV_ETA)


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
    """CV heterodyne, bits/symbol, qosst-skr."""
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
    """CV homodyne, bits/symbol, qosst-skr."""
    T = T_of_L(L_km)
    if T < 1e-12:
        return 0.0
    rate, va, hit_edge = _optimise_va(
        GaussianTrustedHomodyneAsymptotic.skr, T, cv_xi_input(T))
    if hit_edge and rate > 0 and va < VA_LO * 1.5:
        print(f"  [warn] hom Va optimum at LOWER bound (Va={va:.3g}) @ {L_km} km")
    return rate


# ============================================================
# DV ENGINE (TNO) + post-hoc reconciliation DV_BETA
# ============================================================
def dv_skr(L_km):
    """DV decoy BB84 bits/pulse, asymptotic, with reconciliation DV_BETA.

    TNO's engine is ideal-reconciliation (beta=1); DV_BETA<1 is applied post-hoc
    by scaling the error-correction cost:
        rate = rate_ideal - (1/beta - 1) * gain * h(QBER).
    This is the standard way reconciliation efficiency enters the BB84 rate.

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
    if DV_BETA >= 1.0 or r_ideal <= 0.0:
        return r_ideal
    try:
        mu = float(np.atleast_1d(mu_opt["mu"])[0])   # optimize_rate returns {"mu": array}
        gain, err = compute_gain_and_error_rate(det, mu, att)
        ec_extra = (1.0/DV_BETA - 1.0) * float(np.atleast_1d(gain)[0]) * _h(float(np.atleast_1d(err)[0]))
        return max(r_ideal - ec_extra, 0.0)
    except Exception:
        return r_ideal


# ============================================================
# BOUNDS
# ============================================================
def plob_lb(L_km):
    eta = T_of_L(L_km)
    if eta >= 0.999999:
        return float('inf')
    return -np.log2(1 - eta)


def tgw_ub(L_km):
    eta = T_of_L(L_km)
    if eta >= 0.999999:
        return float('inf')
    return np.log2((1 + eta) / (1 - eta))


# ============================================================
# COMPUTE SWEEPS
# ============================================================
print("Computing DV sweep...")
dv_rates = np.array([dv_skr(L) for L in DISTANCES_KM])
print("Computing CV heterodyne sweep...")
cv_het_rates = np.array([cv_het_skr(L) for L in DISTANCES_KM])
print("Computing CV homodyne sweep...")
cv_hom_rates = np.array([cv_hom_skr(L) for L in DISTANCES_KM])
print("Computing bounds...")
plob_rates = np.array([plob_lb(L) for L in DISTANCES_KM])
tgw_rates = np.array([tgw_ub(L) for L in DISTANCES_KM])


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
    print(f"\nNo DV/CV crossover in 0-50 km (DV {side} CV throughout)")


# ============================================================
# FIGURE 2: bits per channel use (protocol comparison)
# ============================================================
fig2, ax2 = plt.subplots(figsize=(9, 6))

ax2.semilogy(DISTANCES_KM[1:], tgw_rates[1:], 'k-', linewidth=1, alpha=0.5,
             label='TGW upper bound')
ax2.semilogy(DISTANCES_KM[1:], plob_rates[1:], 'k--', linewidth=1, alpha=0.5,
             label='PLOB lower bound')

mask_dv = dv_rates > 1e-7
ax2.semilogy(DISTANCES_KM[mask_dv], dv_rates[mask_dv], 'b-', linewidth=2.5,
             label='DV: Decoy BB84 (TNO)')

mask_het = cv_het_rates > 1e-7
ax2.semilogy(DISTANCES_KM[mask_het], cv_het_rates[mask_het], 'r-', linewidth=2.5,
             label='CV: GG02 Heterodyne (qosst-skr)')

mask_hom = cv_hom_rates > 1e-7
ax2.semilogy(DISTANCES_KM[mask_hom], cv_hom_rates[mask_hom], 'm--', linewidth=2.5,
             label='CV: GG02 Homodyne (qosst-skr)')

ax2.set_xlabel('Distance (km)', fontsize=12)
ax2.set_ylabel('Secret key rate (bits per channel use)', fontsize=12)
ax2.set_title('DV vs CV — Protocol Comparison (asymptotic, metro range)', fontsize=12)
ax2.set_xlim(0, 80)
ax2.set_ylim(1e-5, 10)
ax2.legend(fontsize=10, loc='lower left')
ax2.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('figure2_bits_per_channel_use.png', dpi=150, bbox_inches='tight')
print("Saved: figure2_bits_per_channel_use.png")


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
ax4.set_xlim(0, 80)
ax4.set_ylim(1e3, 1e9)
ax4.legend(fontsize=10, loc='lower left')
ax4.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('figure4_bits_per_second.png', dpi=150, bbox_inches='tight')
print("Saved: figure4_bits_per_second.png")


# ============================================================
# FIGURE 3: DV extended range (sanity check — does DV eventually die?)
# ============================================================
DISTANCES_KM_LONG = np.arange(0, 301, 5)   # 0-300 km, 5 km steps

print("Computing DV extended sweep (0-400 km)...")
dv_rates_long = np.array([dv_skr(L) for L in DISTANCES_KM_LONG])
plob_long = np.array([plob_lb(L) for L in DISTANCES_KM_LONG])
tgw_long  = np.array([tgw_ub(L) for L in DISTANCES_KM_LONG])

fig3, ax3 = plt.subplots(figsize=(9, 6))

ax3.semilogy(DISTANCES_KM_LONG[1:], tgw_long[1:], 'k-', linewidth=1, alpha=0.5,
             label='TGW upper bound')
ax3.semilogy(DISTANCES_KM_LONG[1:], plob_long[1:], 'k--', linewidth=1, alpha=0.5,
             label='PLOB lower bound')

mask_dv_long = dv_rates_long > 1e-12
ax3.semilogy(DISTANCES_KM_LONG[mask_dv_long], dv_rates_long[mask_dv_long],
             'b-', linewidth=2.5, label='DV: Decoy BB84 (TNO)')

ax3.set_xlabel('Distance (km)', fontsize=12)
ax3.set_ylabel('Secret key rate (bits per channel use)', fontsize=12)
ax3.set_title('DV Extended Range — Asymptotic Decoy BB84', fontsize=12)
ax3.set_xlim(0, 300)
ax3.set_ylim(1e-12, 10)
ax3.legend(fontsize=10, loc='lower left')
ax3.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('figure3_dv_extended.png', dpi=150, bbox_inches='tight')
print("Saved: figure3_dv_extended.png")

for i, L in enumerate(DISTANCES_KM_LONG):
    if dv_rates_long[i] < 1e-10 and dv_rates_long[max(i-1, 0)] >= 1e-10:
        print(f"DV reaches ~1e-10 bits/pulse at: {L} km")
        break


# ============================================================
# PRINT TABULAR RESULTS
# ============================================================
print(f"\n{'L(km)':>6} {'DV b/pulse':>11} {'CV-het b/sym':>13} {'CV-hom b/sym':>13} "
      f"{'DV bps':>11} {'CV-het bps':>12} {'CV-hom bps':>12}")
for i, L in enumerate(DISTANCES_KM):
    if L % 5 == 0:
        print(f"{L:6d} {dv_rates[i]:11.3e} {cv_het_rates[i]:13.3e} {cv_hom_rates[i]:13.3e} "
              f"{dv_rates[i]*DV_CLOCK_HZ:11.3e} {cv_het_rates[i]*CV_CLOCK_HZ:12.3e} "
              f"{cv_hom_rates[i]*CV_CLOCK_HZ:12.3e}")