"""
Final point-to-point DV vs CV comparison — thesis figures.

CV computed via qosst-skr (consistent framework for both detection schemes):
  - GaussianTrustedHeterodyneAsymptotic
  - GaussianTrustedHomodyneAsymptotic
DV via TNO asymptotic decoy BB84.

EVALUATION PARAMETERS: deployed / commercial-grade hardware (see param.md for
full sourcing). DV uses cooled-SNSPD deployed-network numbers (Tang 2016 Hefei);
CV uses deployed homodyne numbers (Lodewyck 2007 / Zhang 2020). This puts the
two protocols on comparable detector footing (eta 0.65 vs 0.60).

CV EXCESS-NOISE MODEL — Wang et al., Opt. Express 27, 13372 (2019), Sec. 5:
    xi_r(T) = (eps_a + eps_l) + eps_b / (eta * T)          [Eq. 12]
Alice-side + fibre noise (eps_a + eps_l) is ~constant; Bob-side measurement
noise (eps_b) arises AFTER the channel, so referred to the channel input (the
security convention, and what qosst-skr's `xi` argument expects) it is
amplified by 1/(eta*T) and grows with distance. Values are Wang's calibration
against their measured prototype (Table 2 / Fig. 8). This model reproduces the
~100 km realistic CV limit; a constant-xi idealised reference curve (Wang
Fig. 8's fixed-noise line) is plotted alongside for the idealised-vs-realistic
gap. Both engines are ASYMPTOTIC (matched tier); finite-size is out of scope.

DV reconciliation: TNO's BB84FullyAsymptoticKeyRateEstimate hardcodes the
error-correction term at the Shannon limit (beta=1, ideal) with no argument to
change it. DV_BETA = 0.95 is therefore applied POST-HOC in dv_skr (matching CV's
tier and the locked param.md set), by scaling the EC cost term. See dv_skr.

Figure 2: bits per channel use (protocol comparison), DV + CV-het + CV-hom
          + idealised-CV reference
Figure 4: bits/s at realistic clocks (deployment view)

Both figures span 0-320 km, i.e. beyond the DV distance limit under this
parameter set, so that the full reach of each protocol is visible. (The
idealised CV reference extends past the axis, as in Wang Fig. 1.)
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
ALPHA_DB_KM = 0.20         # SMF-28 at 1550 nm; deployed fibre 0.19-0.21 dB/km

# DV (TNO decoy BB84) — deployed SNSPD receiver (Tang 2016 Hefei network)
DV_EFFICIENCY     = 0.65   # SNSPD system efficiency (Tang 2016: 0.64-0.66)
DV_DARK_COUNT     = 1e-7   # per gate: 100 cps at ~1 ns window / 1 GHz clock
DV_QBER           = 0.005  # well-aligned source; field QBER <1%
DV_BETA           = 0.95   # reconciliation, matched to CV tier. NB: TNO engine is
                           # ideal (beta=1); 0.95 applied post-hoc in dv_skr.

# CV (qosst-skr GG02) — deployed homodyne (Lodewyck 2007 / Zhang 2020)
CV_ETA            = 0.60   # homodyne efficiency (Lodewyck 0.606, Zhang 0.613)
CV_VEL            = 0.10   # electronic noise, deployment-grade homodyne (SNU)
CV_BETA           = 0.95   # reconciliation efficiency, deployed CV standard

# CV excess noise — Wang et al. 2019, Eq. 12 (channel-input referred):
#   xi_r(T) = CV_XI_AL + CV_XI_B / (CV_ETA * T)
# Calibrated by Wang against their prototype's measured excess noise
# (Table 2; model line in Fig. 8, "agree approximately"). eps_b does not
# attenuate with the channel, hence the 1/(eta*T) amplification at input.
CV_XI_AL          = 0.005   # eps_a + eps_l: Alice + fibre channel (SNU)
CV_XI_B           = 0.0005  # eps_b: Bob-side measurement noise (SNU)

# Idealised CV reference (Wang Fig. 8 fixed-noise line): constant channel-input
# xi with NO Bob-side growth. Plotted thin on Figure 2 to show the idealised-
# vs-realistic gap; not part of the deployed baseline.
CV_XI_IDEAL       = 0.01    # SNU, constant
INCLUDE_IDEALISED = True

# Clock rates for Figure 4 (realistic deployment)
DV_CLOCK_HZ       = 1e9     # 1 GHz, typical modern DV system
CV_CLOCK_HZ       = 100e6   # 100 MHz, representative deployed CV symbol rate

# Distance range
DISTANCES_KM = np.arange(0, 321, 0.5)

# Plotting floor for bits/channel-use (Figure 2). Rates below this are treated as
# no key and simply not drawn, so the DV tail is visible out to its cutoff.
RATE_FLOOR = 1e-10

# Va search bounds (modulation variance, SNU). Kept unconstrained: the GG02
# asymptotic optimum exceeds ~10 SNU only at sub-km distance (not physically
# achievable on real modulators, but outside the range of interest). The
# crossover and all points >=5 km are unaffected; 0-1 km CV is upper-idealised.
VA_LO, VA_HI = 1e-2, 100.0


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

#===========================================================
# FIGURE A1: repeaterless (PLOB) bound check  [Appendix]
# ============================================================
# Eq. (2): R <= -log2(1 - T). Diverges at L = 0 (T = 1), so L = 0 is skipped.
with np.errstate(divide='ignore'):
    plob_bound = -np.log2(1.0 - T_of_L(DISTANCES_KM))

print("\nRepeaterless (PLOB) bound check:")
for name, rates in [("DV decoy BB84", dv_rates),
                    ("CV heterodyne", cv_het_rates),
                    ("CV homodyne", cv_hom_rates)]:
    m = (rates > RATE_FLOOR) & np.isfinite(plob_bound)
    ratio = plob_bound[m] / rates[m]
    j = int(np.argmin(ratio))
    n_viol = int(np.sum(rates[m] > plob_bound[m]))
    print(f"  {name:14s} closest approach: bound/rate = {ratio[j]:6.2f} "
          f"at {DISTANCES_KM[m][j]:5.1f} km | violations: {n_viol}")

figA, (axA1, axA2) = plt.subplots(1, 2, figsize=(13, 5.5))

# (a) DV
axA1.semilogy(DISTANCES_KM[1:], plob_bound[1:], 'k:', linewidth=2.0,
              label=r'PLOB bound $-\log_2(1-T)$')
axA1.semilogy(DISTANCES_KM[mask_dv], dv_rates[mask_dv], 'b-', linewidth=2.5,
              label='DV: Decoy BB84 (TNO)')
axA1.set_title('(a) DV: decoy-state BB84', fontsize=12)
axA1.set_xlim(0, 320)

# (b) CV
axA2.semilogy(DISTANCES_KM[1:], plob_bound[1:], 'k:', linewidth=2.0,
              label=r'PLOB bound $-\log_2(1-T)$')
axA2.semilogy(DISTANCES_KM[mask_het], cv_het_rates[mask_het], 'r-', linewidth=2.5,
              label='CV: GG02 Heterodyne (qosst-skr)')
axA2.semilogy(DISTANCES_KM[mask_hom], cv_hom_rates[mask_hom], 'm--', linewidth=2.5,
              label='CV: GG02 Homodyne (qosst-skr)')
axA2.set_title('(b) CV: GG02', fontsize=12)
axA2.set_xlim(0, 100)

for ax in (axA1, axA2):
    ax.set_xlabel('Distance (km)', fontsize=12)
    ax.set_ylabel('Secret key rate (bits per channel use)', fontsize=12)
    ax.set_ylim(RATE_FLOOR, 10)
    ax.legend(fontsize=10, loc='lower left')
    ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('appendix_plob_check.png', dpi=150, bbox_inches='tight')
print("Saved: appendix_plob_check.png")

# ============================================================
# PRINT TABULAR RESULTS
# ============================================================
print(f"\n{'L(km)':>6} {'xi_r':>9} {'DV b/pulse':>11} {'CV-het b/sym':>13} "
      f"{'CV-hom b/sym':>13} {'DV bps':>11} {'CV-het bps':>12} {'CV-hom bps':>12}")
for i, L in enumerate(DISTANCES_KM):
    if L % 5 == 0:
        print(f"{L:6d} {cv_xi_input(T_of_L(L)):>9.4g} {dv_rates[i]:11.3e} "
              f"{cv_het_rates[i]:13.3e} {cv_hom_rates[i]:13.3e} "
              f"{dv_rates[i]*DV_CLOCK_HZ:11.3e} {cv_het_rates[i]*CV_CLOCK_HZ:12.3e} "
              f"{cv_hom_rates[i]*CV_CLOCK_HZ:12.3e}")