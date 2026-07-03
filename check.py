"""
Calibrate the CV excess-noise convention against QOSST's OWN documentation example.

The QOSST docs give a reference example (0-20 km) with:
    xi_bob = 0.02, eta = 0.65, Vel = 0.01, beta = 0.95, Va = 5, alpha = 0.2 dB/km
and initialise the channel as:
    GaussianChannel(T, xi_bob / (T * eta))      # <-- the convention, NO factor of 2

This script computes the asymptotic SKR at those exact parameters two ways:
  (A) directly via qosst-skr's GaussianTrustedHeterodyneAsymptotic.skr, using
      xi_input = xi_bob / (T * eta)   [QOSST example convention]
  (B) the same but with the extra factor of 2 (the old suspect form)
so you can see which one behaves sensibly and how much the factor of 2 changes things.

If qosst-sim is installed, it ALSO runs QOSST's own asymptotic calculator on the
identical example as the ground-truth reference, and prints the comparison.

Run in the qkd-env.
"""
import numpy as np
from qosst_skr.gaussian_trusted_heterodyne_asymptotic import (
    GaussianTrustedHeterodyneAsymptotic as HET)

# --- QOSST documentation example parameters (verbatim) ---
XI_BOB = 0.02
ETA    = 0.65
VEL    = 0.01
BETA   = 0.95
VA     = 5.0
ALPHA  = 0.2   # dB/km

def T_of(dist_km):
    return 10 ** (-ALPHA * dist_km / 10.0)

def skr_at(dist_km, factor):
    """SKR via qosst-skr with xi_input = factor * xi_bob / (T*eta)."""
    T = T_of(dist_km)
    xi = factor * XI_BOB / (T * ETA)
    try:
        return HET.skr(Va=VA, T=T, xi=xi, eta=ETA, Vel=VEL, beta=BETA)
    except Exception as e:
        return float('nan')

print("QOSST example params: xi_bob=0.02, eta=0.65, Vel=0.01, beta=0.95, Va=5, "
      "alpha=0.2 dB/km\n")
print(f"{'L(km)':>6} {'xi_input(x1)':>13} {'SKR (x1)':>12} {'xi_input(x2)':>13} {'SKR (x2)':>12}")
for L in range(0, 26, 2):
    T = T_of(L)
    xi1 = 1.0 * XI_BOB / (T*ETA)
    xi2 = 2.0 * XI_BOB / (T*ETA)
    print(f"{L:6d} {xi1:13.4f} {skr_at(L,1.0):12.4e} {xi2:13.4f} {skr_at(L,2.0):12.4e}")

# find cutoffs
def cutoff(factor):
    for L in range(0, 200):
        if skr_at(L, factor) > 1e-9 and skr_at(L+1, factor) <= 1e-9:
            return L+1
    return None
print(f"\nCutoff with x1 (QOSST convention): ~{cutoff(1.0)} km")
print(f"Cutoff with x2 (old suspect form): ~{cutoff(2.0)} km")

# --- Optional: QOSST's own calculator as ground truth, if qosst-sim present ---
print("\n--- Ground truth via qosst-sim (if installed) ---")
try:
    from qosst_sim.modulation.gaussian_qam import GaussianQAM
    from qosst_sim.detector import NoisyHeterodyneDetector
    from qosst_sim.channel import GaussianChannel
    from qosst_sim.simulator.gaussian_channel_asymptotic_calculator import (
        GaussianChannelAsymptoticCalculator)
    mod = GaussianQAM(105, 8, VA, 0.0746269)
    det = NoisyHeterodyneDetector(ETA, VEL)
    print(f"{'L(km)':>6} {'qosst-sim asymptotic SKR':>26} {'our skr (x1)':>14}")
    for L in range(0, 21, 4):
        T = T_of(L)
        ch = GaussianChannel(T, XI_BOB / (T * det.eta))   # QOSST's exact line
        calc = GaussianChannelAsymptoticCalculator(mod, ch, det, BETA)
        print(f"{L:6d} {calc.skr():26.4e} {skr_at(L,1.0):14.4e}")
    print("If the two columns track closely, the x1 (no factor-of-2) convention is confirmed.")
except Exception as e:
    print(f"qosst-sim not available ({e}); rely on the x1 vs x2 comparison above.")