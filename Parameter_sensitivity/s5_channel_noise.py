"""
Group 5 — Channel noise: DV QBER vs CV Bob-side excess noise eps_b. Twin x-axis.

WHAT CHANGED AND WHY. This script previously swept a single total CV excess noise
xi. Under Wang et al. 2019 Eq. 12 the input-referred excess noise is

    xi_r(T) = (eps_a + eps_l) + eps_b / (eta * T)

so "total xi" is not a free parameter — it is an output that depends on distance.
The sweepable physical quantities are eps_b (Bob-side, arises after the channel,
amplified by 1/(eta*T) when referred to the input) and the constant Alice+fibre
floor eps_a + eps_l. eps_b is the informative axis because it alone controls how
fast CV degrades with distance, and therefore sets the CV reach and the trusted-
node budget. The floor is swept as a secondary curve set (see --floor).

The x-axis is eps_b in SNU. The printout reports the implied total xi_r at the
test distance so the figure can be read against measured field values (Wang
Table 2: 0.005 SNU at 0 km rising to ~0.098 SNU at 100 km).

CV eps_b applies to both detection schemes -> CV-het and CV-hom both shown.
DV QBER swept to its ~11% BB84 cutoff. Crossings across the two axes are NOT
physically meaningful. No error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, cv_xi_input, T_of_L, twin_plot3,
                         arg_common, RANGES, BAND)

args = arg_common()
L = args.distance
T = T_of_L(L)

xb = RANGES['cv_xi_b']
print(f"[5 channel noise] DV QBER 0.1-11% / CV eps_b {xb[0]:.0e}-{xb[-1]:.0e} SNU"
      f" | L={L:g} km")
print(f"    implied total xi_r at {L:g} km: "
      f"{cv_xi_input(T, xi_b=xb[0]):.4f} - {cv_xi_input(T, xi_b=xb[-1]):.4f} SNU")
print(f"    deployed eps_b band {BAND['cv_xi_b']} -> xi_r "
      f"{cv_xi_input(T, xi_b=BAND['cv_xi_b'][0]):.4f} - "
      f"{cv_xi_input(T, xi_b=BAND['cv_xi_b'][1]):.4f} SNU")

dv     = np.array([dv_rate(L_km=L, qber=q) for q in RANGES['dv_qber']])
cv_het = np.array([cv_rate(L_km=L, xi_b=x, detection="heterodyne") for x in xb])
cv_hom = np.array([cv_rate(L_km=L, xi_b=x, detection="homodyne") for x in xb])

twin_plot3(RANGES['dv_qber']*100, dv, xb*1e3, cv_het, cv_hom,
           'DV QBER  (%)',
           'CV Bob-side excess noise  $\\epsilon_b$  ($10^{-3}$ SNU)',
           f'Key rate vs channel noise   ($L = {L:.0f}$ km)',
           's5_channel_noise.png',
           band_dv=(BAND['dv_qber'][0]*100, BAND['dv_qber'][1]*100),
           band_cv=(BAND['cv_xi_b'][0]*1e3, BAND['cv_xi_b'][1]*1e3),
           band_dv_label='DV deployed (0.5–2.1%)',
           band_cv_label='CV deployed (0.5–2.0 $\\times10^{-3}$ SNU)',
           legend_loc='lower left', output_dir=args.output_dir)