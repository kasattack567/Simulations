"""
Group 1 — Key rate vs detector efficiency. Detector efficiency is the SAME
physical quantity for all -> shared axis, three curves (DV, CV-het, CV-hom).
Crossings ARE meaningful here. Both deployed bands shaded (DV SNSPD 0.65-0.93,
CV homodyne 0.60-0.72); DV's band sits higher, reflecting the detector gap.
Deterministic engines -> no error bars.

NB under Wang Eq. 12 the CV curves respond to eta TWICE: once through detection
efficiency and once through the input-referred excess noise xi_r = eps_a + eps_l
+ eps_b/(eta*T). The CV slope is therefore steeper than detection efficiency
alone would produce. This is physical, not an artefact.
"""
import numpy as np
from sens_common import (band_label, dv_rate, cv_rate, cv_xi_input, T_of_L, shared_plot3,
                         arg_common, BAND, BAND_DV, BAND_CV)

args = arg_common()
L = args.distance
print(f"[1 detector eff] eta 0.50-1.00 (shared) | L={L:g} km")
print(f"    xi_r at eta=0.50 / 1.00: {cv_xi_input(T_of_L(L), eta=0.50):.5f}"
      f" / {cv_xi_input(T_of_L(L), eta=1.00):.5f} SNU")

eta = np.linspace(0.50, 1.00, 40)
dv     = np.array([dv_rate(L_km=L, eta=e) for e in eta])
cv_het = np.array([cv_rate(L_km=L, eta=e, detection="heterodyne") for e in eta])
cv_hom = np.array([cv_rate(L_km=L, eta=e, detection="homodyne") for e in eta])

shared_plot3(eta, dv, cv_het, cv_hom,
             'Detector efficiency  $\\eta$',
             f'Key rate vs detector efficiency   ($L = {L:.0f}$ km)',
             's1_detector_eff.png',
             bands=[(BAND['dv_eta'], BAND_DV, band_label('dv_eta')),
                    (BAND['cv_eta'], BAND_CV, band_label('cv_eta'))],
             output_dir=args.output_dir)