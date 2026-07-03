"""
Group 1 — Key rate vs detector efficiency, at 15 km. Detector efficiency is the
SAME physical quantity for all -> shared axis, three curves (DV, CV-het, CV-hom).
Crossings ARE meaningful here. Both deployed bands shaded (DV SNSPD 0.65-0.93,
CV homodyne 0.60-0.72); DV's band sits higher, reflecting the detector gap.
Deterministic engines -> no error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, shared_plot3, arg_output_dir,
                         BAND, L_KM, BAND_DV, BAND_CV)

print(f"[1 detector eff] eta 0.50-1.00 (shared) | L={L_KM} km")
eta = np.linspace(0.50, 1.00, 40)
dv     = np.array([dv_rate(eta=e) for e in eta])
cv_het = np.array([cv_rate(eta=e, detection="heterodyne") for e in eta])
cv_hom = np.array([cv_rate(eta=e, detection="homodyne") for e in eta])

shared_plot3(eta, dv, cv_het, cv_hom,
             'Detector efficiency  $\\eta$',
             f'Key rate vs detector efficiency   ($L = {L_KM:.0f}$ km)',
             's1_detector_eff.png',
             bands=[(BAND['dv_eta'], BAND_DV, 'DV deployed SNSPD (0.65–0.93)'),
                    (BAND['cv_eta'], BAND_CV, 'CV deployed homodyne (0.60–0.72)')],
             output_dir=arg_output_dir())
