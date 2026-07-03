"""
Group 4 — Detector noise: DV dark count vs CV electronic noise, 15 km.
Twin x-axis (different units). CV v_el applies to both detection schemes, so
both CV-het and CV-hom are shown on the top axis. DV flat: dark counts negligible
at metro distance (corroborated by independent NetSquid simulation).
Crossings across the two axes are NOT physically meaningful. No error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, twin_plot3, arg_output_dir,
                         RANGES, BAND, L_KM, DARK_CPS_TO_PERGATE)

print(f"[4 detector noise] DV dark 1-1e5 cps / CV v_el 0.0-0.5 | L={L_KM} km")
dv     = np.array([dv_rate(dark_pergate=c*DARK_CPS_TO_PERGATE) for c in RANGES['dv_dark']])
cv_het = np.array([cv_rate(vel=v, detection="heterodyne") for v in RANGES['cv_vel']])
cv_hom = np.array([cv_rate(vel=v, detection="homodyne") for v in RANGES['cv_vel']])

twin_plot3(RANGES['dv_dark'], dv, RANGES['cv_vel'], cv_het, cv_hom,
           'DV dark-count rate  (cps)', 'CV electronic noise  $v_{el}$  (SNU)',
           f'Key rate vs detector noise   ($L = {L_KM:.0f}$ km)',
           's4_detector_noise.png',
           band_dv=BAND['dv_dark'], band_cv=BAND['cv_vel'],
           band_dv_label='DV deployed (1–1000 cps)',
           band_cv_label='CV deployed (0.05–0.11 SNU)',
           logx_dv=True, output_dir=arg_output_dir())
