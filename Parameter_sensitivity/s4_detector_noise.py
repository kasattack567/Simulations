"""
Group 4 — Detector noise: DV dark count vs CV electronic noise. Twin x-axis
(different units). CV v_el applies to both detection schemes, so both CV-het and
CV-hom are shown on the top axis. DV flat: dark counts negligible at this
distance (corroborated by independent NetSquid simulation).

v_el is the one CV noise term that does NOT enter the Wang Eq. 12 excess-noise
expression — qosst-skr's trusted-detector engine takes it separately, so there is
no double-counting between v_el and xi_r.

Crossings across the two axes are NOT physically meaningful. No error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, twin_plot3, arg_common,
                         RANGES, BAND, DARK_CPS_TO_PERGATE)

args = arg_common()
L = args.distance
print(f"[4 detector noise] DV dark 1-1e5 cps / CV v_el 0.0-0.5 SNU | L={L:g} km")

dv     = np.array([dv_rate(L_km=L, dark_pergate=c*DARK_CPS_TO_PERGATE)
                   for c in RANGES['dv_dark']])
cv_het = np.array([cv_rate(L_km=L, vel=v, detection="heterodyne")
                   for v in RANGES['cv_vel']])
cv_hom = np.array([cv_rate(L_km=L, vel=v, detection="homodyne")
                   for v in RANGES['cv_vel']])

twin_plot3(RANGES['dv_dark'], dv, RANGES['cv_vel'], cv_het, cv_hom,
           'DV dark-count rate  (cps)', 'CV electronic noise  $v_{el}$  (SNU)',
           f'Key rate vs detector noise   ($L = {L:.0f}$ km)',
           's4_detector_noise.png',
           band_dv=BAND['dv_dark'], band_cv=BAND['cv_vel'],
           band_dv_label='DV deployed (1–1000 cps)',
           band_cv_label='CV deployed (0.015–0.11 SNU)',
           logx_dv=True, legend_loc='lower left', output_dir=args.output_dir)