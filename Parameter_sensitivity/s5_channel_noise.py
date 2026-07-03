"""
Group 5 — Channel noise: DV QBER vs CV excess noise, 15 km. Twin x-axis.
CV xi applies to both schemes -> CV-het and CV-hom both shown. DV swept to its
~11% BB84 cutoff. Crossings across the two axes are NOT physically meaningful.
No error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, twin_plot3, arg_output_dir,
                         RANGES, BAND, L_KM)

print(f"[5 channel noise] DV QBER 0.1-11% / CV xi 0.001-0.08 | L={L_KM} km")
dv     = np.array([dv_rate(qber=q) for q in RANGES['dv_qber']])
cv_het = np.array([cv_rate(xi_bob=x, detection="heterodyne") for x in RANGES['cv_xi']])
cv_hom = np.array([cv_rate(xi_bob=x, detection="homodyne") for x in RANGES['cv_xi']])

twin_plot3(RANGES['dv_qber']*100, dv, RANGES['cv_xi'], cv_het, cv_hom,
           'DV QBER  (%)', 'CV excess noise  $\\xi$  (SNU)',
           f'Key rate vs channel noise   ($L = {L_KM:.0f}$ km)',
           's5_channel_noise.png',
           band_dv=(BAND['dv_qber'][0]*100, BAND['dv_qber'][1]*100),
           band_cv=BAND['cv_xi'],
           band_dv_label='DV deployed (0.5–2.1%)',
           band_cv_label='CV deployed (0.01–0.03 SNU)',
           output_dir=arg_output_dir())
