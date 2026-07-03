"""
Group 3 — Key rate vs fibre attenuation alpha, at 15 km. Shared identical
parameter -> shared axis, all three curves (DV, CV-het, CV-hom).

JUSTIFICATION FOR THE SWEEP (per supervisor): a sweep around 0.2 dB/km is not
physically meaningful — essentially all deployed silica fibre worldwide sits at
~0.2, and silica loss has been ~flat for 40 years (0.154 dB/km in 1985 ->
0.1396 in 2024). The motivation comes from Hollow-Core Fibre (HCF): the air core
gives both near-vacuum light speed and record-low loss. Petrovich et al. (Nature
Photonics 19, 1203, 2025) demonstrated 0.091 dB/km at 1550 nm, with models
projecting toward 0.01 dB/km. The sweep 0.01-0.25 therefore spans
[HCF projected -> HCF demonstrated -> conventional silica].

Deterministic engines -> no error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, shared_plot3, arg_output_dir,
                         RANGES, BAND, L_KM, BAND_SH)

print(f"[3 fibre atten] alpha 0.01-0.25 dB/km (HCF-motivated) | L={L_KM} km")
a = RANGES['alpha']
dv     = np.array([dv_rate(alpha=x) for x in a])
cv_het = np.array([cv_rate(alpha=x, detection="heterodyne") for x in a])
cv_hom = np.array([cv_rate(alpha=x, detection="homodyne") for x in a])

shared_plot3(a, dv, cv_het, cv_hom,
             'Fibre attenuation  $\\alpha$  (dB/km)',
             f'Key rate vs fibre attenuation   ($L = {L_KM:.0f}$ km)',
             's3_fibre_atten.png',
             bands=[(BAND['alpha'], BAND_SH, 'HCF (0.091) → silica (0.20)')],
             output_dir=arg_output_dir())
