"""
Group 3 — Key rate vs fibre attenuation alpha. Shared identical parameter ->
shared axis, all three curves (DV, CV-het, CV-hom).

JUSTIFICATION FOR THE SWEEP (per supervisor): a sweep around 0.2 dB/km is not
physically meaningful — essentially all deployed silica fibre worldwide sits at
~0.2, and silica loss has been ~flat for 40 years (0.154 dB/km in 1985 ->
0.1396 in 2024). The motivation comes from Hollow-Core Fibre (HCF): the air core
gives both near-vacuum light speed and record-low loss. Petrovich et al. (Nature
Photonics 19, 1203, 2025) demonstrated 0.091 dB/km at 1550 nm, with models
projecting toward 0.01 dB/km. The sweep 0.01-0.25 therefore spans
[HCF projected -> HCF demonstrated -> conventional silica].

NB under Wang Eq. 12 alpha enters the CV rate TWICE: through the channel
transmittance and through xi_r = eps_a + eps_l + eps_b/(eta*T). Lower loss
suppresses excess noise as well as raising transmittance, so HCF helps CV
super-linearly. This is the clearest new result in the sensitivity suite.

Deterministic engines -> no error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, cv_xi_input, T_of_L, shared_plot3,
                         arg_common, RANGES, BAND, BAND_SH)

args = arg_common()
L = args.distance
a = RANGES['alpha']
print(f"[3 fibre atten] alpha {a[0]:.2f}-{a[-1]:.2f} dB/km (HCF-motivated) | L={L:g} km")
print(f"    xi_r at alpha={a[0]:.2f} / {a[-1]:.2f}: "
      f"{cv_xi_input(T_of_L(L, a[0])):.5f} / {cv_xi_input(T_of_L(L, a[-1])):.5f} SNU")

dv     = np.array([dv_rate(L_km=L, alpha=x) for x in a])
cv_het = np.array([cv_rate(L_km=L, alpha=x, detection="heterodyne") for x in a])
cv_hom = np.array([cv_rate(L_km=L, alpha=x, detection="homodyne") for x in a])

shared_plot3(a, dv, cv_het, cv_hom,
             'Fibre attenuation  $\\alpha$  (dB/km)',
             f'Key rate vs fibre attenuation   ($L = {L:.0f}$ km)',
             's3_fibre_atten.png',
             bands=[(BAND['alpha'], BAND_SH, 'HCF (0.091) → silica (0.20)')],
             output_dir=args.output_dir)
