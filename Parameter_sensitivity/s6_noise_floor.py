"""
Group 6 (new) — CV excess-noise floor: key rate vs the constant Alice+fibre term
eps_a + eps_l, at fixed eps_b. Complements s5, which sweeps the Bob-side term.

Under Wang Eq. 12 the two CV noise terms behave completely differently: the floor
eps_a + eps_l shifts the whole rate-distance curve down without changing its
shape, while eps_b steepens the roll-off and sets the reach. Separating them is
the point of the decomposition, so both deserve a sweep. There is no DV analogue
of a distance-independent channel-noise floor, so this is a CV-only figure.

Deterministic engines -> no error bars.
"""
import numpy as np
from sens_common import (cv_rate, cv_xi_input, T_of_L, shared_plot, arg_common,
                         RANGES, BAND, BAND_SH, CV_XI_B, CV_COLOR)

args = arg_common()
L = args.distance
T = T_of_L(L)
xa = RANGES['cv_xi_al']

print(f"[6 noise floor] CV eps_a+eps_l {xa[0]:.3f}-{xa[-1]:.3f} SNU"
      f" (eps_b fixed at {CV_XI_B:.0e}) | L={L:g} km")
print(f"    implied total xi_r at {L:g} km: "
      f"{cv_xi_input(T, xi_al=xa[0]):.4f} - {cv_xi_input(T, xi_al=xa[-1]):.4f} SNU")

cv_het = np.array([cv_rate(L_km=L, xi_al=x, detection="heterodyne") for x in xa])
cv_hom = np.array([cv_rate(L_km=L, xi_al=x, detection="homodyne") for x in xa])

# Both curves are CV here (heterodyne and homodyne), so the labels, colours and
# markers must be overridden — shared_plot's defaults assume one DV and one CV.
shared_plot(xa, cv_het, cv_hom,
            'CV channel noise floor  $\\epsilon_a + \\epsilon_l$  (SNU)',
            f'CV key rate vs excess-noise floor   ($L = {L:.0f}$ km, '
            f'$\\epsilon_b = {CV_XI_B:.0e}$ SNU)',
            's6_noise_floor.png',
            band=BAND['cv_xi_al'], band_label='Deployed floor (0.003–0.008 SNU)',
            labels=('CV — GG02 heterodyne', 'CV — GG02 homodyne'),
            colors=(CV_COLOR, '#7d2d8c'), markers=('-s', '--D'),
            legend_loc='upper right', output_dir=args.output_dir)