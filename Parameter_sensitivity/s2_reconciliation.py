"""
Group 2 — Key rate vs reconciliation efficiency.

The two protocols do NOT share this parameter. CV's beta multiplies the whole of
I_AB; DV's f_EC scales only the error-correction leakage. Mueller 2024 Eq. 6
relates them by beta = (1 - f h(e))/(1 - h(e)), NOT by f = 1/beta, and at DV's
operating QBER of 0.5% the two are numerically miles apart: f = 1.036 is an
effective beta of 0.9975, and a true beta of 0.95 would need f = 2.05.

Plotting them on one axis, as an earlier version of this figure did, therefore
compares two different quantities. The figure now uses TWIN X-AXES: CV beta on
the bottom, DV f_EC on the top, each swept over its own deployed band. Curve
crossings between the blue and red lines are not physically meaningful, exactly
as in the twin-axis panels of s4 and s5.

The finding survives the correction and is in fact its clearest illustration.
DV is nearly flat across its whole band, because leakage is a small part of its
rate; CV climbs steeply, because beta gates everything. Deterministic engines ->
no error bars.

Deployed bands (see BAND in sens_common for full sourcing):
  CV beta   0.90-0.96   Hajomer 2024 (0.9091 at FER=0) to Wang 2017 (0.964)
  DV f_EC   1.036-1.20  Mueller 2024 Cascade to Xu et al. RMP 2020 Cascade
"""
import numpy as np
from sens_common import (band_label, dv_rate, cv_rate, twin_plot3, arg_common,
                         RANGES, BAND, BAND_DV, BAND_CV)

args = arg_common()
L = args.distance
print(f"[2 reconciliation] CV beta 0.80-1.00 | DV f_EC 1.00-1.40 | L={L:g} km")
b = RANGES['beta']          # CV reconciliation efficiency
f = RANGES['dv_f_ec']       # DV error-correction efficiency

dv     = np.array([dv_rate(L_km=L, f_ec=x) for x in f])
cv_het = np.array([cv_rate(L_km=L, beta=x, detection="heterodyne") for x in b])
cv_hom = np.array([cv_rate(L_km=L, beta=x, detection="homodyne") for x in b])

twin_plot3(f, dv, b, cv_het, cv_hom,
           'DV error-correction efficiency  $f_{EC}$',
           'CV reconciliation efficiency  $\\beta$',
           f'Key rate vs reconciliation   ($L = {L:.0f}$ km)',
           's2_reconciliation.png',
           band_dv=BAND['dv_f_ec'], band_cv=BAND['beta'],
           band_dv_label=band_label('dv_f_ec'),
           band_cv_label=band_label('beta'),
           output_dir=args.output_dir)