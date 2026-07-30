"""
Group 2 — Key rate vs reconciliation efficiency beta, at 15 km. Shared axis,
three curves. Finding: DV near-insensitive to beta (only scales the small EC
cost); both CV schemes depend on it strongly (beta gates the whole DW rate).
Deterministic engines -> no error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, shared_plot3, arg_output_dir,
                         RANGES, BAND, L_KM, BAND_SH)

print(f"[2 reconciliation] beta 0.80-1.00 (shared) | L={L_KM} km")
b = RANGES['beta']
dv     = np.array([dv_rate(beta=x) for x in b])
cv_het = np.array([cv_rate(beta=x, detection="heterodyne") for x in b])
cv_hom = np.array([cv_rate(beta=x, detection="homodyne") for x in b])

shared_plot3(b, dv, cv_het, cv_hom,
             'Reconciliation efficiency  $\\beta$',
             f'Key rate vs reconciliation efficiency   ($L = {L_KM:.0f}$ km)',
             's2_reconciliation.png',
             bands=[(BAND['beta'], BAND_SH, 'Deployed $\\beta$ (0.90-0.96)')],
             output_dir=arg_output_dir())
