"""
Group 2 — Key rate vs reconciliation efficiency beta, at 15 km. Shared axis,
three curves. Finding: DV near-insensitive to beta (only scales the small EC
cost); both CV schemes depend on it strongly (beta gates the whole DW rate).
Deterministic engines -> no error bars.
"""
import numpy as np
from sens_common import (dv_rate, cv_rate, shared_plot3, arg_common,
                         RANGES, BAND, BAND_SH)

args = arg_common()
L = args.distance
print(f"[2 reconciliation] beta 0.80-1.00 (shared) | L={L:g} km")
b = RANGES['beta']
dv     = np.array([dv_rate(L_km=L, beta=x) for x in b])
cv_het = np.array([cv_rate(L_km=L, beta=x, detection="heterodyne") for x in b])
cv_hom = np.array([cv_rate(L_km=L, beta=x, detection="homodyne") for x in b])

shared_plot3(b, dv, cv_het, cv_hom,
             'Reconciliation efficiency  $\\beta$',
             f'Key rate vs reconciliation efficiency   ($L = {L:.0f}$ km)',
             's2_reconciliation.png',
             bands=[(BAND['beta'], BAND_SH, 'Deployed $\\beta$ (0.90-0.96)')],
             output_dir=args.output_dir)
