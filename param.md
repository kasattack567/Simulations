# Evaluation Parameters — DV vs CV QKD (Deployed Baseline)

Parameters reflect QKD as fielded today: deployed/field-tested fibre, current commercial-grade detectors, hardware that has run outside the lab. Both protocols sit at the same tier — DV with cooled SNSPD receivers (as deployed QKD networks use), CV with the homodyne receivers real CV deployments use. The same locked baseline is used throughout: point-to-point (0–320 km) and the network-topology study at native link lengths. CV excess noise follows Wang et al. 2019 Eq. 12 (channel-input referred): `xi_r(T) = (eps_a + eps_l) + eps_b / (eta * T)`, so total excess noise is an output of the model, not an input. TNO's BB84 engine is at the Shannon limit (f_EC = 1); a realistic f_EC = 1.036 is applied post-hoc in `dv_skr`. Intrinsic QBER is encoded via `polarization_drift = arcsin(sqrt(QBER))` because TNO's `error_detector` is inert in the asymptotic rate path.

## Locked parameter set

### Shared
| Parameter | Value | Primary reference |
|---|---|---|
| Fibre attenuation alpha | 0.20 dB/km | SMF-28 @ 1550 nm |

### DV (decoy BB84, SNSPD receiver)
| Parameter | Value | Primary reference |
|---|---|---|
| Detector SDE | 0.92 | Clason et al. 2026 (arXiv:2606.06107), ID281 SNSPD |
| Receiver insertion loss | 2.0 dB | cf. Kelsey 2026 (same detector family) |
| Total receiver efficiency | 0.580 | = 0.92 × 10^(−0.2) |
| Dark count (per gate) | 3e-8 | Clason 2026, 30 cps at 1 ns / 1 GHz |
| Intrinsic QBER | 0.005 | Clason 2026 (0.5 ± 0.2%) |
| f_EC | 1.036 | Mueller et al. 2024 (arXiv:2408.15758), Cascade on live industrial system |
| Clock | 1 GHz | Modern DV standard (Clavis XGR) |

### CV (GG02, homodyne + heterodyne)
| Parameter | Value | Primary reference |
|---|---|---|
| Detector efficiency eta | 0.60 | Lodewyck 2007 (0.606), Zhang 2020 (0.613) |
| Electronic noise v_el | 0.10 SNU | Zhang 2020, deployment-grade homodyne |
| Excess noise, Alice + fibre (eps_a + eps_l) | 0.005 SNU | Wang et al. 2019, Eq. 12 / Table 2 |
| Excess noise, Bob-side (eps_b) | 0.0005 SNU | Wang et al. 2019 prototype calibration |
| Reconciliation beta | 0.95 | Zhang 2020 Table I |
| Modulation variance Va | optimised per distance | qosst-skr Va sweep |
| Clock | 100 MHz | QOSST-class deployed CV |

## Primary sources

- Clason et al., "Field trial over 303 km deployed fibre" (arXiv:2606.06107, 2026) — DV detector, dark count, QBER.
- Mueller et al., "Cascade on a live industrial system" (arXiv:2408.15758, 2024) — DV f_EC.
- Wang et al., "High key rate CV-QKD with a real local oscillator," Opt. Express 27, 13372 (2019) — CV excess-noise model (Eq. 12), Table 2 calibration.
- Zhang et al., "Long-Distance CV-QKD over 202.81 km," PRL 125, 010502 (2020) — CV eta, v_el, beta.
- Lodewyck et al., "CV-QKD over 25 km all-fibre," PRA 76, 042305 (2007) — CV homodyne anchor.
- Pietri et al., "QOSST," Quantum 8, 1575 (2024) — qosst-skr engine and xi convention.
- Kelsey 2026 — receiver-node optical budget for ID281-class SNSPD systems.
- LuxQuanta NOVA LQ Gen-2 — commercial CV-QKD (100 km / 20 dB spec), independent hardware check on eps_b.