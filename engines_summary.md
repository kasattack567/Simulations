# Simulation engines: how DV and CV are modelled

Two independent, published key-rate engines are used, one per protocol, held at a
matched level of rigour so the comparison is fair.

| | DV | CV |
|---|---|---|
| Engine | TNO `qkd_key_rate` | `qosst-skr` (QOSST project) |
| Protocol | Decoy-state BB84 | GG02, Gaussian-modulated coherent states |
| Class used | `BB84FullyAsymptoticKeyRateEstimate` | `GaussianTrusted{Homodyne,Heterodyne}Asymptotic` |
| Security | Asymptotic | Asymptotic, collective attacks, reverse reconciliation |
| Output | bits per pulse | bits per symbol |

---

## DV — TNO decoy-state BB84

The engine takes detector efficiency, dark-count rate per gate, polarisation
drift and channel attenuation, then **optimises the mean photon number μ** and
returns the asymptotic secret key rate.

Two implementation details worth recording:

1. **QBER is injected as polarisation drift.** The asymptotic rate path derives
   the channel error rate from dark counts and polarisation drift; it does *not*
   read the `error_detector` argument. The target QBER is therefore encoded as
   `polarisation_drift = arcsin(√QBER)`, verified against the engine's own
   `compute_gain_and_error_rate` to better than 10⁻³. Without this the model runs
   at near-zero misalignment error and is optimistic.

2. **Reconciliation efficiency is applied afterwards.** The engine hardcodes
   error correction at the Shannon limit (β = 1). β = 0.95 is applied by scaling
   the error-correction cost:
   `rate = rate_ideal − (1/β − 1) · Q · h(E)`.
   This is the standard `f(E)` scaling and is exact, not an approximation. One
   caveat: μ remains optimised at β = 1, which makes DV slightly *pessimistic* —
   a conservative direction for a comparative study.

## CV — qosst-skr GG02

The engine takes modulation variance Va, channel transmittance T, excess noise ξ,
detector efficiency η, electronic noise v_el and reconciliation efficiency β, and
returns the asymptotic rate under collective attacks. Va is **optimised at every
distance** (capped at 10 SNU to reflect deployed modulator headroom).

**Trusted-detector model.** Detector loss and electronic noise are treated as
calibrated and are not attributed to Eve. This is the standard realistic model
(Lodewyck et al. 2007; Appendix of Wang et al. 2019) and is what deployed systems
assume. Because v_el is a separate argument, it must not also be folded into ξ —
doing so double-counts detector noise.

**Excess noise follows Wang et al. 2019, Eq. 12:**

> ξ_r(T) = (ε_a + ε_l) + ε_b / (ηT),  with ε_a + ε_l = 0.005 and ε_b = 0.0005 SNU

Bob-side noise (ε_b) arises *after* the channel, so it does not attenuate; when
referred to the channel input — the security convention, and what qosst-skr's ξ
argument expects — it is amplified by 1/(ηT) and therefore grows with distance.
This is why CV has a hard distance ceiling that a constant ξ does not reproduce.

---

## Why the two are comparable

- **Same tier of rigour:** both asymptotic. No finite-key on either side.
- **Shared parameters where the quantity is physically the same:** fibre
  attenuation 0.20 dB/km, reconciliation efficiency β = 0.95.
- **Detector efficiency matched by tier**, not by number: DV 0.65 (deployed
  SNSPD), CV 0.60 (deployed homodyne).
- **Per-hop evaluation in the network model:** relayed links are evaluated at each
  hop's own length, which is essential because ξ_r depends on T.
- **Clock rates applied only at the reporting stage** (DV 1 GHz, CV 100 MHz) to
  convert bits/symbol into bits/second.

## Validation

| Check | Result |
|---|---|
| CV reproduces Wang et al. Fig. 1 (ξ = 0.01) | 404.6 km vs paper's ~400 km at 10⁻⁹ |
| CV optimal Va converges (β = 0.95) | 3.712 vs paper's 3.71 |
| CV reach vs commercial hardware | 94 km / 18.8 dB vs LuxQuanta NOVA LQ Gen-2 spec of 100 km / 20 dB |
| Excess-noise model vs measured prototype | tracks Wang Table 2 across 0–100 km |

## Headline results

- DV/CV crossover: **~50 km** at β = 0.95. Soft — spans 27–57 km across the
  deployed β range (0.90–0.96), so it should be quoted as a range.
- CV reach **94 km**; DV reach **279 km**. Both robust to β (±5%).
- Trusted-node budget: CV needs **2.8–3.8×** more relays than DV across all six
  topologies, when relay spacing is sized to a 10 Mbps service target.

## Known limitations

- **No finite-key analysis.** Note this is a deliberate symmetry choice, not
  purely a tooling limit: TNO *does* provide `BB84FiniteKeyRateEstimate`, but
  qosst-skr ships asymptotic classes only. Using finite-key on DV alone would
  bias the comparison, and CV's parameter-estimation penalty is generally the
  harsher of the two — so the asymptotic choice mildly favours CV. Stated as a
  limitation with its direction.
- **No composable security or coherent-attack analysis** for CV; security is
  against collective attacks.
- **Trusted-detector assumption for CV.** qosst-skr also ships
  `GaussianUntrustedHomodyneAsymptotic`, so the paranoid (detector noise
  attributed to Eve) variant is available as a robustness check for homodyne.
  Not currently used.
- **No per-node insertion loss** in the network model. Real deployments add mux,
  demux and switch losses at each node; the Madrid MadQCI trial reports paths made
  unfeasible by exactly this. Current relay budgets are therefore slightly
  optimistic for both protocols.
