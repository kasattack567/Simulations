# Evaluation Parameters — DV vs CV QKD Comparison (Deployed / Commercial-Grade Baseline)

**Framing adjective: "best demonstrated in deployed QKD systems / current
commercial-grade hardware"**

These parameters reflect QKD systems as fielded today: deployed or field-tested
fibre, current commercial-grade detectors, and hardware that has run outside the
laboratory. The aim is a fair like-for-like comparison at a single consistent
tier: DV evaluated with the cooled SNSPD receivers that deployed QKD networks
field, and CV with the homodyne receivers that real CV deployments use. Neither
protocol is handicapped with cheaper hardware than the other.

The same locked baseline is used throughout the project: in the point-to-point
comparison (simulated over 0-50 km, the range in which both protocols produce
key — CV's reach at this baseline is ~34 km) and in the network-topology study,
where it is applied to real optical network topologies at their native link
lengths. Parameter values are properties of the hardware tier, not of any one
deployment distance, so a single baseline serves both studies.

This corrects an earlier draft that used deployed *SPAD* numbers for DV
(eta ~ 0.2). That is only representative of QKD appliances with low-cost internal
SPADs; when the same field experiments swap in commercial SNSPDs they reach
0.6-0.8 efficiency. For a comparison at the deployed-network tier the SNSPD
figures are the appropriate choice.

> **Framing note (audit trail).** An earlier version of this document was titled
> "Metropolitan Scale" and framed the baseline as "best demonstrated in
> metropolitan deployments". The parameter values and their sources are
> unchanged from that version: every value was already drawn from deployed /
> field-tested / commercial-grade systems, several of which are long-haul
> experiments (Zhang 2020 at 202 km; Huang 2016 at 100 km; the 303 km Swedish
> field trial). The "metropolitan" label described the tier's canonical anchor
> (Tang 2016), not a selection criterion, and has been removed so the baseline
> honestly covers the network study's full range of link lengths. No simulation
> results are affected by this reframing.

---

## 0. FINAL LOCKED PARAMETER SET (one reference per parameter)

These are the fixed values used in all simulations going forward. Each carries a
single primary deployment reference; corroborating sources are in sections 3/6.
Both protocols are fixed at the **same tier** (best-demonstrated deployed /
commercial-grade), per the fairness argument in methodology.md.

### Shared
| Parameter | Symbol | Value | Primary reference |
|---|---|---|---|
| Fibre attenuation | alpha | 0.20 dB/km | Standard SMF-28 @ 1550 nm; deployed fibre 0.19-0.21 (Corning SMF-28 spec; corroborated Tang 2016 installed fibre) |

### DV (decoy BB84, deployed SNSPD receiver)
| Parameter | Symbol | Value | Primary reference |
|---|---|---|---|
| Detector efficiency | eta_d | 0.65 | Tang et al. 2016, PRX 6 011024 — deployed Hefei field network, SNSPD 64-66% |
| Dark count rate | d_c | 100 cps | Tang et al. 2016 (100 Hz); corroborated TFLN KT field test 2025 (120 cps, APL Photonics 10 031301) |
| Intrinsic QBER | e_mis | 0.005 | TFLN KT commercial-node field test 2025 (QBER 0.58%, APL Photonics 10 031301); Avesani/1.5 GHz field ~0.4% |
| Reconciliation | beta_DV | 0.95 | Matched to CV tier; deployed BB84 LDPC reaches ~0.95 (see NB — TNO engine is ideal beta=1, applied post-hoc) |
| Clock | f_DV | 1 GHz | Commercial Clavis XGR (1 GHz); field systems 0.625-1.5 GHz (npj QI 7 2021; arXiv:2602.08908) |

### CV (GG02 homodyne, deployed fibre receiver)
| Parameter | Symbol | Value | Primary reference |
|---|---|---|---|
| Detector efficiency | eta | 0.60 | Lodewyck et al. 2007, PRA 76 042305 (0.606); corroborated Zhang 2020 (0.613) |
| Electronic noise | v_el | 0.10 | Zhang et al. 2020, PRL 125 010502 (deployment-grade homodyne ~0.1 SNU); corroborated arXiv:2305.03419 (0.1) |
| Excess noise | xi_ch (channel) | 0.01 | Channel-only excess noise (detector noise carried separately via Vel in the trusted-detector model). "Typical" channel value 0.01 SNU (Lin/Lütkenhaus dimension-reduction, arXiv:2101.05799; Ghorai et al., arXiv:0807.3751). Deployed field range: 5G fronthaul trusted xi_T = 0.0009–0.0094 (arXiv:2104.04360); Zhang 2020 channel-input 0.0015–0.008 (arXiv:2001.02555) |
| Reconciliation | beta | 0.95 | Deployed CV LDPC standard 0.95 (Shen 2020, Chin. Phys. B 29 040301, 95.42%); 0.98 optimistic ceiling |
| Modulation variance | Va | optimised (cap 10 SNU) | Operating-point dial; deployed Va~3-6 (Jouguet; arXiv:2305.03419 Va=5.8). Numerically optimised per distance |

NB on DV reconciliation: TNO's BB84FullyAsymptoticKeyRateEstimate hardcodes the
error-correction term at the Shannon limit (beta_DV = 1, ideal); it exposes no
reconciliation argument. beta_DV = 0.95 is applied post-hoc to match the CV tier.
Whether to keep DV at engine-default 1.0 or matched 0.95 is a fairness decision
(see methodology.md) — flagged for supervisor.

---

## 1. Summary table

The locked evaluation baseline is given in **Section 0** above (one primary
reference per parameter). Detailed per-parameter sourcing — including the full
list of corroborating experiments and the sensitivity-sweep bands — is in
**Section 3**. This section previously duplicated the Section 0 table; it has been
removed to keep a single authoritative parameter table.

---

## 2. Why deployed SNSPD numbers for DV (not SPAD)

Deployed DV-QKD networks — at every scale from city field networks to the
longest field trials — use **cooled SNSPDs**, now available as turnkey
closed-cycle commercial units (no liquid helium). Reading the actual deployed
and field experiments:

| Experiment | Type | eta_d | DCR | Notes |
|---|---|---|---|---|
| Tang 2016 (Hefei) | Deployed field network, 17-30 km | 0.64-0.66 | 100 Hz | SNSPDs at 2.1 K |
| Intermodal 2025 (arXiv:2602.16680) | Field, IDQ commercial SNSPD | 0.80 | (QBER<1%) | vs internal SPADs at 0.15 |
| Boston (arXiv:1708.00434) | Silicon photonics field test | >0.85 (WSi) / 0.30 (NbN) | ~1000 Hz | Greater Boston field test |
| MMF field (Toshiba 2024) | Multimode-fibre field trial | 0.50 | <10 Hz | Decoy BB84 |
| 200 km MDI | Long-range | 0.56 | 1 Hz | Cold-filtered SNSPD |
| Swedish 303 km field trial (arXiv:2606.06107) | Long-haul deployed fibre | 0.92-0.93 | <=1 Hz | Field SNSPDs; bounds sweep top end |

**eta_d = 0.65 sits squarely in the deployed range** and matches the canonical
Tang 2016 field network — a *system* efficiency (including coupling losses),
which is the correct quantity for a network study; detector-level figures run
higher (commercial IDQ 0.80; Swedish field 0.92+), and the sensitivity band
0.65-0.93 covers them. DCR = 100 cps likewise matches Tang 2016 and the centre
of the 1-1000 Hz deployed range.

The SNSPD choice is consistent across every scale this project studies: the
same receiver class serves the short-link networks where CV also operates and
the long-haul links where only DV survives, so a single DV receiver baseline is
valid throughout the network-topology analysis.

The contrast: the same systems' **internal SPADs** give 0.15 efficiency (IDQ
receiver) — that is the "QKD appliance" reading, not the "deployed network"
reading. We use the latter because it is the fairer comparison against CV's
homodyne receiver and reflects what deployed QKD networks actually field.

---

## 3. Detailed sourcing, parameter by parameter

### DV detector efficiency = 0.65

- Tang et al. 2016 (PRX 6, 011024), deployed three-user MDI network in Hefei
  over 17-30 km installed fibre: SNSPDs at 2.1 K, **system efficiency 64% and
  66%**, DCR 100 Hz. This is the canonical deployed-network paper.
- Intermodal QKD 2025 (arXiv:2602.16680): IDQ **commercial SNSPDs at 80%**
  efficiency in field operation; internal SPADs only 15%.
- Boston field-test QKD (arXiv:1708.00434): WSi SNSPDs **>85%** efficiency.
- Multimode-fibre field trial (Toshiba, arXiv:2410.18646): SNSPD **50%**, DCR <10 Hz.
- **0.65 = deployed-network baseline; sensitivity band 0.65-0.93 (lo: Tang 2016; hi: Swedish 303km field trial, arXiv:2606.06107). NbN 0.30 excluded as low outlier.**

### DV dark count = 100 cps

- Tang 2016 deployed network: **100 Hz**.
- Range across deployed SNSPD experiments: 1 Hz (200 km MDI, cold-filtered) to
  ~1000 Hz (Boston WSi). 100 Hz is the centre.
- Convert to per-gate: at 1 GHz with ~1 ns window, p_dc ~ 1e-7 per gate.
- **100 cps = deployed baseline; sensitivity band 1-1000 cps (lo: Swedish field <=1 cps, arXiv:2606.06107; hi: Boston ~1000 cps, arXiv:1708.00434).**

### DV intrinsic QBER = 0.005

- Intermodal 2025 IDQ SNSPD field run: **QBER < 1%**.
- Boaron 2018 (ULL fibre): 0.5% at short distance.
- Clavis XGR commercial field deployment: <1% average.
- **0.005 (0.5%) = baseline; sensitivity band 0.5-2.1% (lo: Swedish SNSPD link 0.5%; hi: same trial InGaAs 2.1%, arXiv:2606.06107). Swept to ~11% BB84 cutoff in the figure.**
- **Engine note:** TNO's `BB84FullyAsymptoticKeyRateEstimate` derives the channel
  error rate from dark counts and `polarization_drift`, NOT from the detector's
  `error_detector` field (which is inert in the asymptotic rate path). Intrinsic
  QBER is therefore encoded as a drift angle, `polarization_drift = arcsin(sqrt(QBER))`,
  which reproduces the target QBER to <1e-3 across 0-11% (verified against the
  library's `compute_gain_and_error_rate`). An earlier version passed `error_detector`
  and ran DV at ~zero misalignment error; this was corrected.

### CV detector efficiency eta = 0.60

- Lodewyck 2007 (all-fibre): 0.606. Zhang 2020 (record): 0.6134.
- Deployed/commercial-fibre CV (Zhang 2017 50 km): comparable.
- CV homodyne efficiency has not improved in 13 years (coupling + photodiode QE
  limit, not technology gap). Best lab bulk Si-photonics TBHD: 0.72.
- **0.60 = deployed baseline; sensitivity band 0.60-0.72 (lo: Lodewyck 2007 / standard CV 0.6; hi: Si-photonic TBHD 0.72, arXiv:2305.03419). PSA-enhanced 0.965 and integrated 0.37 excluded (different tier).**

### CV electronic noise v_el = 0.10

- Zhang 2020 (deployment-grade homodyne): 0.12-0.27 SNU.
- Best lab balanced homodyne (Chi 2010): ~0.05 SNU (13 dB below shot).
- **0.10 = deployed baseline; sensitivity band 0.05-0.11 (lo: best balanced homodyne ~0.05; hi: deployed field BPD 0.11 SNU, npj QI 2025).**

### CV excess noise xi = 0.015 SNU (Bob-referred)

- Huang 2016 (100 km, controlled excess noise): 0.015 SNU.
- Zhang 2020: 0.0015 (short) to 0.008 (200 km), channel-input referred.
- Commercial-fibre field tests (Zhang 2017): ~0.01-0.03 SNU.
- **0.015 SNU = deployed baseline; sensitivity band 0.01-0.03 (lo: 5G fronthaul field 0.0146, arXiv:2104.04360; hi: commercial-fibre field ~0.03).**
- **NB code conversion (CORRECTED):** the qosst-skr `xi` argument is the
  CHANNEL-INPUT-referred excess noise; the Bob-referred value is converted as
  `xi_input = xi_bob / (T * eta)`. This matches QOSST's own documentation example
  (`GaussianChannel(T, xi_bob / (T * eta))`) and was calibrated against QOSST's
  own asymptotic calculator (qosst-sim) to <15% (residual = ideal-Gaussian vs
  discrete-QAM modulation, expected). An earlier version used
  `2 * xi_bob / (T * eta)` — the extra factor of 2 is NOT in the QOSST reference,
  made the excess noise 2x too high, and collapsed the CV range to ~23 km. With
  the factor of 2 removed, CV reaches ~34 km at baseline (see caveat 7).

### CV reconciliation beta = 0.95

- Deployed CV: 0.95-0.96 standard (Zhang 2020 uses 0.95 short, 0.98 only at
  extreme 200 km low-SNR).
- **0.95 = deployed standard; sensitivity band 0.90-0.96 (lo: Lodewyck/SECOQC field 0.9; hi: rate-adaptive MET-LDPC 0.964, arXiv:1703.04916). Older homodyne 0.80 noted as conservative floor.**

### CV modulation variance Va = optimised (capped at 10 SNU)

Va is **not a fixed system parameter** like the rows above. It is the
modulation strength Alice chooses, and the standard is to report the secret
key rate at its optimal value. The code therefore optimises Va numerically at
each distance (max over Va of the qosst-skr rate) rather than fixing one number.

- Experimental GG02 systems operate at Va ~ 3-6 SNU: optimisation studies find
  Va_opt ~ 2.8 (asymptotic) and ~3.7 at 50 km; deployed/integrated systems
  calibrate to ~5.8-8 SNU; some integrated receivers run lower (~0.5).
- The **unconstrained asymptotic optimum** only exceeds ~10 SNU at sub-km
  distance, where it runs to >100 SNU — not physically achievable on real
  modulators and below the distances where the comparison is contested.
- **Cap Va <= 10 SNU** to reflect deployed hardware. This affects only the
  0-1 km points (lowering them slightly); the DV/CV crossover (~18 km) and all
  points from 5 km out are unchanged, as their optima already sit below 10 SNU.

### Fibre attenuation = 0.20 dB/km (baseline); sweep motivated by Hollow-Core Fibre

**Baseline value.** Standard SMF-28 at 1550 nm; deployed field fibre 0.19-0.21
dB/km. **0.20 = deployed standard.**

**Why a sweep needs a different justification (per supervisor).** Sweeping
attenuation around 0.20 is not physically meaningful: essentially all deployed
silica fibre worldwide sits at ~0.2 dB/km, and silica loss has been effectively
flat for 40 years (0.154 dB/km in 1985 -> 0.1396 in 2024 — a hard material
floor). There is no real deployed variation to explore.

**The sweep is instead motivated by Hollow-Core Fibre (HCF).** HCF guides light
through an air core rather than glass, which removes most Rayleigh scattering and
breaks the silica loss floor (it also propagates near vacuum speed, ~45% faster,
though that affects latency not key rate). This gives a genuine, current reason
for attenuation to span a wide range:

- **Demonstrated:** 0.091 dB/km at 1550 nm — Petrovich et al., "Broadband optical
  fibre with an attenuation lower than 0.1 dB/km," Nature Photonics 19, 1203-1208
  (2025) (Southampton/Microsoft, DNANF design). First fibre of any kind to beat
  the silica floor.
- **Projected:** models suggest ~0.01 dB/km is achievable (same line of work).
- **Conventional silica:** 0.20 dB/km (status quo).

**Sweep 0.01-0.25 dB/km** spans [HCF projected -> HCF demonstrated (0.091) ->
conventional silica (0.20)]; sensitivity band shaded 0.091-0.20 (HCF demonstrated
to silica). HCF has substantial operational challenges (splicing, bending,
connectorisation) and is not yet a deployed QKD medium — so this sweep is a
forward-looking "what if QKD ran on HCF" study, framed as such, not a claim
about present deployment. Because CV is loss-limited, lower attenuation is
expected to benefit CV more than DV.

### Clock rates: DV = 1 GHz, CV = 100 MHz  (section 4a)

Clock rate is a **hardware-speed assumption, separate from the protocol physics.**
It does not affect bits-per-channel-use (Figure 2) at all — it only scales
bits/s (Figure 4). It is documented here for completeness and flagged as a
distinct choice rather than folded into the physics parameters.

Deployed clock rates span a wide range on both sides:

| System | Type | Clock |
|---|---|---|
| Tang 2016 (Hefei) | Deployed DV research network | 75 MHz |
| Boston silicon photonics | DV field test | 625 MHz |
| Clavis XGR | Commercial DV product | 1 GHz |
| Boaron 2018 | Lab DV (ULL fibre) | 2.5 GHz |
| Lodewyck 2007 | CV all-fibre | 0.35 MHz (effective) |
| Zhang 2020 | CV record | 5 MHz |
| QOSST / modern CV | High-rate CV | 100 MHz |
| Newest integrated CV | Chip-based CV | >1 GBaud |

**Baseline choice: DV 1 GHz / CV 100 MHz** — both at the "current
commercial-grade" end of their deployed ranges (DV matches the Clavis XGR
commercial product; CV matches modern high-rate / QOSST-class systems). This
is the fair "what you can buy today" reading.

**Alternative for sensitivity / discussion:** the *deployed research network*
clocks are much lower — DV 75 MHz (Tang 2016) and CV 5 MHz (Zhang 2020). Using
those would widen the DV:CV clock ratio from 10x to ~15x and lower all of
Figure 4 by 1-2 orders of magnitude, without changing Figure 2 at all. Worth
noting in the report as the conservative deployed-field reading.

**Why this is kept separate:** clock speed reflects engineering maturity and
cost, not the security physics of the protocol. Figure 2 (bits/channel use)
deliberately removes it to isolate the protocol comparison; Figure 4 (bits/s)
deliberately includes it to show the deployment reality. The DV clock advantage
in Figure 4 is a statement about current DV hardware being faster, NOT about DV
being a more efficient protocol per channel use.

---

## 4. Honest caveats (for the report's critical-thinking section)

1. **DV detector choice is the key modelling decision.** Deployed DV-QKD — from
   city field networks (Tang 2016) to long-haul field trials (303 km, Swedish
   2026) — fields cooled SNSPDs (eta 0.5-0.93). If instead one assumed the
   low-cost SPADs built into commercial QKD appliances (eta ~0.15), DV would
   look far worse. We use SNSPD because (a) it is what deployed QKD networks
   field at every scale this project studies, and (b) it is the fair match to
   CV's good homodyne receiver. The sweep shows both regimes.

2. **CV detector efficiency has stagnated** at ~0.6 for 13 years — a physics
   limit, not a closing gap.

3. **Excess-noise convention:** xi is Bob-referred; code converts to channel
   input via `xi_input = xi_bob / (T * eta)` (no factor of 2), matching QOSST's
   documentation example and calibrated against QOSST's own calculator. See the
   corrected NB under section 3 (CV excess noise) and caveat 7.

4. **Clock rates are not "fair" across protocols:** DV ~1 GHz, deployed CV
   ~5-100 MHz. This is hardware speed, not protocol. Bits-per-channel-use
   (Fig 2) isolates the protocol; bits/s (Fig 4) shows deployment reality.
   Baseline uses commercial-grade clocks (DV 1 GHz, CV 100 MHz); the deployed
   research-network reading (DV 75 MHz, CV 5 MHz) would lower Fig 4 by 1-2
   orders of magnitude but leave Fig 2 unchanged. See section 4a.

5. **These are baseline values for the sensitivity sweep**, not final claims.
   The sweep (at ~15 km) varies each parameter around these baselines.

6. **DV engine corrections (affects all DV numbers).** Three issues in the TNO
   path were found and fixed: (a) intrinsic QBER was not entering the rate
   (`error_detector` is inert; now encoded via `polarization_drift`, see section 3);
   (b) the reconciliation correction silently failed because `optimize_rate`
   returns a dict (`mu_opt["mu"]`); (c) reconciliation is hardcoded ideal in the
   engine and applied post-hoc. With these fixed, DV carries its genuine 0.5% QBER
   and 0.95 reconciliation cost. All pre-fix DV figures are superseded.

7. **CV excess-noise convention correction (affects all CV numbers).** The
   Bob-to-channel-input conversion previously carried an erroneous factor of 2
   (`2*xi_bob/(T*eta)`), making the excess noise 2x too high and collapsing the
   CV range to ~23 km. The correct conversion is `xi_bob/(T*eta)` — matching
   QOSST's own documented example and calibrated against QOSST's asymptotic
   calculator (qosst-sim). With this fixed, CV reaches ~34 km at baseline and
   the DV/CV crossover (bits/channel use) sits at ~18 km. All pre-fix CV figures
   are superseded. NB: the CV engine assumes ideal Gaussian modulation (qosst-skr),
   so it sits ~10-15% above QOSST's discrete-QAM reference — expected, and the
   correct protocol-level idealisation (mirrors idealised decoy-BB84 for DV).

8. **Corrected headline results (post-fix, asymptotic, bits/channel use):**
   DV/CV crossover ~18 km (CV wins below, DV above); CV reach ~34 km;
   DV reach ~275 km; homodyne and heterodyne track within a few percent.

9. **Single baseline across two studies.** The point-to-point figures evaluate
   this baseline over 0-50 km; the network-topology study applies the same
   baseline to real topologies whose links span tens to thousands of km. The
   parameters are hardware-tier properties and carry across scales; what
   changes with scale is the *outcome* (which links each protocol can serve),
   which is precisely the quantity under study.

---

## 5. Comparison to a peer's DV preset (sanity check)

A peer's independent DV parameter table (different codebase) used:
eta_d = 0.65 (labelled "SNSPD metropolitan" in their naming), d_c = 100 cps,
alpha = 0.20, source error 0.005, f_EC = 1.2 (beta ~ 0.83). These match the DV
detector and channel values derived here from Tang 2016 and the deployed SNSPD
literature, providing independent corroboration of the DV baseline. (Their
beta ~ 0.83 differs from our 0.95, so the corroboration is of the
detector/channel values, not reconciliation.)

---

## 6. Primary sources

### Deployed / field-tested DV
- Tang et al., "MDI-QKD over untrustful metropolitan network," PRX 6, 011024
  (2016), arXiv:1509.08389. **Canonical deployed field network.** eta 0.64-0.66,
  DCR 100 Hz.
- Intermodal QKD, arXiv:2602.16680 (2025). IDQ commercial SNSPD 80%, QBER <1%.
- "Metropolitan QKD with silicon photonics" (Boston), arXiv:1708.00434. WSi
  SNSPD >85%.
- "Metro-scale QKD using multimode fibre" (Toshiba), arXiv:2410.18646. SNSPD
  50%, DCR <10 Hz.
- Clavis XGR commercial deployment, arXiv:2604.16236 (2026). (SPAD appliance
  reading: eta ~0.15, QBER <1%.)
- TFLN on-chip BB84 over field-deployed Korea Telecom fibre, APL Photonics 10,
  031301 (2025). 32 km commercial-node channel; **dark count 120 cps, QBER 0.58%**.
- GHz-rate polarization BB84, arXiv:2602.08908 (2026). 1.5 GHz clock, intrinsic
  QBER ~0.4%.
- Detector-integrated on-chip QKD receiver, npj QI 7 (2021). Waveguide SNSPD
  >90% efficiency, <1 cps dark count, 2.6 GHz clock (bounds the sweep top end).
- Deployed field trial over 303 km, arXiv:2606.06107 (2026). Field SNSPDs
  0.92-0.93, dark count <=1 cps; SNSPD link QBER 0.5%, InGaAs link QBER 2.1%.
  Anchors the high ends of the DV efficiency and QBER sensitivity bands.

### Deployed / field-tested CV
- Zhang et al., "CV-QKD over 50 km commercial fiber," Quantum Sci. Technol. 4,
  035006 (2019), arXiv:1709.04618. Field test on commercial fibre.
- Zhang et al., "Long-Distance CV-QKD over 202.81 km," PRL 125, 010502 (2020),
  arXiv:2001.02555. eta 0.6134, full parameter table.
- Lodewyck et al., "QKD over 25 km all-fibre CV system," PRA 76, 042305 (2007).
  eta 0.606 — the homodyne validation anchor. Defines the noise convention:
  chi_tot = chi_line + chi_het/T referred to channel input (matches qosst-skr).
- Huang et al., "Long-distance CV-QKD by controlling excess noise," Sci. Rep. 6,
  19201 (2016). Excess noise 0.015 at 100 km.
- Telefonica/Huawei/UPM Madrid CV field trial on commercial fibre (2018).

### CV engine / QOSST
- Pietri et al., "QOSST: A Highly-Modular Open Source Platform for Experimental
  CV-QKD," Quantum 8, 1575 (2024), arXiv:2404.18637. The qosst-skr package used
  here. Documents excess-noise convention: skr takes channel-input xi; example
  channel init is GaussianChannel(T, xi_bob/(T*eta)) — the source of the
  corrected conversion (no factor of 2).

### Context (lab records — used only to bound the sweep ranges)
- 5G fronthaul CV-QKD field study, arXiv:2104.04360. 13.2 km urban link, total
  excess noise xi = 0.0146 SNU at 250 MHz — corroborates xi = 0.015 baseline.
- Silicon-photonics TBHD CV-QKD, arXiv:2305.03419. V_A = 5.8, v_el = 0.1,
  beta = 0.98; fibre TBHD efficiency 0.53, on-chip 0.28.
- Shen et al., GC-LDPC reconciliation, Chin. Phys. B 29, 040301 (2020).
  Reconciliation efficiency 95.42% — corroborates beta = 0.95.
- Integrated silicon CV-QKD over 50.4 km, arXiv:2508.08722. Integrated receiver
  efficiency 0.37 (chip-scale, below deployed-fibre 0.60 tier).
- SNSPD 98% efficiency: Reddy et al., Optica 7, 1649 (2020).
- SNSPD 0.01 cps dark count (at 5.6% eta): Shibata et al., Opt. Lett. 39, 5078
  (2014).
- Best lab homodyne (0.72 eta): Wang et al., arXiv:2305.03419 (2023).

### Hollow-Core Fibre (motivates the attenuation sweep, section 3)
- Petrovich et al., "Broadband optical fibre with an attenuation lower than 0.1
  decibel per kilometre," Nature Photonics 19, 1203-1208 (2025), arXiv:2503.21467.
  DNANF hollow-core fibre: 0.091 dB/km at 1550 nm (record, first to beat silica),
  <0.2 dB/km across 66 THz; models project toward 0.01 dB/km. Pilot involved
  >1200 km of installed fibre carrying live traffic.
- Context: silica minimum loss 0.154 dB/km (1985) -> 0.1396 (2024), ~flat for
  40 years — the reason a sweep around 0.2 silica is not meaningful.

### Detection schemes (per supervisor steer)
Both homodyne and heterodyne CV detection are included throughout. A commercial
CV-QKD system tested in the group's lab used a homodyne detector; other systems
use both. qosst-skr provides asymptotic trusted calculators for each
(GaussianTrustedHomodyneAsymptotic / GaussianTrustedHeterodyneAsymptotic), so
both are evaluated on a single consistent framework in the main comparison and
the sensitivity analysis.