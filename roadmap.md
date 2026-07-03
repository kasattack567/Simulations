# Project Roadmap — DV vs CV QKD Comparison

MSc dissertation (UCL Quantum Technologies). Supervisor: Alejandra Beghelli.

---

## 1. The project, in one line

Compare **discrete-variable (DV)** and **continuous-variable (CV)** quantum key
distribution head-to-head for **point-to-point metropolitan links (0-50 km)**,
using validated simulators, and report which protocol performs better, where,
and why.

### Why it matters
DV and CV are the two main families of QKD. DV (single-photon, decoy BB84) is
mature and tolerates high loss; CV (coherent/squeezed states, homodyne or
heterodyne detection) uses cheaper telecom-style hardware but historically
struggles with noise and distance. Which to deploy in a metro network is an
open, practical question — and the honest answer is "it depends on the metric,"
which this project makes precise.

### Protocols and engines
- **DV: decoy-state BB84**, computed with the **TNO** library
  (`tno.quantum.communication.qkd_key_rate`), built from Attema 2021 and
  Ma-Fung-Lo 2007.
- **CV: GG02 Gaussian-modulated coherent states**, homodyne and heterodyne,
  computed with **qosst-skr** (the QKD secret-key-rate package of the QOSST
  ecosystem).
- Both run **asymptotic** (infinite-key) so the comparison is fair; CV has no
  mature finite-size proof in these tools, and mixing asymptotic-CV with
  finite-DV would unfairly penalise DV. The asymptotic choice is discussed as
  a limitation in the report.

### Two reporting metrics (deliberately)
- **Bits per channel use** — strips out hardware speed, asks which protocol is
  *fundamentally* better. (Figure 2)
- **Bits per second** — includes realistic per-protocol clock rates, asks which
  system to *deploy today*. (Figure 4)
See `methodology.md` for why these answer different questions.

---

## 2. What we have done

### Simulator validation (against published experiments)
Each engine has been validated against the papers it is built from and against
independent experiments:

- **DV / TNO vs Attema 2021** (its source paper): reproduces their Figure 1a —
  same slope, cutoff at the stated ~40 dB. Primary asymptotic DV validation.
- **DV / TNO vs Boaron 2018** (experiment): agreement 0.50x-1.12x across
  252-405 km. Differences traced to three causes (tighter Attema finite-key
  bound, 4-state vs 3-state, fixed vs swapped detectors) — none are bugs. Tests
  the finite-key path, so secondary to Attema.
- **DV / TNO vs MQZL 2005**: constant ~2.67x offset — TNO uses tighter modern
  bounds; documents that no single universal asymptotic decoy formula exists.
- **CV / qosst-skr vs Lodewyck 2007** (homodyne source paper): reproduces
  published mutual information and Holevo bound to ~0.1%, key rate to ~1%, and
  the full Figure 5 distance curve. Strongest CV validation.
- **CV vs QOSST 2024** (heterodyne): qosst-sim matched the paper's 25 km rate to
  2% (the paper was made with qosst-sim, so partly circular).
- **CV vs Kish 2024 / Takeoka (TGW/PLOB bounds)**: curves sit correctly relative
  to the channel-capacity bounds.

### Key technical settlements
- **Factor-of-2 excess-noise correction** for CV: `xi_input = 2*xi_bob/(T*eta)`,
  fixing an earlier ~3x inflation.
- **CV engine choice: qosst-skr, not qosst-sim.** qosst-sim has a `nu=1`
  simplification that drops detector electronic noise from Eve's information,
  making it slightly optimistic near the CV cutoff. qosst-skr is the faithful
  implementation (validated to Lodewyck 1%) and lets homodyne and heterodyne
  sit in one consistent framework.
- **Homodyne added.** qosst-skr provides an asymptotic trusted-homodyne
  calculator. Finding: at realistic parameters homodyne ~ heterodyne (within
  ~6% across most of the metro range; heterodyne better at short range, homodyne
  near the cutoff, crossover ~10 km) — the textbook 3 dB homodyne advantage does
  not survive realistic noise. Confirms rather than overturns the heterodyne
  results.

### Evaluation parameters (separated from validation parameters)
Per supervisor's steer, validation and evaluation parameters are now separate.
Evaluation uses a common standard — **"best demonstrated in metropolitan
deployments / commercial-grade hardware"** — so DV and CV compete on a fair
footing rather than on whatever each validation paper used. Full sourcing in
`param.md`. Headline baseline:
- DV (SNSPD): eta 0.65, dark count 100 cps, QBER 0.5%, clock 1 GHz
- CV (homodyne): eta 0.60, v_el 0.10, excess noise 0.015 SNU, beta 0.95, clock 100 MHz
- Shared fibre 0.20 dB/km
Key correction this stage: DV uses **deployed SNSPD** numbers (Tang 2016 metro
network), not low-cost SPADs — putting DV (0.65) and CV (0.60) on comparable
detector footing.

### Deliverables produced
`param.md` (parameters + sourcing + clock section), `methodology.md` (two-figures
reasoning), validation scripts/plots (Attema, Boaron, MQZL, Lodewyck, Takeoka,
Kish-format), and the main comparison code producing Figures 2/3/4.

---

## 3. What we are aiming to do next

### Immediate (before final point-to-point results)
1. **Supervisor sign-off on evaluation parameters** (`param.md`) — the
   SNSPD-vs-SPAD choice for DV is the main judgement call.
2. **Re-run Figure 2 on the locked deployed parameters** and confirm where the
   DV/CV crossover lands. CV got worse at this stage (v_el and excess noise up),
   so the earlier "~10 km crossover / CV wins short range" claim must be
   re-verified before it goes in the report.
3. **Decide homodyne placement** — on the main plots as a third CV curve, or as
   a separate validation result (awaiting supervisor view).

### Parameter sensitivity analysis (next major piece)
At a fixed distance (~15 km, near the crossover), vary each parameter in turn to
see which improvements matter most to each protocol. Built on the locked
evaluation baseline as a parameterised template. Functionally-analogous
parameters grouped (pending supervisor approval of the looser pairings):
1. Detector efficiency (clean pair)
2. Reconciliation efficiency (clean pair)
3. Fibre attenuation (identical, shared)
4. Detector noise — DV dark count vs CV electronic noise (looser; twin axes)
5. Channel noise — DV QBER vs CV excess noise (looser; twin axes)

### After point-to-point is finalised
- **Network topology stage** — extend from single links to network structures
  (ring / star / mesh / trusted-node). Not yet started.

### Reporting
- Draft to supervisor in the **first days of August** (her leave: 8-12 Aug,
  18 Aug-1 Sep).
- Report must include explicit **critical-thinking / limitations** discussion:
  asymptotic-vs-finite-size choice, the deployed-parameter assumptions, the
  clock-rate framing, and the homodyne "same hardware, different measurement"
  caveat.

---

## 4. Open questions / watch points

- **Why qosst-skr heterodyne runs ~25% above the QOSST paper** beyond the `nu=1`
  attribution — worth confirming for bulletproof numbers, though the Lodewyck
  validation already anchors qosst-skr independently.
- **DV dark-count unit** (per-gate 1e-7 vs per-window) — sanity-check via the DV
  cutoff distance in the extended-range figure.
- **Crossover location on new parameters** — the single most important thing to
  re-verify before building the narrative.
