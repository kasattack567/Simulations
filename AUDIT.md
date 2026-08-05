# Full Project Audit — DV vs CV QKD Comparison

**Date:** August 2026. **Method:** every script in the project was executed
against the installed engines (TNO `qkd_key_rate`, `qosst-skr`); the source code
of both engines was read line-by-line; the CV secret-key-rate formulas were
re-implemented independently from the textbook covariance-matrix expressions and
compared numerically; all saved figures and CSVs were regenerated and compared
against the project copies; the topology data files were audited edge-by-edge
against great-circle geometry.

**Bottom line:** the science is sound. The two engines are used correctly, the
post-hoc β correction is mathematically exact, the headline numbers reproduce,
the topology selection and the cost table reproduce *exactly*, and the
LAYER42/RNPBRAZIL zero-coverage result is genuine geography. What the audit
found is a handful of stale docstrings, one crash-after-figures bug, two claims
whose validation conditions need stating precisely, and two upstream dataset
quirks worth a line in the report.

---

## 1. Theory verification — all sound

### 1.1 DV: TNO `BB84FullyAsymptoticKeyRateEstimate` (source read)

The engine computes, per pulse at optimised intensity μ:

    R = P0(μ)·Y0 + P1(μ)·Y1·[1 − h(e1)] − Q(μ)·h(E)

with Y0 = 1 − (1−p_dc)², Y1 = 1 − (1−p_dc)²(1−η_sys), and
e1 = [Y1 − (1−p_dc)·η_sys·cos(2θ)]/(2·Y1). This is the standard fully-asymptotic
decoy bound (with the Lo–Ma–Chen vacuum credit and no sifting factor, i.e.
efficient BB84 in the asymptotic limit).

- **The post-hoc β correction is exact, not approximate.** The engine's EC term
  is literally `gain·h(E)` (Shannon limit, f=1). Subtracting
  (1/β − 1)·gain·h(E) therefore produces exactly the rate with f = 1/β. Verified
  at source level.
- **`error_detector` is confirmed inert** in this rate path — the misalignment
  error enters only through `polarization_drift`. The workaround is necessary,
  not optional.
- **drift = arcsin(√QBER) is exact**, not a fit: at p_dc → 0 the engine's
  e1 = (1 − cos 2θ)/2 = sin²θ identically.
- **Dark-count units confirmed**: `dark_count_rate` is a per-gate probability
  per detector (it enters as (1−p_dc)², two detectors), so 100 cps at a 1 GHz /
  ~1 ns gate → 1e-7 per gate is correct. *Caveat:* the constant
  `DARK_CPS_TO_PERGATE = 1e-9` in `sens_common.py` hard-codes the 1 GHz gate; it
  would not follow a changed DV clock. Fine as used, worth one sentence.
- **New caveat found:** `optimize_rate` bounds μ ∈ [0.1, 0.9] by default and the
  project never overrides it. Measured: μ* runs 0.899 (0 km) → 0.500 (275 km),
  so the bound never binds in the region of interest — but at sub-km distances
  μ* presses against 0.9, so the very-short-range DV rate is marginally
  bound-limited. Negligible; state if pressed.

### 1.2 CV: `qosst-skr` trusted asymptotic calculators (source read + independent check)

The heterodyne and homodyne classes implement the Fossier 2009 / Lodewyck 2007
covariance-matrix formulas verbatim:

- χ_hom = (1+V_el)/η − 1, χ_het = (2 − η + 2V_el)/η, χ_tot = χ_line + χ_det/T
- I_AB: full log₂ ratio for heterodyne; ½·log₂ for homodyne (verified
  algebraically identical to ½log₂(1 + ηTV_a/(1 + V_el + ηTξ)))
- Holevo via the standard symplectic eigenvalues (A, B, C, D expressions match
  the published forms exactly; the heterodyne class subtracts g(0)=0 for a fifth
  unit eigenvalue, harmless).

**Independent verification:** the heterodyne SKR was re-implemented from the
textbook expressions with no reference to the qosst-skr source and compared at
three (Va, T, ξ, η, V_el, β) points spanning the operating range. Agreement to
**1e-12** — exact. The CV engine is therefore validated independently of the
QOSST project itself.

### 1.3 The Wang Eq. 12 excess-noise assembly

`cv_xi_input(T) = (ε_a+ε_l) + ε_b/(η·T)` matches Wang Opt. Express 27, 13372
(2019), Eq. 12, and is passed as qosst-skr's channel-input ξ with V_el carried
separately — no double counting. The double dependence of the CV rate on η and α
(through detection *and* through ξ_r) is real and correctly implemented.

---

## 2. Headline numbers — all reproduce

| Claim | Re-derived | Status |
|---|---|---|
| Crossover 49.9–50.6 km at β=0.95 | 49.9 (bisection) / 50.6 (cv_hetro 0.5 km grid) | ✓ |
| Crossover soft, 27.2–56.8 km over β 0.90–0.96 | 27.2 / 33.9 / 49.9 / 56.8 | ✓ |
| CV reach 94.4 km (hom = het) | 94.4 both | ✓ |
| DV reach 279.1 km | 279.1 | ✓ |
| Six-topology selection reproduces exactly | same six; CANARIE dedup warning fires; PCA var 49.5%+26.3% | ✓ |
| Relay ratio 2.8–3.8× (deployed clocks, 10 Mbps) | 2.8–3.8 across all six | ✓ |
| Relay ratio 2.7–7.0× / inf (matched, reach criterion) | 2.7–7.0; SAGO & CESNET: DV=0 | ✓ (this is the saved bar chart) |
| Cost table_summary.csv | regenerated **identical to 0.000%** on every matched grid point | ✓ |
| Wang Fig. 1 reach "404.6 km at 1e-9" | 404.6 km — condition pinned: const-ξ=0.01 homodyne rate falls below **1e-9 bits/symbol** | ✓ with condition |
| "Optimal Va 3.712 vs Wang's 3.71" | 3.712 — condition pinned: **long-distance limit** of const-ξ=0.01 homodyne (converges from ~200 km; at 50 km under the Wang model it is 4.6) | ✓ with condition |
| LuxQuanta 94 km / 18.8 dB vs 100 km / 20 dB spec | 94.4 km ⇒ 18.9 dB | ✓ |

All saved figures in the project (figure2/3/4, s1–s6, topo maps, topo_direct,
pca_selection, cost_summary, relay bar chart) were regenerated and are
**current** — pixel differences are DPI/matplotlib-version rendering only. No
stale pre-Wang figure survives in the project folder.

---

## 3. Errors found (fix before submission)

### Code

1. **`cv_hetro.py` crashes at the final results table.** `DISTANCES_KM` is a
   float array and the print uses `{L:6d}` → `ValueError` at the first row. All
   three figures save *before* the crash, so outputs are unaffected, but the
   tabular results never print. Fix: `{int(L):6d}` (and gate the `L % 5 == 0`
   test on the half-km grid, which currently also fires only on integers by
   luck of float modulo).
2. **`network.py` uses `cm.get_cmap`,** deprecated since matplotlib 3.7 and
   removed in 3.11 — currently a warning, will become a crash on newer
   environments. Same pattern likely in the relayed variant. One-line fix.

### Documentation / comments that contradict the code or the measurements

3. **`net_common.py` floor comment overstates DV's floored reach.** It says the
   1 kbps floor moves DV "from 279.1 km to about 265 km — roughly 5%". Measured:
   the rate crosses 1e-6 b/pulse at **251.9 km** (span_km bisection agrees:
   252.1), i.e. a ~10% shrink. The qualitative point stands; the number is wrong.
4. **`topo_loader.py` docstring lists the wrong six topologies**
   (DARKSTRAND, USA100, REDIRIS, LAMBDARAIL, NETRAIL, HIBERNIAUK — an earlier
   selection). The `TOPOLOGY_FILES` dict is correct.
5. **`select_topologies.py` docstring contradicts its implementation.** The
   header says "2 MEDOIDS per cluster … typical members chosen over extremes";
   the code (and its inline comment) implements **farthest-point / maximin
   sampling on the PCA plane, seeded at each cluster's |PCA1| extreme**, per the
   paper's Sec. 6.A.4 wording. The methods chapter must describe the maximin
   version. (The selection itself reproduces exactly either way.)
6. **`roadmap.md` "Key technical settlements" still stated the factor-of-2
   conversion as the fix** — superseded twice over; already corrected in the
   updated roadmap.
7. **`network.py` docstring mixes old and new reach figures** ("CV reach
   ~40 km" alongside "~94 km") — cosmetic.
8. **`cost_common.fd` sign convention:** the code defines
   FD = (C_DV − C_CV)/C_CV, positive ⇒ CV cheaper. Earlier project notes state
   the opposite orientation ((C_CV − C_DV)/C_DV). If the report quotes FD,
   match the code's convention or flip it knowingly. (The summary figure only
   uses R*, which is unaffected.)

### Claims to soften / state precisely

9. **"Tracks Wang Table 2 across 0–100 km" is too strong.** Measured against
   Wang's reported input-referred excess noise: 0 km 0.0058 vs 0.0052 ✓;
   50 km 0.0133 vs 0.0158 ✓; 100 km 0.0883 vs 0.0979 ✓; but **20 km 0.0071 vs
   0.0134 and 80 km 0.0382 vs 0.0721 — the model is ~½ the measured value at the
   intermediate points.** Wang's own wording is "agree approximately". Phrase
   as: "matches at 0/50/100 km within ~15% and is conservative by roughly a
   factor of two at the intermediate points, consistent with Wang's own
   'approximate' characterisation of the fit."
10. **Homodyne ≈ heterodyne needs a range qualifier.** Measured hom/het: 0.77 at
    1 km, 0.88 at 5 km, 0.95 at 10 km, ~1.01–1.08 from 20 km to the cutoff.
    So: "within ~5% beyond ~15 km (homodyne marginally ahead); heterodyne up to
    ~25% ahead below 5 km; the hom/het crossover sits near 15–20 km" — not the
    roadmap's "~10 km", and not "within a few percent" globally.
11. **State the validation conditions** for the two Wang checks (row table
    above): 404.6 km is at a 1e-9 b/sym floor; Va→3.712 is the long-distance
    limit under const-ξ. Both true; neither is meaningful without its condition.

---

## 4. Dataset findings (report in the data/methods section)

- **CANARIE duplication** in the official `mega_graph_metrics.csv` (CANARIE19
  duplicated, CANARIE24 missing vs the paper's Table 9): confirmed live — the
  dedup warning fires and 105 → 104 unique rows. Already handled by the script.
- **LAYER42 / RNPBRAZIL zero direct coverage is genuine geography, resolved.**
  LAYER42's six nodes are Seattle, San Francisco, Los Angeles, Chicago, New
  York, Washington DC — a US coast-to-coast network whose *shortest* edge
  (NYC–DC) is 491.8 km and longest (SF–DC) 4,898 km; every edge exceeds DV's
  279 km unaided reach. RNPBRAZIL's shortest edge is 376.5 km. The claim will be
  challenged and now has a one-line answer with the node list.
- **New: two CESNET edges are shorter than the great-circle distance** between
  their endpoints (edge 9–7: 74.7 km stored vs 157.3 km haversine, ratio 0.475;
  edge 11–9: ratio 0.814) — physically impossible for fibre, an upstream
  TopologyBench data quirk (likely a node-coordinate or mapping error). Effect:
  those two edges are slightly flattered, identically for both protocols; no
  ranking changes. Worth one sentence next to the CANARIE note — it strengthens
  the "dataset was audited" narrative.
- **Detour factors:** the files' "Computed Length" embeds routing detours over
  great-circle of ×1.25–1.5 (TATANID and SAGO exactly 1.5 throughout;
  GERMANY50 mostly ×2.3–2.5). Meanwhile the *synthetic* studies use Euclidean
  distances (DETOUR = 1.0, a stated simplification). The real-topology results
  therefore already include realistic fibre routing; the synthetic ones are a
  lower bound on fibre length. Say this once — it pre-empts an examiner
  question about mixing conventions.

---

## 5. Design quirks (not bugs, but know they exist)

- **Va cap differs between the two engine wrappers**: `cv_hetro.py` searches Va
  up to 100 SNU (documented as deliberately unconstrained for the sub-km
  region); `sens_common.py` — which feeds the sensitivity, network, relay and
  cost layers — caps at 10 SNU per param.md. Verified irrelevant beyond ~5 km
  (all optima < 10 there), so no result is affected, but the report should
  quote **one** convention (the 10 SNU cap) and note Figure 2's 0–1 km segment
  is upper-idealised.
- **`topo_direct.py` pins a 1 GHz clock locally** for its bits/s panel but takes
  its viability floor from `net_common.reachable_threshold`, which follows the
  env-configurable clocks. With defaults (matched) the two agree; if
  `QKD_CV_CLOCK_HZ=100e6` were exported, the CV floor would tighten by 10×
  while the displayed clock stayed at 1 GHz. Under your teacher-approved
  matched-clock policy this never triggers — but don't export that env var when
  regenerating this particular figure.
- **`REACHABLE_BPS` floor asymmetry**, as documented: any positive floor costs
  DV proportionally more (gentle tail) than CV (cliff). Measured: DV 279→252 km
  (−10%), CV 94.4→94.3 km (−0.1%). Quote 252 km wherever "usable DV reach"
  appears.
- **Matched-clock cost surface is *not* globally degenerate.** The ~1.06 R* is
  specific to (c_min = 10 Mbps, s = 0). At matched clocks, R* = 0.66 at
  100 Mbps (CV cheaper even at equal system prices — its higher bits/channel-use
  means fewer parallel systems per link) and R* = 1.3–7.0 at 1 Mbps as site cost
  grows (DV's longer reach saves relays). Since the whole report is going
  matched-clock per the supervisor, this *becomes the cost chapter's story*:
  the winner flips with the service target, and the surface shows exactly where.

---

## 6. Everything executed, at a glance

| Script | Ran | Output vs saved |
|---|---|---|
| cv_hetro.py | ✓ (crashes at final print, figures fine) | figures 2/3/4 current |
| s1–s6 + run_all_sensitivity | ✓ | all six current |
| select_topologies.py | ✓ | CSV + PCA plot reproduce exactly |
| topo_loader self-check | ✓ | lengths audited vs haversine |
| topo_maps.py | ✓ | current |
| topo_direct.py | ✓ | current; LAYER42/RNPBRAZIL zeros genuine |
| topo_realscale_relays.py (4 clock/criterion variants) | ✓ | saved chart = reach/matched |
| network.py / _sweep / _users / _topology / _shapes / _tradeoff | ✓ | sensible; one deprecation warning |
| cost_common + topo_cost_summary | ✓ | table_summary.csv identical to 0.000% |

## 7. Recommended actions (30 minutes of fixes)

1. `cv_hetro.py`: `{int(L):6d}`.
2. `network.py`: `plt.get_cmap` (×2 call sites).
3. `net_common.py` floor comment: 265 → 252 km.
4. `topo_loader.py` docstring: correct the six names.
5. `select_topologies.py` docstring: replace the medoid paragraph with the
   maximin description (methods chapter must match).
6. Reconcile the FD sign convention wherever the report mentions it.
7. In the report, use the softened/conditioned phrasings of §3.9–3.11 and the
   dataset sentences of §4.
