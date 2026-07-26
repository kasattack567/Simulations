# Methodology Notes

## The clock-rate question

The two figures answer different questions and must not be mixed. Figure 2 (bits
per channel use, equal clock) asks which protocol is *fundamentally* better and
strips out hardware speed — CV leads at short metro range, DV past the crossover.
Figure 4 (bits/s, realistic clocks) asks which system to *deploy today* and
includes DV's faster clock — DV wins almost everywhere. Equalising the clock in
Fig 2 is fair because, unlike detector efficiency or noise, the clock is not
intrinsic to the protocol (both could run at any rate; integrated CV is reaching
GHz) — so removing it isolates the underlying per-channel-use physics. Together
they give the real conclusion: CV's deployment disadvantage is contingent (clock
speed, improving), not fundamental (per-channel-use, where it leads at short
range). Keep each claim matched to its figure. NB: re-run Fig 2 on the new
deployed parameters before committing to "CV wins under ~10 km" — CV got worse
(v_el 0.01->0.10, xi 0.0072->0.015), so the crossover may have moved.

## Fairness of the parameter comparison

DV and CV parameters cannot be made equal — they measure different physical
quantities (e.g. DV dark count vs CV excess noise have no common unit), so
forcing numerical equality would hide the comparison inside a false equivalence.
Nor can each protocol simply be taken at face value, since DV is a more mature
technology and its deployed numbers are naturally higher — that would conflate
protocol merit with engineering maturity. The resolution is to fix each protocol
at the same tier of its own demonstrated range — best-demonstrated metropolitan
deployment for both, not lab records for one and field values for the other.
Fairness here is matched-maturity-tier on a common standard, not matched numbers.
This answers the deployment question ("given the best each can field today, which
wins where?"); Fig 2 partially isolates protocol merit by removing the most
maturity-dependent factor (clock), and the sensitivity sweep quantifies the
maturity gap directly by showing how the crossover moves as each protocol's
parameters improve.

## Parameter senstivity sweep findings

These are all correct now, and each one *does* say something worth extracting. Let me go through them one at a time — what the plot shows, and the one-sentence takeaway you'd write for each. I'll be honest where a plot says less than you'd hope.

First, one thing that changed and you should note: these are now at **L = 20 km**, not 15. That's fine (20 km is closer to your ~18 km crossover, so it's a sensible operating point to probe), just make sure your text says 20 km consistently.

**s1 — Detector efficiency.** Both protocols rise with η, roughly parallel on the log axis. The key extract: **CV sits above DV across the entire efficiency range at 20 km** (both CV curves are ~2× above DV throughout), and **CV is steeper** — it gains more per unit of efficiency. Takeaway: *"At 20 km, CV out-rates DV at any given detector efficiency, and CV benefits more from efficiency improvements — but note the two operate in different deployed η bands (CV ~0.60, DV ~0.65–0.93), so at their real operating points DV's higher achievable efficiency partly compensates."* The shaded bands make that point visually — CV lives in the pink band, DV can reach further right into the blue.

**s2 — Reconciliation efficiency.** This is your strongest single plot. DV is dead flat; both CV curves climb steeply and only overtake DV around β ≈ 0.93. Takeaway: *"DV is essentially insensitive to reconciliation efficiency (it only scales DV's small error-correction term), whereas CV depends on β critically — β gates the entire CV rate. CV only beats DV at 20 km once β exceeds ~0.93, which is why deployed CV needs excellent (0.95+) reconciliation to be competitive."* That's a genuine, quotable finding about a real asymmetry between the protocols.

**s3 — Fibre attenuation (HCF).** Both fall with attenuation; CV is well above DV across the whole range and CV is steeper. The extract that matters: **as you move left toward HCF values (0.091 → 0.01), CV's rate climbs much faster than DV's** — the gap widens dramatically at low loss. Takeaway: *"Lower-loss fibre disproportionately benefits CV: at conventional silica (0.2) CV leads DV modestly, but toward hollow-core values (0.091 and below) CV's advantage grows sharply, because CV is more loss-limited and so gains more when loss falls."* This directly ties your sweep to the HCF motivation Alejandra gave — it answers "would HCF help, and who does it help more?" Answer: yes, and CV more than DV.

**s4 — Detector noise.** DV flat across dark count; CV declines gently with electronic noise. The honest extract: this is the **least informative** of the five, because at 20 km neither protocol is very sensitive over its *deployed* band. Takeaway: *"Within deployed ranges, neither protocol is strongly limited by detector noise at metro distance — DV is essentially immune to dark counts here (corroborated by an independent NetSquid simulation), and CV degrades only mildly with electronic noise."* That's still worth saying — "this parameter doesn't matter much at metro range" is a legitimate, useful result — but don't over-claim significance. Note the y-axis is very zoomed (7×10⁻² to 1.3×10⁻¹), so the CV decline looks steeper than it is; over the *deployed* band the change is small.

**s5 — Channel noise.** Both fall as noise rises; DV (QBER) survives to ~9%, CV (ξ) dies at ~0.055 SNU. The extract: **within their deployed bands both are healthy, and both have comfortable margin before their cutoffs** — DV's cutoff (~11% QBER) and CV's (~0.05 SNU) are both far outside the deployed shaded regions. Takeaway: *"At 20 km both protocols operate with substantial margin against channel noise; CV's excess-noise cutoff (~0.05 SNU) and DV's QBER cutoff (~9–11%) both sit well beyond deployed values, so channel noise is not the binding constraint at metro range."* The crossings here are meaningless (separate axes — your caption already says so).

**The overall story across all five** (this is what to write in the section intro): *At the ~20 km metro operating point, CV out-rates DV per channel use, but its advantage is conditional — it depends critically on good reconciliation (s2) and low loss (s3), the two parameters CV is most sensitive to, while DV is remarkably robust to all of them. Detector and channel noise (s4, s5) are not binding for either protocol at metro range. This explains why CV wins in the short-metro regime only when its supporting conditions (β ≥ 0.95, low-loss fibre) are met.*

That gives each plot a clear job:
- s1: CV leads at equal η, but DV reaches higher η in deployment.
- s2: CV needs good β to win; DV doesn't care about β. **(headline)**
- s3: HCF helps CV more than DV. **(ties to supervisor's HCF steer)**
- s4: detector noise irrelevant at metro range.
- s5: channel noise not binding; both have margin.

Two are strong findings (s2, s3), two are "this doesn't matter much" results (s4, s5) which are still worth stating, and s1 is a solid supporting comparison. That's a legitimate, honest sensitivity chapter — not every parameter has to be dramatic; showing which ones *don't* matter is part of the analysis.


# The Equal Trusted-Node-Budget Frontier — model, usage, and thesis positioning

## What this analysis is

Given a budget of K trusted relay nodes, which protocol (DV decoy-BB84 or
CV GG02 heterodyne) converts them into more network key rate — and at what
trust-exposure cost per delivered key? Each protocol receives its OWN greedy
allocation of the same budget, so the frontier compares each protocol's best
use of identical resources. This is the precise novelty claim verified by the
prior-art sweep: no prior work compares DV and CV **on a common trusted-node-
budget axis** (do NOT claim the broader "no one has compared DV and CV on
networks" — that is falsified by Mariani et al. 2025).

## Model definition (each bullet is a defended methodological choice)

- **Relays live on fibre edges** (midpoint-chain model, consistent with
  `topo_realscale_relays.py`): r relays on an edge of length L subdivide it
  evenly into r+1 hops of L/(r+1); the edge's effective rate is the
  single-hop rate (bottleneck of identical hops).
- **Greedy allocation, coverage-first lexicographic** objective
  (coverage, then survivorship-corrected average pair rate) — the network
  chapter's coverage-first framing applied to resource allocation.
- **Dead-edge revival lookahead**: at backbone scale no single relay can
  revive a 300 km edge for CV, so candidates include "spend the exact batch
  of r relays that brings edge e under the protocol's reach", scored per
  relay spent. During a batch, intermediate budgets hold the pre-batch state
  (a partially subdivided dead edge delivers nothing). This prevents the
  myopic failure mode of one-at-a-time greedy.
- **Routing**: widest-path (max-min bottleneck) over effective edge rates;
  ties broken by fewer trusted intermediaries. Average rate counts ALL pairs
  with unreachable pairs as zero (survivorship correction).
- **Trust exposure counts ALL trusted intermediaries** on the delivered
  path: intermediate topology nodes (which relay keys in the clear under
  the trusted-node model) plus inserted budget relays. Three metrics:
  - `frac_exposed`: fraction of covered pairs transiting >= 1 intermediary
  - `mean_interm`: mean intermediaries per covered pair
  - `interm_per_bit`: rate-weighted intermediaries per delivered key bit,
    sum(rate_p * m_p) / sum(rate_p) — the headline security-cost metric.
- **Units**: allocations are computed once (within-protocol greedy choices
  are invariant to constant rate scaling); bits/channel-use and bits/s
  figures are produced from the same run via per-protocol clocks
  (DV 1 GHz, CV 100 MHz).
- **Engine caching**: TNO / qosst-skr are called only to build a
  distance->rate lookup table (0.5 km grid for DV to 300 km, 0.25 km grid
  for CV to 45 km; log-linear interpolation). The greedy search never
  re-invokes the engines.

## How to run

Place `topo_relay_frontier.py` in `Networks/` alongside `topo_loader.py`
(and with `net_common.py` importable). Then:

    # primary: metro-rescaled frontier at a chosen span
    python Networks/topo_relay_frontier.py --scale metro --span 100 --kmax 20

    # spans worth running: 60 (CV mostly healthy), 100 (contested),
    # 150 (CV coverage-limited) — the frontier story differs across them
    python Networks/topo_relay_frontier.py --scale metro --span 150 --kmax 20

    # extension: real backbone scale (revival regime)
    python Networks/topo_relay_frontier.py --scale real --kmax 30

    # subset for quick iteration
    python Networks/topo_relay_frontier.py --scale metro --span 100 \
        --topologies HIBERNIAUK NETRAIL --kmax 10

Outputs per run: `relay_frontier_<tag>_rate_channel.png`,
`_rate_second.png`, `_exposure.png`, `_elasticity.png`, plus summary tables
and the decay-slope printout.

Expected runtime: table build is the engine-bound part (~600 DV + ~180 CV
calls, once per run). The greedy search itself is seconds for the five
smaller topologies; USA100 (100 nodes) is the slow one — expect minutes,
scale `--kmax` down first if iterating.

## The elasticity mechanism (the claim that makes this a law, not a plot)

A protocol's marginal gain per relay is governed by the steepness of its
rate-distance decay: halving a hop recovers more rate for the
faster-decaying protocol. The script prints both fitted decay slopes
(log10 rate per km, DV over 5-150 km, CV over 5-35 km) and their ratio;
the elasticity figure shows marginal gain per relay for both protocols.
The testable prediction: CV's marginal gain per relay exceeds DV's roughly
in proportion to the slope ratio, wherever both protocols are operating
within reach. At backbone scale the picture inverts on COVERAGE: DV's ~7x
reach advantage means each DV relay revives far more of the network per
relay than CV's (revival cost per edge scales with ceil(L/reach) - 1).

## Thesis positioning — the three nearest prior works and your delta

1. **Mariani et al., arXiv:2504.02372 (2025)** — DV/CV/hybrid on trusted-
   node complex networks. Their axis is node DENSITY and percolation; no
   relay-count budget, no marginal-value-per-relay analysis, no quantified
   trust metric (one qualitative sentence). Cite as nearest neighbour;
   your delta = relay-budget axis + elasticity + trust-exposure metrics.
2. **Karavias et al., ONDM 2025** — ILP cost comparison of PM vs EB QKD on
   switched metro networks. Protocol-ARCHITECTURE comparison on a MONEY
   axis; not DV-vs-CV, not trusted-node-denominated. Your delta =
   re-denominating the budget in trusted nodes and switching the compared
   variable to protocol family.
3. **Selentis-Boulntadakis et al., ONDM 2025 / Makris et al.
   arXiv:2310.17262** — relayed vs switched architectures, single-protocol
   (DV), using rate-decay-function reasoning. Closest to "rate vs hops";
   your delta = protocol family as the variable, plus exposure accounting.

Terminology: "relay elasticity" and "equal trusted-node-budget frontier"
appear unused in the literature — coin them, with a footnote distinguishing
from adjacent existing terms (percolation threshold, path coverage,
cost-per-key, bits-per-relay-use, and the qualitative marginal-trusted-node
observation in Amer, Krawec & Wang, arXiv:2005.12404).

## Caveats for the Limitations section

- Greedy + revival lookahead is a heuristic, not an ILP optimum; the
  frontier is therefore a LOWER bound on each protocol's best use of K
  (state this; it biases neither protocol systematically).
- Even subdivision of edges is optimal for identical hops under bottleneck
  routing but ignores real site availability along fibre routes.
- Intermediate topology nodes counted as trusted intermediaries assumes
  the standard trusted-node network model (keys transit PoPs in the clear);
  MDI/TF or zero-trust relay architectures would change the exposure
  accounting — cite as future work.
- Asymptotic rates, Euclidean-derived edge lengths at metro rescale, fixed
  clock assumptions — inherited from the main study's stated limitations.
- Rate lookup interpolation error is negligible relative to parameter
  uncertainty, but state grid resolutions for reproducibility.