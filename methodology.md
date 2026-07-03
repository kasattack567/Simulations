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