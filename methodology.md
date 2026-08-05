# Methodology Notes

> **Currency note.** All numbers in this document were re-derived by direct
> execution of the engines against the locked baseline in `param.md` under the
> Wang Eq. 12 excess-noise model. The sensitivity operating point is **L = 50 km**
> (it was 15 km, then 20 km, in earlier drafts, tracking the crossover as the CV
> noise model was corrected). Any text quoting 15 km, 20 km, a ~18 km crossover,
> or a ~34 km CV reach is stale and superseded.

---

## The clock-rate question

The two figures answer different questions and must not be mixed. Figure 2 (bits
per channel use, equal clock) asks which protocol is *fundamentally* better and
strips out hardware speed — CV leads below the crossover, DV above it. Figure 4
(bits/s, deployed clocks) asks which system to *deploy today* and includes DV's
faster clock — DV wins almost everywhere.

Equalising the clock in Fig 2 is fair because, unlike detector efficiency or
noise, the clock is not intrinsic to the protocol: both could in principle run at
any rate, and integrated CV is already reaching GHz symbol rates. Removing it
isolates the underlying per-channel-use physics. Together the two figures give
the real conclusion: **CV's deployment disadvantage is contingent (clock speed,
improving), not fundamental (per-channel-use, where it leads at short range).**
Keep each claim matched to its figure.

On the current baseline the crossover sits at **49.9 km** at beta = 0.95, and it
is soft — **27.2 km at beta = 0.90, 56.8 km at beta = 0.96**. Quote it as
"roughly 50 km, spanning 27-57 km across the deployed reconciliation band",
never as a single number. CV reach is **94.4 km**, DV reach **279.1 km**
(zero-rate); under the 1 kbit/s reachability floor at 1 GHz the usable figures
are **CV 94.3 km, DV 251.9 km** — the floor costs DV ~10% and CV almost
nothing, because CV's rate cliffs where DV's tail decays gently. (This corrects
a stale net_common comment claiming ~265 km.)

Homodyne vs heterodyne, measured precisely: heterodyne leads by up to ~25%
below 5 km (hom/het 0.77 at 1 km, 0.88 at 5 km), the two cross near 15-20 km,
and homodyne stays 1-8% ahead from there to the common 94.4 km cutoff. Say
"within ~5% beyond 15 km", not "within a few percent" globally.

### The clock is not confined to Figure 4

This is the trap worth flagging explicitly, because it bit the relay analysis.
The clock choice does not stay politely inside the bits/s figure: **any criterion
expressed in bits/s inherits it.** The trusted-node budget is sized by bisecting
against a bits/s goal, so the threshold each protocol faces is `goal / clock`.

- At **deployed clocks** (DV 1 GHz, CV 100 MHz) with a 10 Mbps service target,
  hop spans are DV 65.0 km vs CV 23.7 km, and CV needs **2.8-3.8x** more relays.
- At **matched clocks** (both 1 GHz) with the same target, hop spans are 65.0 vs
  61.6 km and the ratio collapses to **~1.1x**.

The second result is an artefact, not a finding. With a common clock, a common
bits/s target maps to the same per-channel-use threshold for both protocols, and
10 Mbps happens to sit near where the two rate curves cross — so the criterion
lands both spans a few km apart and never looks at the region where the protocols
actually differ. Cost is a deployment question, so the cost model is pinned at
deployed clocks; the relay-count claim must always be quoted with its clock
condition attached. Note that `net_common.py` currently *defaults* to matched
clocks, so the deployed run needs `QKD_CV_CLOCK_HZ=100e6` set explicitly.

---

## Fairness of the parameter comparison

DV and CV parameters cannot be made equal — they measure different physical
quantities (DV dark count and CV excess noise have no common unit), so forcing
numerical equality would hide the comparison inside a false equivalence. Nor can
each protocol simply be taken at face value, since DV is the more mature
technology and its deployed numbers are naturally higher — that would conflate
protocol merit with engineering maturity.

The resolution is to fix each protocol at the **same tier of its own demonstrated
range**: best-demonstrated deployed / commercial-grade hardware for both, not lab
records for one and field values for the other. Fairness here is
matched-maturity-tier on a common standard, not matched numbers.

This answers the deployment question ("given the best each can field today, which
wins where?"). Figure 2 partially isolates protocol merit by removing the most
maturity-dependent factor (the clock), and the sensitivity sweep quantifies the
maturity gap directly by showing how far the crossover moves as each protocol's
parameters improve.

The asymmetry in the deployed bands is real and should be stated rather than
smoothed over: DV detector efficiency spans 0.65-0.93 while CV spans 0.60-0.72,
because deployed SNSPDs genuinely outperform deployed homodyne receivers. That
asymmetry *is* the maturity gap, and showing it is part of the argument.

---

## Parameter sensitivity sweep — findings at L = 50 km

Each parameter is varied in turn around the locked baseline at the crossover
distance, where both protocols are competitive and the sweeps are therefore
informative. Because the engines are deterministic analytic formulas, there is
no run-to-run variance and the plots carry **no error bars** (unlike a
Monte-Carlo simulator such as NetSquid).

One structural consequence of the Wang Eq. 12 model worth stating in the report:
**eta and alpha now enter twice for CV** — once through detection, once through
`xi_r(T) = (eps_a + eps_l) + eps_b/(eta*T)`. That coupling is physical and
intended, and it is why the CV curves in s1 and s3 are steeper than they were
under the old constant-xi model.

### s1 — Detector efficiency

Both protocols rise with eta, roughly parallel on the log axis, with CV modestly
above DV at equal efficiency (ratio 1.06 at eta = 0.50 rising to 1.13 at
eta = 1.00) and slightly steeper.

The headline is the comparison **at the deployed operating points, not at equal
eta**: DV at its deployed 0.65 gives 2.00e-2 bits/channel use, CV at its deployed
0.60 gives 2.00e-2 — they are level to three significant figures, which is simply
another statement of the fact that 50 km is the crossover.

> *Takeaway:* At equal detector efficiency CV holds a small (~10%) edge and gains
> marginally more per unit of improvement, but the two operate in different
> deployed bands (CV 0.60-0.72, DV 0.65-0.93). DV's higher *achievable*
> efficiency cancels CV's per-unit advantage exactly at the crossover, and
> outweighs it beyond.

**MUST SAY THIS WHEN WRITING THE RESULTS SECTION.** The s1 figure shows the CV
curve sitting above the DV curve at every eta, which looks like it contradicts
the claim that the two are level at the crossover distance. It does not, and the
reason has to be stated explicitly or a reader will think the figure and the
crossover result disagree:

- The crossover is **not** "CV equals DV at every eta". It is "CV at its baseline
  eta = 0.60 equals DV at its baseline eta = 0.65". Two different points on the
  axis, both giving 2.009e-2 bits/channel use.
- The s1 plot puts both protocols on one shared eta axis, so at any single eta
  one protocol is being read away from its baseline. At eta = 0.65 CV is being
  read *better* than its deployed value; at eta = 0.60 DV is being read *worse*
  than its. Either way CV sits above.
- So the curves cross in *parameter space*, not on the plot. The place the two
  meet is where the two shaded band edges line up, not where the lines do.

Also worth a sentence: the CV homodyne and heterodyne curves look coincident in
s1 because they differ by only 0.8-1.3% across the whole sweep, and the y-axis
spans a factor of 100. Same reason the DV/CV gap (6-13%) looks smaller than it
is. If the separation needs to be visible, tighten the y-limits (data spans
~1.5e7 to 3.5e7 in the matched-clock panel), at the cost of the two panels no
longer sharing a y-axis.

### s2 — Reconciliation efficiency **(headline plot)**

DV is essentially flat (1.95e-2 at beta = 0.80 through 2.01e-2 at beta = 1.00, a
3% span). CV climbs steeply across the same range — 6.87e-3 to 3.31e-2, a factor
of 4.8 — and crosses DV at **beta = 0.9502**.

> *Takeaway:* DV is almost insensitive to reconciliation efficiency, because beta
> only scales its small error-correction term. CV depends on beta critically:
> beta gates the entire CV rate. At 50 km, CV beats DV only once beta exceeds
> ~0.95, which is exactly the deployed standard — so CV's competitiveness at
> metro range rests on reconciliation being at, and staying at, the top of its
> demonstrated band. This is a genuine structural asymmetry between the protocol
> families, not a parameter artefact.

This also explains why the crossover is soft: beta moves the CV curve and barely
touches the DV curve, so the deployed beta band alone swings the crossover from
27 to 57 km.

### s3 — Fibre attenuation (hollow-core fibre motivation)

Both fall as attenuation rises, but at very different rates. CV/DV ratio by
attenuation: 4.63 at 0.01 dB/km, 2.19 at 0.05, 1.52 at 0.091 (HCF demonstrated),
1.16 at 0.15, 1.00 at 0.20 (silica, the crossover by construction), 0.84 at 0.25.

> *Takeaway:* Lower-loss fibre benefits CV disproportionately. At conventional
> silica the two are level at 50 km; at the demonstrated hollow-core value of
> 0.091 dB/km CV leads by ~1.5x, and at the projected 0.01 dB/km by ~4.6x. The
> mechanism is the double dependence noted above — reducing alpha raises T, which
> both improves detection *and* suppresses the eps_b/(eta*T) noise term — so CV
> gains twice where DV gains once.

This answers the supervisor's HCF steer directly: yes, HCF would help, and it
would help CV substantially more than DV. Frame it as a forward-looking "what if
QKD ran on HCF" study — HCF is not yet a deployed QKD medium.

### s4 — Detector noise (weak-response result)

DV is flat across the entire deployed dark-count band and well beyond: 1.9990e-2
at 1 cps, 1.9981e-2 at 1000 cps, and still 1.912e-2 at 100,000 cps — a 4% drop
over five orders of magnitude. CV declines gently and monotonically with
electronic noise: 2.194e-2 at v_el = 0, 1.996e-2 at the deployed 0.10, 1.466e-2
at 0.50, and does not reach zero anywhere in the physically plausible range.

> *Takeaway:* Within deployed ranges, neither protocol is limited by detector
> noise at metro distance. DV is effectively immune to dark counts here — the
> signal simply dominates at 50 km — which is corroborated by an independent
> NetSquid simulation. CV degrades only mildly with electronic noise (~1% per
> 0.01 SNU around the baseline).

State this as a genuine result rather than padding: identifying which parameters
*don't* bind is part of a sensitivity analysis. Do not over-read the CV slope —
the y-axis is heavily zoomed, and across the deployed band 0.05-0.11 the change
is about 5%.

### s5 — Channel noise

Twin axes; the two curves have no common unit and crossings between them are
**not** physically meaningful (the caption says so).

DV falls with QBER: 2.28e-2 at 0.1%, 2.00e-2 at the deployed 0.5%, 1.24e-2 at
2.1% (the top of the deployed band), and reaches zero at **9.8%**, close to the
textbook ~11% BB84 bound.

CV falls with eps_b: 2.20e-2 at 1e-4, 2.00e-2 at the deployed 5e-4, 1.35e-2 at
2e-3 (top of the deployed band), and dies at **eps_b = 5.8e-3**, corresponding to
an assembled xi_r of 0.101 SNU at 50 km.

> *Takeaway:* At 50 km both protocols carry substantial margin against channel
> noise — DV sits at 0.5% against a 9.8% cutoff, CV at 5e-4 against a 5.8e-3
> cutoff, roughly an order of magnitude of headroom each. Channel noise is not
> the binding constraint at metro range for either family.

The more useful thing eps_b controls is not the rate at a fixed distance but the
**reach**: eps_b = 5e-4 gives 94.4 km, eps_b = 2e-3 gives 68.8 km. Since reach
sets the trusted-node budget, eps_b is the CV parameter with the largest
downstream effect on the network and cost results, even though its effect on the
rate at 50 km looks unremarkable. Say this explicitly — it links the sensitivity
chapter to the network chapter.

### The overall story (section intro text)

*At the ~50 km crossover, DV and CV are level per channel use, but they arrive
there by different routes and with different fragilities. CV's parity is
conditional: it requires reconciliation at the top of the deployed band (s2) and
degrades faster than DV as loss rises (s3). DV is remarkably robust across every
axis tested — flat in beta, flat in dark count, and tolerant of QBER to nearly
10%. Detector and channel noise bind neither protocol at metro range (s4, s5).
The parameters that matter for the network-level result are therefore beta and
eps_b: the first sets where the crossover falls, the second sets CV's reach and
hence its trusted-node budget.*

### Job of each plot

| Plot | Job | Strength |
|---|---|---|
| s1 | CV leads at equal eta; DV reaches higher eta in deployment | supporting |
| s2 | CV needs beta >= 0.95 to compete; DV is indifferent to beta | **headline** |
| s3 | HCF helps CV far more than DV | **headline** (ties to supervisor's steer) |
| s4 | Detector noise does not bind at metro range | weak-response result |
| s5 | Channel noise does not bind; eps_b matters via reach, not rate | supporting |

Two strong findings (s2, s3), two "this doesn't matter much" results (s4, s5)
that are still worth stating, and s1 as a solid supporting comparison. Not every
parameter has to be dramatic; showing which ones do not matter is part of the
analysis.

---

## Choices to defend in the viva

1. **Asymptotic on both sides.** A deliberate symmetry choice, not purely a
   tooling limit: TNO does ship `BB84FiniteKeyRateEstimate`, but qosst-skr ships
   asymptotic classes only. Applying finite-key to DV alone would bias the
   comparison, and CV's parameter-estimation penalty is generally the harsher of
   the two — so the asymptotic choice **mildly favours CV**. State the direction
   of the bias, not just its existence.

2. **DV reconciliation applied post-hoc.** The TNO engine hardcodes error
   correction at the Shannon limit and exposes no beta argument, so beta = 0.95
   is applied by scaling the EC cost. This is the standard f(E) scaling and is
   exact. One residual: mu remains optimised at beta = 1, making DV very slightly
   **pessimistic** — a conservative direction for a comparative study.

3. **Trusted-detector model for CV.** Detector loss and electronic noise are
   treated as calibrated and not attributed to Eve. This is the standard
   realistic model and what deployed systems assume, but qosst-skr also ships
   `GaussianUntrustedHomodyneAsymptotic`, so the paranoid variant is available as
   a robustness check and is currently unused. Flag as a limitation.

4. **No per-node insertion loss** in the network model. Real deployments add mux,
   demux and switch losses at every node; the MadQCI trial reports paths made
   unfeasible by exactly this. Relay budgets are therefore slightly optimistic
   for both protocols — but for both, so the ratio is more robust than the
   absolute counts.

5. **Ideal Gaussian modulation for CV.** qosst-skr assumes it, so the CV rate
   sits ~10-15% above a discrete-QAM reference. This is the correct
   protocol-level idealisation and mirrors the idealised decoy-BB84 treatment on
   the DV side — the two idealisations are matched, which is what fairness
   requires here.