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