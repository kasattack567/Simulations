# Project Roadmap — DV vs CV QKD Comparison

MSc dissertation (UCL). Supervisor: Alejandra Beghelli.
Secondary supervisor: Emilio Hugues-Salas.

**Last updated: August 2026.** Numbers re-verified by direct execution of the
engines. Supersedes the version that described the project as a
"point-to-point metropolitan (0-50 km)" study — the scope has grown.

---

## 1. The project, in one line

Compare **discrete-variable (DV)** and **continuous-variable (CV)** quantum key
distribution head-to-head — point to point, then across **real optical network
topologies at native geographic scale**, then on **cost** — using validated
simulators, and report which protocol performs better, where, why, and at what
trusted-node budget.

### Why it matters
DV and CV are the two main families of QKD. DV (single-photon, decoy BB84) is
mature and tolerates high loss; CV (coherent states, homodyne or heterodyne
detection) uses telecom-style hardware but is more sensitive to noise and
reconciliation. Which to deploy is an open, practical question, and it is
increasingly a *per-link* question rather than a one-off procurement choice:
hybrid encoders and SDN-controlled networks (MadQCI) already switch between DV
and CV on a link-by-link basis. The honest answer is "it depends on the metric
and the topology" — this project makes that precise, and supplies the crossover
and relay budget those hybrid architectures currently assume.

### Protocols and engines
- **DV: decoy-state BB84**, computed with the **TNO** library
  (`tno.quantum.communication.qkd_key_rate`), built from Attema 2021 and
  Ma-Fung-Lo 2007.
- **CV: GG02 Gaussian-modulated coherent states**, homodyne and heterodyne,
  computed with **qosst-skr** (QOSST ecosystem).
- Both run **asymptotic** (infinite-key) so the comparison is symmetric. Note
  this is a deliberate choice, not purely a tooling limit: TNO *does* ship a
  finite-key class, but using it on DV alone would bias the comparison, and
  CV's parameter-estimation penalty is generally harsher — so the asymptotic
  choice mildly favours CV. Stated as a limitation *with its direction*.

### Reporting metrics (deliberately separated)
- **Bits per channel use** — strips out hardware speed, asks which protocol is
  *fundamentally* better.
- **Bits per second** — includes deployed per-protocol clock rates (DV 1 GHz,
  CV 100 MHz), asks which system to *deploy today*.
- **Trusted-node count** — the network-level currency: what each protocol costs
  in relays to light up a given real topology.
- **Break-even cost ratio** — normalised, so the reader supplies their own cost
  assumptions.

See `methodology.md` for why these answer different questions, and for the trap
that a bits/s *criterion* silently inherits the clock choice.

---

## 2. What we have done

### Simulator validation (against published experiments)
- **DV / TNO vs Attema 2021** (its source paper): reproduces Figure 1a — same
  slope, cutoff at the stated ~40 dB. Primary asymptotic DV validation.
- **DV / TNO vs Boaron 2018** (experiment): agreement 0.50x-1.12x across
  252-405 km. Differences traced to three causes (tighter Attema finite-key
  bound, 4-state vs 3-state, fixed vs swapped detectors) — none are bugs.
- **DV / TNO vs MQZL 2005**: constant ~2.67x offset — TNO uses tighter modern
  bounds; documents that no single universal asymptotic decoy formula exists.
- **CV / qosst-skr vs Lodewyck 2007** (homodyne source paper): reproduces
  published mutual information and Holevo bound to ~0.1%, key rate to ~1%, and
  the full Figure 5 distance curve. Strongest CV validation.
- **CV vs Wang et al. 2019** (the excess-noise model's source): reproduces
  Fig. 1 reach (404.6 km vs ~400 km at 1e-9 for xi = 0.01) and the optimal
  modulation variance (Va = 3.712 vs 3.71); tracks Table 2 across 0-100 km.
- **CV vs commercial hardware**: 94 km / 18.8 dB against the LuxQuanta NOVA LQ
  Gen-2 spec of 100 km / 20 dB. Independent, and not a paper.
- **CV vs Kish 2024 / Takeoka-Pirandola bounds**: curves sit correctly relative
  to channel-capacity bounds.

### Key technical settlements
- **CV excess noise now follows Wang Eq. 12**, a two-parameter decomposition
  (`eps_a + eps_l` constant floor, `eps_b` amplified by 1/(eta*T)) rather than a
  single constant xi. This supersedes *both* earlier forms — the erroneous
  factor-of-2 conversion and the corrected single-parameter one. A constant xi
  cannot reproduce CV's hard distance ceiling. **Consequence: CV reach moved
  ~34 -> 94 km and the crossover ~18 -> ~50 km.** See `param.md` caveat 7.
- **DV engine corrections**: intrinsic QBER now enters via `polarization_drift`
  (the `error_detector` field is inert in the asymptotic path); the
  reconciliation correction no longer silently fails on `optimize_rate`'s dict
  return; beta = 0.95 applied post-hoc since the engine is ideal.
- **CV engine choice: qosst-skr, not qosst-sim.** qosst-sim's `nu=1`
  simplification drops detector electronic noise from Eve's information, making
  it optimistic near the cutoff.
- **Homodyne added.** At realistic parameters homodyne ~ heterodyne (within a few
  percent; ~0.9% apart at 50 km, identical reach to 0.1 km). The textbook 3 dB
  homodyne advantage does not survive realistic noise.
- **Survivorship-bias fix**: network rates average over ALL pairs, with
  unreachable pairs counted as zero, not over surviving links only.
- **Rescaling topologies rejected as unsound**: national backbones are used at
  native scale, because rescaling destroys the link-length distribution the
  topology encodes.

### Evaluation parameters
Validation and evaluation parameters are kept separate per supervisor's steer.
Evaluation uses a common standard — **"best demonstrated in deployed QKD systems
/ current commercial-grade hardware"** — so DV and CV compete on a fair footing.
Full sourcing in `param.md`. Headline baseline:
- DV (SNSPD): eta 0.65, dark count 100 cps, QBER 0.5%, beta 0.95, clock 1 GHz
- CV (homodyne/heterodyne): eta 0.60, v_el 0.10, eps_a+eps_l 0.005 SNU,
  eps_b 0.0005 SNU, beta 0.95, clock 100 MHz
- Shared: fibre 0.20 dB/km
DV uses **deployed SNSPD** numbers (Tang 2016 field network), not low-cost
SPADs — putting DV (0.65) and CV (0.60) on comparable detector footing.

### Point-to-point results (locked)
- **Crossover 49.9 km** at beta = 0.95; soft, spanning **27.2-56.8 km** across
  the deployed beta band 0.90-0.96. Quote as a range.
- **CV reach 94.4 km; DV reach 279.1 km.**
- Sensitivity analysis complete on five axes at L = 50 km. Headlines: beta gates
  the CV rate entirely while DV is flat (CV overtakes at beta = 0.9502);
  low-loss/hollow-core fibre benefits CV disproportionately (CV/DV 1.00 at
  0.20 dB/km rising to 4.63 at 0.01); detector and channel noise bind neither
  protocol at metro range. Full write-up in `methodology.md`.

### Network / topology stage (complete)
- **Topology selection**: six topologies chosen from TopologyBench's 105 real
  networks by the dataset's own prescribed PCA + k-means procedure (no
  survivability pre-filter, which would have silently deleted the metro tail):
  **TATANID, SAGO, GERMANY50, CESNET, LAYER42, RNPBRAZIL**. Reproducible with
  fixed seeds. A data-integrity issue in the official `mega_graph_metrics.csv`
  (CANARIE19 duplicated, CANARIE24 missing) is handled with a dedup warning.
- **Direct-link (no relay) coverage at native scale** — the result that motivates
  the whole relay analysis. DV link coverage 0-100%, CV 0-24%; and on the two
  long-span topologies (LAYER42, RNPBRAZIL) **neither protocol carries a single
  edge unaided**. Real backbones are simply not reachable without trusted nodes.
- **Relay / trusted-node budget** with midpoint-chain relays sized to a service
  target. At deployed clocks and 10 Mbps: **CV needs 2.8-3.8x more relays than
  DV** across all six topologies.
- Supporting suite: random-network sweeps, canonical topology shapes
  (ring/star/tree/mesh), user-count and trade-off studies.

### Cost stage (complete)
Follows Karavias et al. (ONDM 2025, the supervisor's own paper): normalise one
component to 1, sweep the ratio, report the fractional difference with the
zero-crossing as break-even. Output is a **break-even cost-ratio surface** over
(service target `c_min`, trusted-node site cost `c_site`), for both deployed and
matched clocks. Absolute unit costing is deliberately avoided — a cryogenic
SNSPD receiver and a room-temperature coherent receiver are not comparable units.

### Deliverables produced
`param.md`, `methodology.md`, `engines_summary.md`, validation scripts/plots,
the point-to-point comparison (Figures 2/3/4), the five sensitivity figures,
the full network suite, the cost model, and the drafted dissertation
Introduction with `refs.bib`.

---

## 3. What is left

### Immediate
1. **Supervisor sign-off on the Wang Eq. 12 model change** — this moved the
   headline numbers substantially, so it should be flagged explicitly rather
   than folded quietly into the results.
2. **Fix the clock default in `net_common.py`.** It currently defaults to
   *matched* clocks (`CV_CLOCK_HZ = 1e9`), so `topo_realscale_relays.py
   --criterion rate` run bare reproduces a 1.1x relay ratio, not the 2.8-3.8x
   headline. Either flip the default to deployed or make the script refuse to
   run without an explicit choice. This is the single most likely source of an
   inconsistent number reaching the report.
3. **Purge superseded figures.** Any plot generated before the Wang model change
   is invalid. Regenerate or delete — do not let a ~34 km CV curve survive into
   the report.

### Writing
- Introduction: **drafted** (Context -> Problem -> Gaps -> Proposal, per the
  supervisor's prescribed structure), with hybrid DV/CV prior work as the framing
  for the gap statement.
- Remaining: Abstract (Context -> Problem -> Solution method -> Results),
  Methods, Results, Discussion, and Conclusions (Past -> Present -> Future, per
  the prescribed structure).
- Report must include an explicit **critical-thinking / limitations** section:
  asymptotic-vs-finite-size (with the direction of its bias), the deployed
  parameter assumptions, the clock-rate framing, the trusted-detector
  assumption, the absence of per-node insertion loss, and the homodyne
  "same hardware, different measurement" caveat.

### Optional if time permits
- **Untrusted-detector CV robustness check.** qosst-skr ships
  `GaussianUntrustedHomodyneAsymptotic`; running it would bound how much of CV's
  standing depends on the trusted-detector assumption. Cheap to do, and it
  pre-empts an obvious viva question.
- **Per-node insertion loss** in the network model.
- **Hybrid allocation model**: assign CV below the crossover and DV above it on
  the same topology, and cost the result against single-technology builds. This
  is the natural next step given the framing and would directly test the
  allocation rule the hybrid literature assumes.

---

## 4. Open questions / watch points

- **Verify the Attema 2021 and Wang 2019 citations** against publisher records
  before submission (both are used as primary sources; both were reconstructed
  rather than read from the original PDF).
- **Small fixes from the audit** (full detail in AUDIT.md): cv_hetro.py
  final-table `{L:6d}` crash (figures save first, so outputs unaffected);
  network.py deprecated `cm.get_cmap`; net_common floor comment 265 -> 252 km;
  topo_loader docstring lists an old selection of six; select_topologies
  docstring says medoids but the code does MAXIMIN sampling (methods chapter
  must describe maximin); cost_common's FD is (C_DV - C_CV)/C_CV with positive
  meaning CV cheaper — match this sign convention in the report.
- **Why qosst-skr heterodyne runs ~25% above the QOSST paper** beyond the `nu=1`
  attribution — the Lodewyck validation anchors qosst-skr independently, so this
  is a nice-to-have rather than a blocker.
- **RESOLVED (audit): LAYER42/RNPBRAZIL zero direct coverage is genuine.**
  LAYER42's six nodes are Seattle, San Francisco, Los Angeles, Chicago, New
  York and Washington DC; its shortest edge (NYC-DC) is 491.8 km, longest
  (SF-DC) 4,898 km — every edge exceeds DV's 279 km unaided reach. RNPBRAZIL's
  shortest edge is 376.5 km. Edge lengths audited against great-circle
  distances (file lengths embed a x1.25-1.5 routing detour, as expected).
- **NEW dataset note (audit): two CESNET edges are shorter than great-circle**
  (edge 9-7: 74.7 km stored vs 157.3 km haversine, ratio 0.475; edge 11-9,
  ratio 0.81) — physically impossible for fibre; an upstream TopologyBench
  quirk affecting both protocols identically. Report beside the CANARIE issue.
- **Relay ratios under `--criterion reach`** are undefined for SAGO and CESNET
  (DV needs zero relays). Decide how to present an infinite ratio honestly;
  reporting the absolute counts alongside is probably the cleanest answer.
