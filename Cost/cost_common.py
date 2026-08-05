"""
cost_common.py — DV vs CV network cost model (Karavias-style), shared machinery.

Follows the framework of Karavias, Hugues-Salas, Beghelli, Lord, Payne,
"Comparative Cost Analysis of Prepare and Measure and Entanglement-Based QKD
Networks" (ONDM 2025): provision the minimum hardware that meets a key-rate
demand c_min, normalise one component's cost to 1, sweep the other's relative
cost, and report the fractional difference in total network cost; the zero
crossing is the break-even cost ratio.

MAPPING TO DV vs CV. In the paper the asymmetry is where sources and detectors
sit; here it is reach -> relay chains and per-hop rate -> parallel systems. The
cost unit is therefore a SYSTEM: one transmitter + one receiver operating one
hop. A chain of n hops contains exactly n Tx and n Rx whatever the protocol, so
a two-ratio (Tx, Rx) model collapses to a single system-level ratio with no
loss of generality for chain provisioning.

PROVISIONING RULE (per edge, per protocol). An edge of native length L split
into n equal hops of length l = L/n needs, per hop, m = ceil(c_min / r_bps(l))
parallel systems, where r_bps is the per-system secret key rate at that
protocol's clock. The edge then uses n*m systems and n-1 trusted-node SITES
(one site hosts all m parallel systems' Rx+Tx; site cost covers secure housing
/ operations, NOT hardware, which is already counted in n*m — no double count).
n is chosen to MINIMISE edge cost

    cost(n) = n*m(n) * c_sys  +  (n-1) * c_site        [per protocol]

so each protocol picks its own cost-optimal relay spacing, and that optimum
shifts with the site-to-system cost ratio q = c_site / c_sys. Hops must clear
the suite's reachability floor (net_common.REACHABLE_BPS, 1 kbit/s) to count
at all — same validity floor as every other script.

Dividing by c_sys, provisioning depends only on q, and

    C_P = c_sys(P) * [ N_sys(P; q_P)  +  q_P * N_TN(P; q_P) ].

NORMALISATION. c_sys(CV) = 1. Sweep R = c_sys(DV)/c_sys(CV) (DV, with its
cryocooled SNSPD receiver, is the plausibly expensive system — consistent with
param.md's SNSPD-tier DV baseline; the cost ratio must describe the same
hardware tier as the performance parameters). Site cost s = c_site in CV-system
units. Then q_CV = s, q_DV = s/R, and

    C_CV = N_sys_cv + s*N_TN_cv                (constant in R)
    C_DV(R) = R*N_sys_dv(s/R) + s*N_TN_dv(s/R) (concave, increasing in R)

so FD(R) = (C_DV - C_CV)/C_CV crosses zero exactly once; the crossing R* is the
break-even. With s = 0 provisioning is R-independent and R* = N_sys_cv/N_sys_dv
analytically (the FD curves are straight lines); s > 0 makes DV re-optimise its
relay spacing as R changes, which is where the curvature comes from.

CLOCKS. Cost is a deployment question (a Figure-4 question, not Figure-2), so
the model is pinned at DEPLOYED clocks: DV 1 GHz, CV 100 MHz. This also avoids
the matched-clock degeneracy documented in net_common (at matched clocks a
common bits/s target maps to the same bits/channel-use threshold for both
protocols, near where the rate curves cross). Override per run:
    QKD_COST_CV_CLOCK_HZ=1e9 python topo_cost.py     # matched-clock variant

CONTINUUM INTUITION (useful for the write-up). Ignoring ceilings and site
costs, an edge needs ~ L * c_min / max_l[ l * r_bps(l) ] systems, so total
hardware scales with (total fibre length x c_min) and the protocol comparison
is governed by the ratio of the protocols' DISTANCE-RATE PRODUCTS l*r(l);
topology enters through quantisation (short edges, ceil effects) and through
the trusted-node term once s > 0. topo_cost.py plots l*r(l) as a supporting
figure.

SCOPE / ASSUMPTIONS (state in the report): every physical edge is provisioned
to carry c_min (per-edge service level, not routed end-to-end demands); relays
sit mid-span per edge and are not shared between edges; fibre is assumed
available with unlimited wavelengths (hardware + sites only — no fibre lease,
no quantum-channel co-existence penalty); rates are the suite's asymptotic
engines with the 1 kbit/s validity floor.
"""
import math
import os
import sys
import numpy as np


# ============================================================
# Locate the project modules (net_common, sens_common, topo_loader).
# Works when cost/ sits next to the flat project folder, or next to
# Networks/ + Parameter_sensitivity/. Override with QKD_NET_DIR /
# QKD_ENGINE_DIR if the layout differs.
# ============================================================
def _add_candidates():
    here = os.path.dirname(os.path.abspath(__file__))
    cands = [os.environ.get("QKD_NET_DIR"), os.environ.get("QKD_ENGINE_DIR"),
             here,
             os.path.join(here, ".."),
             os.path.join(here, "..", "Networks"),
             os.path.join(here, "..", "Parameter_sensitivity")]
    for c in cands:
        if c and os.path.isdir(c):
            a = os.path.abspath(c)
            if a not in sys.path:
                sys.path.insert(0, a)
    # net_common has its own locator for sens_common; point it at whichever
    # candidate actually holds sens_common.py.
    if "QKD_ENGINE_DIR" not in os.environ:
        for c in cands:
            if c and os.path.isfile(os.path.join(c, "sens_common.py")):
                os.environ["QKD_ENGINE_DIR"] = os.path.abspath(c)
                break


_add_candidates()

from net_common import REACHABLE_BPS                      # noqa: E402  1 kbit/s validity floor
from sens_common import dv_rate, cv_rate                  # noqa: E402  locked-baseline engines
from topo_loader import (TOPOLOGY_FILES, load_topology,   # noqa: E402
                         edge_lengths)


# ============================================================
# CLOCKS — pinned DEPLOYED (see module docstring). Env-overridable.
# ============================================================
DV_CLOCK_HZ = float(os.environ.get("QKD_COST_DV_CLOCK_HZ", 1e9))
CV_CLOCK_HZ = float(os.environ.get("QKD_COST_CV_CLOCK_HZ", 100e6))

PROTO_LABEL = {"dv": "DV — decoy BB84",
               "cv_het": "CV — GG02 heterodyne",
               "cv_hom": "CV — GG02 homodyne"}


def clock_hz(protocol):
    return DV_CLOCK_HZ if protocol == "dv" else CV_CLOCK_HZ


def set_clocks(dv_hz=None, cv_hz=None):
    """Override the clocks before build_tables(). Callers must read the clocks
    as cost_common.DV_CLOCK_HZ / .CV_CLOCK_HZ (module attributes) rather than
    importing them by value, or labels will report the pre-override figures."""
    global DV_CLOCK_HZ, CV_CLOCK_HZ
    if dv_hz:
        DV_CLOCK_HZ = float(dv_hz)
    if cv_hz:
        CV_CLOCK_HZ = float(cv_hz)
    return DV_CLOCK_HZ, CV_CLOCK_HZ


# ============================================================
# RATE TABLES — one engine sweep per protocol, then interpolation.
# The engines are deterministic, so a dense table + log-linear interpolation
# reproduces them to well under the ceil() granularity the cost model sees.
# Hops shorter than the first grid point are credited the first grid point's
# rate (conservative: short hops are never credited more than the 0.5 km rate).
# ============================================================
L_MIN, L_MAX, L_STEP = 0.5, 300.0, 0.5


class RateTable:
    def __init__(self, protocol, l_min=L_MIN, l_max=L_MAX, step=L_STEP,
                 verbose=True, clock_override_hz=None):
        self.protocol = protocol
        self.clock_hz = float(clock_override_hz or clock_hz(protocol))
        self.L = np.arange(l_min, l_max + step / 2, step)
        fn = {"dv": lambda l: dv_rate(l),
              "cv_het": lambda l: cv_rate(l, detection="heterodyne"),
              "cv_hom": lambda l: cv_rate(l, detection="homodyne")}[protocol]
        if verbose:
            print(f"  building {protocol} rate table "
                  f"({len(self.L)} points, clock {self.clock_hz:.0e} Hz)...",
                  flush=True)
        self.r_bps = np.array([fn(float(l)) for l in self.L]) * self.clock_hz
        pos = self.r_bps > 0.0
        if not pos.any():
            raise RuntimeError(f"{protocol}: engine returned no positive rate")
        self.L_pos = self.L[pos]
        self.log_r = np.log10(self.r_bps[pos])
        self.l_reach = float(self.L_pos.max())      # zero-rate table limit
        self.rmax_bps = float(self.r_bps.max())

    def batch(self, l_km):
        """Vectorised per-system rate (bits/s) at hop length(s) l_km."""
        l = np.asarray(l_km, dtype=float)
        out = np.zeros_like(l)
        inr = l <= self.l_reach
        if inr.any():
            lc = np.clip(l[inr], self.L_pos[0], self.l_reach)
            out[inr] = 10.0 ** np.interp(lc, self.L_pos, self.log_r)
        return out

    def __call__(self, l_km):
        return float(self.batch(np.array([l_km]))[0])

    def l_at_rate(self, r_bps_target):
        """Largest hop length whose per-system rate >= target (0 if none)."""
        ok = self.r_bps >= r_bps_target
        return float(self.L[ok].max()) if ok.any() else 0.0

    def span_floor(self, floor_bps=REACHABLE_BPS):
        """Longest valid hop (rate >= the suite's reachability floor)."""
        return self.l_at_rate(floor_bps)

    def distance_rate_product(self):
        """(L, L*r_bps(L)) arrays and the (l*, value*) maximiser — the
        continuum hardware figure of merit (bit*km/s per system)."""
        prod = self.L_pos * (10.0 ** self.log_r)
        i = int(np.argmax(prod))
        return self.L_pos, prod, float(self.L_pos[i]), float(prod[i])


def build_tables(cv_detection="heterodyne", verbose=True,
                 dv_hz=None, cv_hz=None, dv_table=None):
    """Rate tables for DV and the chosen CV detection. ~30 s (DV dominates).

    dv_hz / cv_hz override the module clocks for this table set only, so
    several clock regimes can coexist in one process. Pass an already-built
    dv_table to reuse it (the DV clock is the same in both regimes, and the
    DV engine is ~25x slower than the CV one)."""
    cv_key = "cv_hom" if cv_detection == "homodyne" else "cv_het"
    dv = dv_table or RateTable("dv", verbose=verbose, clock_override_hz=dv_hz)
    cv = RateTable(cv_key, verbose=verbose, clock_override_hz=cv_hz)
    return {"dv": dv, cv_key: cv}, cv_key


# ============================================================
# PER-EDGE PROVISIONING — joint optimisation over the hop count n.
# ============================================================
def provision_edge(L_km, rt, cmin_bps, q, floor_bps=REACHABLE_BPS):
    """Cost-optimal midpoint-chain provisioning of one edge.

    Minimises cost(n) = n*ceil(c_min/r(L/n)) + (n-1)*q over the hop count n,
    in units of this protocol's system cost (q = c_site/c_sys). Returns a dict
    (n_hops, m_per_hop, systems, relays, hop_km, cost) or None if no hop count
    yields hops clearing the reachability floor (cannot happen for L > 0 while
    the table's shortest point clears the floor).

    Search bound: m(n) is non-increasing in n (shorter hops -> higher rate) and
    bottoms out at m_floor = ceil(c_min/r_max); once m has reached m_floor,
    cost(n) = n*m_floor + (n-1)*q increases strictly with n, so no n beyond the
    first m_floor-achieving hop count needs checking.
    """
    if L_km <= 0:
        return dict(n_hops=0, m_per_hop=0, systems=0, relays=0,
                    hop_km=0.0, cost=0.0)
    l_valid = rt.span_floor(floor_bps)
    if l_valid <= 0:
        return None
    m_floor = max(1, math.ceil(cmin_bps / rt.rmax_bps))
    l_mfloor = rt.l_at_rate(max(cmin_bps / m_floor, floor_bps))
    n_cap = int(math.ceil(L_km / max(l_mfloor, 1e-9))) + 1

    n = np.arange(1, n_cap + 1, dtype=float)
    l = L_km / n
    r = rt.batch(l)
    valid = r >= floor_bps
    if not valid.any():
        return None
    m = np.where(valid, np.ceil(cmin_bps / np.where(valid, r, 1.0)), np.inf)
    cost = n * m + (n - 1.0) * q
    i = int(np.argmin(np.where(valid, cost, np.inf)))
    return dict(n_hops=int(n[i]), m_per_hop=int(m[i]),
                systems=int(n[i] * m[i]), relays=int(n[i]) - 1,
                hop_km=float(l[i]), cost=float(cost[i]))


def provision_network(topo, rt, cmin_bps, q, floor_bps=REACHABLE_BPS):
    """Provision every edge of a topology. Totals in this protocol's
    system-cost units: cost = systems + q*relays."""
    systems = relays = 0
    cost = 0.0
    per_edge = []
    for L in edge_lengths(topo, "real"):
        e = provision_edge(float(L), rt, cmin_bps, q, floor_bps)
        if e is None:
            raise RuntimeError(
                f"{topo['name']}: edge of {L:.0f} km infeasible for "
                f"{rt.protocol} even fully subdivided")
        systems += e["systems"]
        relays += e["relays"]
        cost += e["cost"]
        per_edge.append(e)
    return dict(systems=systems, relays=relays, cost=cost, per_edge=per_edge)


# ============================================================
# NETWORK COSTS in CV-system units, and the break-even solver.
#   c_sys(CV) = 1,  c_sys(DV) = R,  c_site = s.
# ============================================================
class CostModel:
    """Caches provisioning per (topology, protocol, c_min, q) so sweeps and
    bisections stay fast. Topologies are loaded once."""

    def __init__(self, tables, cv_key, floor_bps=REACHABLE_BPS):
        self.rt = tables
        self.cv_key = cv_key
        self.floor = floor_bps
        self.topos = {n: load_topology(n) for n in TOPOLOGY_FILES}
        self._cache = {}

    def provision(self, name, proto, cmin_bps, q):
        key = (name, proto, float(cmin_bps), round(float(q), 6))
        if key not in self._cache:
            p = provision_network(self.topos[name], self.rt[proto],
                                  cmin_bps, q, self.floor)
            self._cache[key] = (p["systems"], p["relays"], p["cost"])
        return self._cache[key]

    # ---- total network costs (CV-system units) ----
    def cost_cv(self, name, cmin_bps, s):
        _, _, c = self.provision(name, self.cv_key, cmin_bps, q=s)
        return c                                   # = N_sys + s*N_TN

    def cost_dv(self, name, cmin_bps, R, s):
        R = max(float(R), 1e-9)
        _, _, c = self.provision(name, "dv", cmin_bps, q=s / R)
        return R * c                               # = R*N_sys + s*N_TN

    def fd(self, name, cmin_bps, R, s):
        """Fractional difference (C_DV - C_CV)/C_CV. Positive => the CV
        network is the cheaper one at this cost ratio."""
        ccv = self.cost_cv(name, cmin_bps, s)
        return (self.cost_dv(name, cmin_bps, R, s) - ccv) / ccv

    def break_even(self, name, cmin_bps, s, lo=0.02, hi=10.0, hi_max=500.0,
                   tol=1e-3):
        """R* with C_DV(R*) = C_CV. C_DV is strictly increasing in R (envelope:
        dC/dR = N_sys > 0), so the crossing is unique; bracket-expand then
        bisect. Returns np.inf if no crossing below hi_max."""
        ccv = self.cost_cv(name, cmin_bps, s)
        g = lambda R: self.cost_dv(name, cmin_bps, R, s) - ccv
        if g(lo) >= 0:
            return lo
        while g(hi) < 0:
            hi *= 2.0
            if hi > hi_max:
                return float("inf")
        a, b = lo, hi
        while b - a > tol * max(1.0, a):
            mid = 0.5 * (a + b)
            if g(mid) < 0:
                a = mid
            else:
                b = mid
        return 0.5 * (a + b)