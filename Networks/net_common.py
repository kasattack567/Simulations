"""
Network foundation for the DV vs CV comparison on real optical backbones.

Model (BB84-style, matching the peer's bb84_network.py, NOT MDI):
  - N users placed uniformly at random in a square area (synthetic studies only;
    the real-topology scripts use native link lengths from TopologyBench).
  - Every user pair is a candidate DIRECT link (all N(N-1)/2 pairs).
  - Each link's key rate is the point-to-point rate at the pair's fibre distance,
    computed by the SAME corrected engines used in the point-to-point and
    sensitivity stages (dv_rate / cv_rate).
  - Distances are straight-line (Euclidean). Real fibre runs longer (detour
    factor ~1.5 in the literature); we use Euclidean as a stated lower-bound
    simplification, appropriate for an MSc-scale study. Set DETOUR>1 to explore.

Trusted-node relays (added in a later script) split long links into shorter
hops; an end-to-end path rate is then the MIN of its hop rates (bottleneck).

Engines are imported from the sensitivity module so there is ONE source of truth
for the locked parameters and the corrected excess-noise / QBER conventions.
Adjust ENGINE_DIR below to point at the folder containing sens_common.py.
"""
import os
import sys
import numpy as np

# --- locate the engine module (sens_common.py) ---------------------------------
# Point this at the folder that contains sens_common.py. Default assumes a sibling
# "Parameter_sensitivity" folder next to this Networks folder.
ENGINE_DIR = os.environ.get(
    "QKD_ENGINE_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "Parameter_sensitivity"),
)
if os.path.isdir(ENGINE_DIR):
    sys.path.insert(0, os.path.abspath(ENGINE_DIR))

try:
    from sens_common import dv_rate, cv_rate, ALPHA_DB_KM  # corrected engines
except ImportError as e:
    raise ImportError(
        "Could not import sens_common. Set the QKD_ENGINE_DIR environment "
        "variable (or edit ENGINE_DIR in net_common.py) to the folder that "
        f"contains sens_common.py. Original error: {e}"
    )


# ============================================================
# POINT-TO-POINT THRESHOLDS
# Refreshed for the Wang et al. 2019 Eq. 12 CV excess-noise model AND the
# Clason-anchored DV baseline with a 2.0 dB receiver loss budget (see
# sens_common.py). Superseded value sets (23/40/275 and 50/94.4/279.1) must not
# be reused.
# ============================================================
CROSSOVER_KM = 59.0    # CV wins (bits/channel use) below this, DV above.
                       # SOFT: moves with CV beta and with the DV receiver loss
                       # budget (26 km at 0 dB RX loss, 71 km at 3 dB). Quote as
                       # a range, not a point. This is the 2.0 dB / beta=0.95 case.
# NB these are ZERO-RATE distances: the point at which the key rate reaches
# zero. They are NOT the usable reach under the 1 kbit/s reachability floor —
# for that use span_km(protocol, 'reach'), which is shorter. Kept for reference
# and for the plotting/annotation call sites that want the asymptotic limit.
CV_REACH_KM  = 94.4    # CV hits zero here (agrees with the 100 km / 20 dB spec of
                       # the LuxQuanta NOVA LQ Gen-2 commercial CV-QKD system)
DV_REACH_KM  = 303.4   # DV hits zero here

DETOUR = 1.0           # 1.0 = Euclidean (stated simplification); ~1.5 = real fibre


# ============================================================
# RELAY SPACING CRITERION — floor vs service rate
# ============================================================
# Relay spacing can be sized two ways and the choice dominates every downstream
# trusted-node-budget result. Both are derived by bisection against a bits/s
# goal, so both are CLOCK-DEPENDENT (threshold = goal / clock):
#
#   'reach' — each hop must clear the reachability floor, REACHABLE_BPS
#             (1 kbit/s). "Relay only where the link would otherwise fail."
#             This is the sizing the relay-count comparison wants, because it
#             lets each protocol's reach show through.
#   'rate'  — each hop must meet a service key rate, TARGET_BPS (10 Mbps), as
#             the cost model requires (c_min). Note what this does to a DV/CV
#             comparison: with matched clocks both protocols face the same
#             per-channel-use threshold, and 10 Mbps sits near the rate at which
#             the two curves cross, so the spans land a few km apart and the
#             relay counts come out nearly equal — not because the protocols are
#             equally capable at distance, but because the criterion never looks
#             at the region where they differ. Use it for costing, not for reach.
#
# Whichever is chosen must be applied to BOTH protocols.
SPACING_CRITERION = "reach"   # floor-sized; see above
TARGET_BPS = 10e6           # 10 Mbps service target (cost-model c_min)


# ============================================================
# CLOCK RATES — bits/channel use -> bits/s. SINGLE SOURCE OF TRUTH.
# ============================================================
# MATCHED CLOCK: DV and CV are both run at 1 GHz. This is the like-for-like
# protocol-physics view — any DV/CV separation in a bits/s figure is then
# intrinsic to the protocols, not to the electronics.
#
# It is NOT the deployed view. Fielded CV symbol rates are ~100 MHz (DAC/ADC, DSP
# chain and shot-noise-limited receiver bandwidth), so a matched-clock figure
# credits CV with a decade of headroom it does not currently have. Say which view
# a figure uses; CLOCK_LABEL below is provided for exactly that.
#
# Every downstream script imports these rather than defining its own, so the two
# views cannot drift apart again. Switch back per run without editing code:
#     QKD_CV_CLOCK_HZ=100e6 python topo_direct.py
DV_CLOCK_HZ = float(os.environ.get("QKD_DV_CLOCK_HZ", 1e9))
CV_CLOCK_HZ = float(os.environ.get("QKD_CV_CLOCK_HZ", 1e9))   # matched (was 100e6)
CV_CLOCK_DEPLOYED_HZ = 100e6   # fielded CV symbol rate, kept for reference

CLOCK_HZ = {"dv": DV_CLOCK_HZ, "cv_het": CV_CLOCK_HZ, "cv_hom": CV_CLOCK_HZ}


def clock_hz(protocol):
    """Repetition rate (Hz) for `protocol`; use this instead of a local table."""
    try:
        return CLOCK_HZ[protocol]
    except KeyError:
        raise ValueError(f"unknown protocol {protocol!r}")


def to_bps(rate, protocol):
    """bits/channel use -> bits/s at `protocol`'s clock.

    Every script that reports a rate goes through this, so no module multiplies
    by a clock of its own. Accepts a scalar or an array.
    """
    clock = clock_hz(protocol)
    if np.isscalar(rate):
        return float(rate) * clock
    return np.asarray(rate, dtype=float) * clock


def _fmt_hz(f):
    return f"{f/1e9:g} GHz" if f >= 1e9 else f"{f/1e6:g} MHz"


CLOCK_MATCHED = (DV_CLOCK_HZ == CV_CLOCK_HZ)
CLOCK_LABEL = (f"DV {_fmt_hz(DV_CLOCK_HZ)}, CV {_fmt_hz(CV_CLOCK_HZ)}"
               + (" — matched" if CLOCK_MATCHED else ""))


# HOP SPANS. Both criteria are now derived by bisection against a rate floor,
# because a span must be consistent with the reachability test applied to the
# resulting hop — otherwise relays get placed at a spacing whose segments then
# fail is_reachable, and the script reports full coverage while leaving dead
# edges in the graph. That was the case while 'reach' returned DV_REACH_KM /
# CV_REACH_KM (279.1 / 94.4 km): those are ZERO-RATE distances, and a hop of
# that length carries no usable key.
#
#   'reach' -> longest hop still clearing REACHABLE_BPS (the reachability floor,
#              1 kbit/s). "Relay only where the link would otherwise fail."
#   'rate'  -> longest hop still clearing TARGET_BPS (a service level, 10 Mbps).
#              For costing a network at a usable rate.
#
# So the two differ only in the target; 'reach' is the floor-sized case. Spans
# are memoised per (criterion, protocol, target, clock). Import stays fast; the
# first call costs one bisection per protocol.
#
# CONSEQUENCE: relay counts move with REACHABLE_BPS as well as with the clock.
# Any previously quoted CV/DV relay ratio must be re-derived, not reused.
_SPAN_CACHE = {}


def span_km(protocol, criterion=None, target_bps=None, recompute=False):
    """Maximum hop length for `protocol` under the chosen sizing criterion.

    criterion='reach' -> longest hop still clearing the reachability floor
                         (REACHABLE_BPS). Use this for "relay only where the
                         link would otherwise fail".
    criterion='rate'  -> longest hop still delivering target_bps bits/second at
                         that protocol's clock (a service level, for costing).
    Both are clock-dependent, since the threshold is bits/s divided by clock.
    recompute=True forces re-derivation by bisection, ignoring the memo.
    """
    criterion = SPACING_CRITERION if criterion is None else criterion
    target_bps = TARGET_BPS if target_bps is None else target_bps
    # 'reach' is sized at the reachability floor, 'rate' at the service target.
    goal_bps = REACHABLE_BPS if criterion == "reach" else target_bps
    thresh = goal_bps / clock_hz(protocol)

    key = (criterion, protocol, goal_bps, clock_hz(protocol))
    if not recompute and key in _SPAN_CACHE:
        return _SPAN_CACHE[key]

    lo, hi = 0.5, 400.0
    if link_rate(lo, protocol) < thresh:
        _SPAN_CACHE[key] = 0.0
        return 0.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if link_rate(mid, protocol) >= thresh:
            lo = mid
        else:
            hi = mid
    _SPAN_CACHE[key] = lo
    return lo


# ============================================================
# GEOMETRY
# ============================================================
def place_users(n, area_km=30.0, seed=None):
    """N users placed uniformly at random in an area_km x area_km square.
    Returns an (n, 2) array of (x, y) coordinates in km."""
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, area_km, size=(n, 2))


def pair_distance_km(users, i, j):
    """Straight-line fibre distance between users i and j (km), x DETOUR."""
    return DETOUR * float(np.hypot(*(users[i] - users[j])))


def all_pairs(n):
    """Yield all unordered index pairs (i, j) with i < j."""
    for i in range(n):
        for j in range(i + 1, n):
            yield i, j


def all_pair_distances(users):
    """Return a flat array of every pair's fibre distance (km)."""
    n = len(users)
    return np.array([pair_distance_km(users, i, j) for i, j in all_pairs(n)])


# ============================================================
# PER-PAIR KEY RATES (direct links)
# ============================================================
# Cache point-to-point rates by rounded distance to avoid recomputing the CV Va
# optimisation for every pair (many pairs share very similar distances).
_RATE_CACHE = {}

def _cached(fn_key, L_km, fn):
    key = (fn_key, round(L_km, 2))
    if key not in _RATE_CACHE:
        _RATE_CACHE[key] = fn(L_km)
    return _RATE_CACHE[key]


def link_rate(L_km, protocol):
    """Point-to-point key rate (bits/channel use) at distance L_km.
    protocol in {'dv', 'cv_het', 'cv_hom'}."""
    if protocol == "dv":
        return _cached("dv", L_km, lambda L: dv_rate(L))
    if protocol == "cv_het":
        return _cached("cv_het", L_km, lambda L: cv_rate(L, detection="heterodyne"))
    if protocol == "cv_hom":
        return _cached("cv_hom", L_km, lambda L: cv_rate(L, detection="homodyne"))
    raise ValueError(f"unknown protocol {protocol!r}")


def pair_rates(users, protocol):
    """Direct-link key rate for every user pair, given protocol.
    Returns (distances_km, rates) as flat arrays aligned by pair."""
    dists = all_pair_distances(users)
    rates = np.array([link_rate(d, protocol) for d in dists])
    return dists, rates


## ============================================================
# REACHABILITY — the single definition of "this link carries a key"
# ============================================================
# A link counts as REACHABLE if SKR >= REACHABLE_BPS. One definition, used
# everywhere, so coverage figures across the suite mean the same thing.
#
# WHY 1 kbps, NOT ZERO. The engines are asymptotic and return positive rates
# down to the rate-distance limit, including a few bits/s that no finite-size
# analysis would certify. A non-zero floor is the honest reading of "works".
# 1 kbps sits below every deployed system in the literature:
#   - CN-QCN backbone: 9.75-360 kbps across 64 links over 10 weeks (npj QI 11, 2025)
#   - Tokyo metro: 2.8-141 kbps; Shanghai CV metro: 0.25-10 kbps (ACM CSUR 53, 2020)
#   - Application demand: 7.4 kbps for quantum-safe IPsec at 46 km (arXiv:2405.04415);
#     DCI trial 2.392 kbps (arXiv:2410.10245). At 1 kbps ~4 fresh AES-256 keys/s.
#   - Finite-size: no key from fewer than ~10^5-10^6 signals (Scarani, arXiv:1010.0521);
#     Wang 2019 shows finite-size reach ~200 km vs asymptotic ~500 km.
#
# ASYMMETRY. A positive floor costs DV more than CV (Wang Eq. 12's 1/(eta*T)
# makes CV's rate fall vertically at its limit while DV's tail decays gently),
# but the effect is small: DV reach moves ~279.1 km -> ~265 km. Not a large
# correction. Sensitivity runs:
#     QKD_REACHABLE_BPS=256    # one AES-256 key/s
#     QKD_REACHABLE_BPS=1e3    # default
#     QKD_REACHABLE_BPS=1e4    # carrier-grade

REACHABLE_BPS = float(os.environ.get("QKD_REACHABLE_BPS", 256))


def reachable_threshold(protocol):
    """The REACHABLE_BPS floor expressed in bits per channel use, for the
    protocol's clock. Compare raw engine output against this."""
    return REACHABLE_BPS / clock_hz(protocol)


def is_reachable(rate, protocol, units="channel"):
    """True if `rate` clears the reachability floor.

    units='channel' (default) for raw engine output in bits per channel use;
    units='second' if the rate has already been converted with to_bps.
    Accepts a scalar or an array; returns a bool or a boolean array.
    """
    thresh = REACHABLE_BPS if units == "second" else reachable_threshold(protocol)
    if np.isscalar(rate):
        return float(rate) >= thresh
    return np.asarray(rate, dtype=float) >= thresh


# Retained so older call sites keep working; prefer is_reachable(). Note this is
# a bits/channel-use quantity and assumes the DV clock.
POSITIVE = REACHABLE_BPS / DV_CLOCK_HZ


def coverage(rates, protocol="dv"):
    """Fraction of pairs clearing the reachability floor. `rates` in bits per
    channel use. Pass ALL pairs, not just the connected ones."""
    if len(rates) == 0:
        return 0.0
    return float(np.mean(is_reachable(rates, protocol)))


def aggregate_rate(rates, protocol="dv"):
    """Sum of reachable key rates across all pairs (bits/channel use)."""
    r = np.asarray(rates)
    return float(r[is_reachable(r, protocol)].sum())


def rate_stats(rates, protocol="dv"):
    """Summary of the achievable-rate distribution over reachable pairs."""
    r = np.asarray(rates)
    live = r[is_reachable(r, protocol)]
    if len(live) == 0:
        return dict(n_connected=0, mean=0.0, median=0.0, min=0.0, max=0.0)
    return dict(
        n_connected=int(len(live)),
        mean=float(live.mean()),
        median=float(np.median(live)),
        min=float(live.min()),
        max=float(live.max()),
    )


def summarise(users, protocols=("dv", "cv_het", "cv_hom")):
    """Compute coverage + rate stats for each protocol on one user layout."""
    n = len(users)
    n_pairs = n * (n - 1) // 2
    out = {"n_users": n, "n_pairs": n_pairs}
    for p in protocols:
        dists, rates = pair_rates(users, p)
        out[p] = {
            "coverage": coverage(rates, p),
            "aggregate": aggregate_rate(rates, p),
            **rate_stats(rates, p),
        }
    return out


# ============================================================
# TRUSTED-NODE RELAYS — midpoint-chain rescue of failed links
# ============================================================
def relays_needed(L_km, reach_km):
    """Number of equally-spaced relays splitting a length-L link into segments
    each <= reach_km. 0 if the link already qualifies.

    reach_km is REQUIRED and should come from span_km(protocol, criterion): it is
    protocol-dependent, so defaulting it to CV's value (as this function used to)
    silently applied CV's spacing to DV as well.
    """
    if reach_km <= 0 or L_km <= reach_km:
        return 0
    import math
    return math.ceil(L_km / reach_km) - 1


def relay_positions(users, i, j, reach_km):
    """(x,y) positions of the relays that evenly split link i-j into segments
    each < reach_km. Empty if none needed. Placed on the straight line."""
    L = pair_distance_km(users, i, j)
    n_r = relays_needed(L, reach_km)
    if n_r == 0:
        return []
    a, b = users[i], users[j]
    # relays at fractions 1/(n_r+1), 2/(n_r+1), ... along the segment
    return [tuple(a + (b - a) * (k / (n_r + 1))) for k in range(1, n_r + 1)]


def relayed_link_rate(users, i, j, protocol, reach_km=None):
    """End-to-end rate of link i-j WITH midpoint-chain relays: the bottleneck
    (minimum) rate over the equal sub-segments. Returns (rate, n_relays).
    If no relays needed, returns the direct rate and 0."""
    reach_km = span_km(protocol) if reach_km is None else reach_km
    L = pair_distance_km(users, i, j)
    n_r = relays_needed(L, reach_km)
    if n_r == 0:
        return link_rate(L, protocol), 0
    seg_len = L / (n_r + 1)
    seg_rate = link_rate(seg_len, protocol)   # all segments equal length
    return seg_rate, n_r   # bottleneck = the (identical) segment rate


def relayed_network(users, protocol, reach_km=None):
    """Apply midpoint-chain relays to every failed link. Returns a dict with
    per-pair rescued rates, total relays used, and relay positions."""
    reach_km = span_km(protocol) if reach_km is None else reach_km
    rates, total_relays, relay_pts = {}, 0, []
    for i, j in all_pairs(len(users)):
        r, n_r = relayed_link_rate(users, i, j, protocol, reach_km)
        if is_reachable(r, protocol):
            rates[(i, j)] = r
        if n_r > 0 and is_reachable(r, protocol):
            total_relays += n_r
            relay_pts.extend(relay_positions(users, i, j, reach_km))
    return dict(rates=rates, total_relays=total_relays, relay_pts=relay_pts)