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
# Refreshed for the Wang et al. 2019 Eq. 12 CV excess-noise model (see
# sens_common.py). The previous values (23 / 40 / 275) came from the superseded
# xi_bob/(T*eta) parameterisation and must not be reused.
# ============================================================
CROSSOVER_KM = 50.0    # CV wins (bits/channel use) below this, DV above.
                       # SOFT: 27.2 km at beta=0.90, 56.8 km at beta=0.96, so the
                       # deployed beta band alone spans 27-57 km. Quote as a range,
                       # not a point. This value is the beta=0.95 case.
CV_REACH_KM  = 94.4    # CV hits zero here (agrees with the 100 km / 20 dB spec of
                       # the LuxQuanta NOVA LQ Gen-2 commercial CV-QKD system)
DV_REACH_KM  = 279.1   # DV hits zero here

DETOUR = 1.0           # 1.0 = Euclidean (stated simplification); ~1.5 = real fibre


# ============================================================
# RELAY SPACING CRITERION — reach vs service rate
# ============================================================
# Relay spacing can be sized two ways and the choice dominates every downstream
# trusted-node-budget result:
#
#   'reach' — each hop must yield a non-zero key. Maximally optimistic, and
#             DEGENERATE for this comparison: DV's 279 km reach means it needs
#             ZERO relays on four of the six topologies, so the CV/DV relay ratio
#             is undefined.
#   'rate'  — each hop must still meet a service key rate TARGET_BPS. This is what
#             the cost model requires (c_min), it is non-degenerate, and it yields
#             a well-defined CV/DV relay ratio of 2.5-4.0x on all six topologies.
#
# 'rate' is the default. Whichever is chosen must be applied to BOTH protocols.
# Spans are hardcoded so importing this module stays fast; re-derive with
# span_km(..., recompute=True) after any engine change.
SPACING_CRITERION = "rate"
TARGET_BPS = 10e6           # 10 Mbps service target (cost-model c_min)
DV_CLOCK_HZ = 1e9           # DV symbol rate (param.md)
CV_CLOCK_HZ = 100e6         # CV symbol rate (param.md)

# (criterion, protocol) -> max hop length, km. 'rate' entries are at 10 Mbps.
# At 1 Mbps: dv 114.9, cv_het 61.6. At 100 Mbps CV cannot meet the target at any
# distance (cv_het 1.6 km, cv_hom 0.0 km) — CV is not a candidate for 100 Mbps.
_SPANS = {
    ("reach", "dv"): 279.1, ("reach", "cv_het"): 94.4, ("reach", "cv_hom"): 94.4,
    ("rate",  "dv"):  65.0, ("rate",  "cv_het"): 23.7, ("rate",  "cv_hom"): 24.0,
}


def span_km(protocol, criterion=None, target_bps=None, recompute=False):
    """Maximum hop length for `protocol` under the chosen sizing criterion.

    criterion='reach' -> longest hop with any positive key.
    criterion='rate'  -> longest hop still delivering target_bps bits/second,
                         using that protocol's clock rate.
    recompute=True re-derives by bisection instead of using the cached table
    (slow; use after changing engine parameters).
    """
    criterion = SPACING_CRITERION if criterion is None else criterion
    target_bps = TARGET_BPS if target_bps is None else target_bps
    cached = (criterion == "reach") or (criterion == "rate" and target_bps == TARGET_BPS)
    if not recompute and cached:
        return _SPANS[(criterion, protocol)]

    clock = DV_CLOCK_HZ if protocol == "dv" else CV_CLOCK_HZ
    thresh = 1e-10 if criterion == "reach" else target_bps / clock
    lo, hi = 0.5, 400.0
    if link_rate(lo, protocol) < thresh:
        return 0.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if link_rate(mid, protocol) >= thresh:
            lo = mid
        else:
            hi = mid
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


# ============================================================
# NETWORK METRICS
# ============================================================
POSITIVE = 1e-9   # a pair "connects" if its rate exceeds this

def coverage(rates):
    """Fraction of pairs achieving a positive key rate."""
    if len(rates) == 0:
        return 0.0
    return float(np.mean(np.asarray(rates) > POSITIVE))


def aggregate_rate(rates):
    """Sum of positive key rates across all pairs (bits/channel use)."""
    r = np.asarray(rates)
    return float(r[r > POSITIVE].sum())


def rate_stats(rates):
    """Summary of the achievable-rate distribution over connected pairs."""
    r = np.asarray(rates)
    live = r[r > POSITIVE]
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
            "coverage": coverage(rates),
            "aggregate": aggregate_rate(rates),
            **rate_stats(rates),
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
        if r > POSITIVE:
            rates[(i, j)] = r
        if n_r > 0 and r > POSITIVE:
            total_relays += n_r
            relay_pts.extend(relay_positions(users, i, j, reach_km))
    return dict(rates=rates, total_relays=total_relays, relay_pts=relay_pts)