"""
Network foundation for the DV vs CV metropolitan comparison.

Model (BB84-style, matching the peer's bb84_network.py, NOT MDI):
  - N users placed uniformly at random in a square metro area.
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
# POINT-TO-POINT THRESHOLDS (from the corrected P2P run; for annotation only)
# ============================================================
CROSSOVER_KM = 23.0    # CV wins (bits/channel use) below this, DV above
CV_REACH_KM  = 40.0    # CV hits zero around here
DV_REACH_KM  = 275.0   # DV hits zero around here

DETOUR = 1.0           # 1.0 = Euclidean (stated simplification); ~1.5 = metro fibre


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
def relays_needed(L_km, reach_km=CV_REACH_KM):
    """Number of equally-spaced relays to split a length-L link into segments
    each < reach_km. 0 if the link already succeeds."""
    if L_km <= reach_km:
        return 0
    import math
    n_segments = math.ceil(L_km / reach_km)
    return n_segments - 1


def relay_positions(users, i, j, reach_km=CV_REACH_KM):
    """(x,y) positions of the relays that evenly split link i-j into segments
    each < reach_km. Empty if none needed. Placed on the straight line."""
    L = pair_distance_km(users, i, j)
    n_r = relays_needed(L, reach_km)
    if n_r == 0:
        return []
    a, b = users[i], users[j]
    # relays at fractions 1/(n_r+1), 2/(n_r+1), ... along the segment
    return [tuple(a + (b - a) * (k / (n_r + 1))) for k in range(1, n_r + 1)]


def relayed_link_rate(users, i, j, protocol, reach_km=CV_REACH_KM):
    """End-to-end rate of link i-j WITH midpoint-chain relays: the bottleneck
    (minimum) rate over the equal sub-segments. Returns (rate, n_relays).
    If no relays needed, returns the direct rate and 0."""
    L = pair_distance_km(users, i, j)
    n_r = relays_needed(L, reach_km)
    if n_r == 0:
        return link_rate(L, protocol), 0
    seg_len = L / (n_r + 1)
    seg_rate = link_rate(seg_len, protocol)   # all segments equal length
    return seg_rate, n_r   # bottleneck = the (identical) segment rate


def relayed_network(users, protocol, reach_km=CV_REACH_KM):
    """Apply midpoint-chain relays to every failed link. Returns a dict with
    per-pair rescued rates, total relays used, and relay positions."""
    rates, total_relays, relay_pts = {}, 0, []
    for i, j in all_pairs(len(users)):
        r, n_r = relayed_link_rate(users, i, j, protocol, reach_km)
        if r > POSITIVE:
            rates[(i, j)] = r
        if n_r > 0 and r > POSITIVE:
            total_relays += n_r
            relay_pts.extend(relay_positions(users, i, j, reach_km))
    return dict(rates=rates, total_relays=total_relays, relay_pts=relay_pts)