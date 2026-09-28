"""
Network map — DV vs CV side by side on the same user layout.

Places N users uniformly at random in a square area, forms all direct
user pairs (BB84-style, no relays), and draws TWO maps sharing the identical
layout: one for DV, one for CV. On each map, a link is drawn only if that
protocol yields a key over that pair's distance, and coloured by rate. The CV
map visibly loses its longer links while DV keeps them (zero-rate reach: DV
~303 km vs CV ~94 km) — the core network finding, shown as a map.

REACHABILITY FLOOR. A pair is drawn only if its rate clears
net_common.REACHABLE_BPS, currently 256 bit/s = one AES-256 key per second
(Nweke et al., arXiv:2306.15031, define the required rate as key length times
refresh rate). The floor is not neutral between the protocols: CV's curve falls
vertically at its limit because eps_b/(eta*T) diverges, so CV's map barely
changes with it, while DV's tail decays gently and DV loses links faster as the
floor rises. The floor is printed, put in the figure subtitle and written into
the output filename, so two runs at different floors cannot be confused or
silently overwrite each other. Override per run:

    QKD_REACHABLE_BPS=20   python network.py     # Jouguet CV field test
    QKD_REACHABLE_BPS=1e4  python network.py     # carrier-grade

Vary the scenario from the command line:
    python network.py --n 20 --area 30 --seed 7
    python network.py --n 40 --area 15           # denser, shorter links
    python network.py --n 15 --area 60 --save out.png

Engines (dv_rate / cv_rate) are imported from the sensitivity module so the
network stage uses the SAME parameters and conventions as the point-to-point and
sensitivity stages. Set QKD_ENGINE_DIR if needed.
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import networkx as nx

from net_common import (place_users, all_pairs, pair_distance_km,
                        link_rate, to_bps, is_reachable, REACHABLE_BPS,
                        CV_REACH_KM, relayed_network, relay_positions,
                        CLOCK_LABEL)

# UNITS: every rate on this figure — edge colours, colourbar, metrics table and
# the printed summary — is in BITS/S, converted with to_bps at the protocol's
# clock (DV and CV matched at 1 GHz). Viability is decided by net_common's
# single reachability floor (REACHABLE_BPS), applied via is_reachable.

USER_C = "#333333"
RELAY_C = "#ff7f0e"   # trusted-node relays

# matplotlib 3.7 deprecated cm.get_cmap and 3.11 removes it. Resolve once here.
CMAP = plt.get_cmap("viridis")


def fmt_floor(bps):
    """Format a rate for a label without forcing it into kbit/s.

    A plain divide-by-1e3 renders the 256 bit/s default as '0.256 kbit/s'.
    """
    if bps >= 1e6:
        return f"{bps/1e6:g} Mbit/s"
    if bps >= 1e3:
        return f"{bps/1e3:g} kbit/s"
    return f"{bps:g} bit/s"


def build_graph(users, protocol):
    """Graph with all users; edges only where `protocol` clears the floor.
    Edge weight = key rate. Returns (G, pos, rates dict, distances dict)."""
    G = nx.Graph()
    for i in range(len(users)):
        G.add_node(i, pos=tuple(users[i]))
    rates, dists = {}, {}
    for i, j in all_pairs(len(users)):
        d = pair_distance_km(users, i, j)
        r = link_rate(d, protocol)          # bits/channel use — viability test
        dists[(i, j)] = d
        if is_reachable(r, protocol):
            r_bps = to_bps(r, protocol)     # stored and plotted in bits/s
            G.add_edge(i, j, rate=r_bps)
            rates[(i, j)] = r_bps
    return G, nx.get_node_attributes(G, "pos"), rates, dists


def draw_map(ax, users, protocol, title, vmin, vmax):
    G, pos, rates, dists = build_graph(users, protocol)
    n = len(users)
    n_pairs = n * (n - 1) // 2

    # colour edges by rate (log scale), thicker = higher rate
    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    if G.number_of_edges() > 0:
        edges = list(G.edges())
        ecolors = [CMAP(norm(G[u][v]["rate"])) for u, v in edges]
        widths = [0.4 + 2.2 * norm(G[u][v]["rate"]) for u, v in edges]
        nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=ecolors,
                               width=widths, alpha=0.7, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=USER_C, node_size=55, ax=ax)

    # networkx disables ticks; re-enable so the km axes are visible
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    ax.set_axis_on()

    # coverage = reachable pairs / ALL pairs. (This previously called
    # net_common.coverage on the surviving links only, which is the fraction of
    # connected links that are connected — always 1.0.)
    n_conn = len(rates)
    cov = n_conn / n_pairs if n_pairs else 0.0
    rvals = np.array(list(rates.values())) if rates else np.array([])
    metrics = dict(
        success_rate=cov,                                   # fraction of pairs with a key
        n_connected=n_conn,
        # avg over ALL pairs (failed links count as zero) — avoids survivorship bias
        avg_key_rate=float(rvals.sum()) / n_pairs if n_pairs else 0.0,
        min_key_rate=float(rvals.min()) if len(rvals) else 0.0,  # over connected links
        max_key_rate=float(rvals.max()) if len(rvals) else 0.0,  # over connected links
        total_key_rate=float(rvals.sum()) if len(rvals) else 0.0,
    )
    ax.set_title(f"{title}\n{n_conn}/{n_pairs} links "
                 f"({cov*100:.0f}% reachable)", fontsize=11)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2, lw=0.5)
    return metrics


def draw_map_relayed(ax, users, protocol, title, vmin, vmax, reach_km):
    """Like draw_map but rescues failed links with midpoint-chain relays,
    drawing the relays and split links. Returns metrics incl. relay count."""
    net = relayed_network(users, protocol, reach_km)
    # relayed_network returns bits/channel use; convert to bits/s to match draw_map
    rates = {k: to_bps(v, protocol) for k, v in net["rates"].items()}
    relay_pts = net["relay_pts"]
    n = len(users)
    n_pairs = n * (n - 1) // 2

    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)

    # draw each rescued/direct link along its relay chain
    for (i, j), r in rates.items():
        pts = ([tuple(users[i])] + relay_positions(users, i, j, reach_km)
               + [tuple(users[j])])
        c = CMAP(norm(r))
        w = 0.4 + 2.2 * norm(r)
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, color=c, lw=w, alpha=0.7, zorder=2)

    ax.scatter(users[:, 0], users[:, 1], c=USER_C, s=55, zorder=3)
    if relay_pts:
        rp = np.array(relay_pts)
        ax.scatter(rp[:, 0], rp[:, 1], marker="s", c=RELAY_C, s=45,
                   zorder=4, edgecolors="k", linewidths=0.4,
                   label=f"{net['total_relays']} relays")
        ax.legend(loc="upper left", fontsize=9, framealpha=0.9)

    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    cov = len(rates) / n_pairs if n_pairs else 0.0
    ax.set_title(f"{title}  (+relays)\n{len(rates)}/{n_pairs} links "
                 f"({cov*100:.0f}% reachable), {net['total_relays']} relays",
                 fontsize=11)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2, lw=0.5)
    rvals = np.array(list(rates.values())) if rates else np.array([])
    return dict(
        success_rate=cov, n_connected=len(rates),
        # avg over ALL pairs (failed links count as zero) — avoids survivorship bias
        avg_key_rate=float(rvals.sum()) / n_pairs if n_pairs else 0.0,
        min_key_rate=float(rvals.min()) if len(rvals) else 0.0,
        max_key_rate=float(rvals.max()) if len(rvals) else 0.0,
        total_key_rate=float(rvals.sum()) if len(rvals) else 0.0,
        total_relays=net["total_relays"],
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=20, help="Number of users")
    p.add_argument("--area", type=float, default=30.0, help="Area side (km)")
    p.add_argument("--seed", type=int, default=None, help="Random seed")
    p.add_argument("--cv", choices=["heterodyne", "homodyne"],
                   default="heterodyne", help="CV detection scheme to map")
    p.add_argument("--save", type=str, default=None, help="Save path (png/pdf)")
    p.add_argument("--relays", action="store_true",
                   help="Rescue failed CV links with midpoint-chain trusted relays")
    args = p.parse_args()

    if args.seed is None:
        args.seed = int.from_bytes(os.urandom(4), "big") % 100000
        print(f"Seed: {args.seed}")

    print(f"reachability floor: {fmt_floor(REACHABLE_BPS)} "
          f"({REACHABLE_BPS/256:.2f} AES-256 keys/s)   clocks: {CLOCK_LABEL}")

    users = place_users(args.n, area_km=args.area, seed=args.seed)
    cv_proto = "cv_hom" if args.cv == "homodyne" else "cv_het"

    # shared colour scale across both maps (so colours are comparable)
    all_rates = []
    for proto in ("dv", cv_proto):
        for i, j in all_pairs(args.n):
            r = link_rate(pair_distance_km(users, i, j), proto)
            if is_reachable(r, proto):
                all_rates.append(to_bps(r, proto))
    if all_rates:
        vmin, vmax = min(all_rates), max(all_rates)
    else:
        # Fully disconnected layout. Anchor the low end on the floor rather than a
        # hardcoded 1e3, so a link sitting just above a sub-kbit/s floor cannot
        # fall off the bottom of the colour scale.
        vmin, vmax = REACHABLE_BPS, 1e9

    plt.rcParams.update({"font.family": "serif", "font.size": 12})
    fig = plt.figure(figsize=(18, 11))
    # gridspec: top row = two maps + slim colourbar; bottom row = metrics panel
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.045],
                          height_ratios=[1, 0.32],
                          wspace=0.18, hspace=0.28,
                          left=0.06, right=0.94, top=0.90, bottom=0.06)
    axDV = fig.add_subplot(gs[0, 0])
    axCV = fig.add_subplot(gs[0, 1], sharex=axDV, sharey=axDV)
    cax = fig.add_subplot(gs[0, 2])           # dedicated colourbar axis
    axT = fig.add_subplot(gs[1, :])           # metrics table panel
    axT.axis("off")

    statsDV = draw_map(axDV, users, "dv", "DV — decoy BB84", vmin, vmax)
    if args.relays:
        statsCV = draw_map_relayed(axCV, users, cv_proto,
                                   f"CV — GG02 {args.cv}", vmin, vmax, CV_REACH_KM)
    else:
        statsCV = draw_map(axCV, users, cv_proto, f"CV — GG02 {args.cv}", vmin, vmax)

    # shared colourbar in its own axis (no overlap with the plots)
    sm = cm.ScalarMappable(norm=mcolors.LogNorm(vmin=vmin, vmax=vmax),
                           cmap=CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("Link key rate (bits / s)")
    # Mark the floor on the colourbar: everything below it was dropped, not drawn
    # faintly, so the reader can see where the cut sits relative to the scale.
    if vmin <= REACHABLE_BPS <= vmax:
        cbar.ax.axhline(REACHABLE_BPS, color="k", lw=1.2, ls=":")

    # network geometry (same for both)
    dists = np.array([pair_distance_km(users, i, j)
                      for i, j in all_pairs(args.n)])
    avg_dist = float(dists.mean())
    total_fibre = float(dists.sum())

    # ---- which network performed better here? Total key rate over ALL pairs.
    # (Failed links count as zero, so a protocol cannot win by dropping hard
    #  links — this avoids the survivorship bias of averaging over survivors.)
    winner = "DV" if statsDV["total_key_rate"] > statsCV["total_key_rate"] else "CV"
    basis = "total key rate"

    # ---- network performance metrics as a proper table in the bottom panel ----
    def fmt(v, e=False):
        return f"{v:.2e}" if e else f"{v:.0f}"
    col_labels = ["Metric", "DV (BB84)", f"CV ({args.cv})"]
    cell_rows = [
        ["Reachable pairs", f"{statsDV['n_connected']}", f"{statsCV['n_connected']}"],
        ["Reachability", f"{statsDV['success_rate']*100:.0f}%", f"{statsCV['success_rate']*100:.0f}%"],
        ["Avg rate, all pairs (bit/s)", fmt(statsDV['avg_key_rate'], True), fmt(statsCV['avg_key_rate'], True)],
        ["Min rate, connected (bit/s)", fmt(statsDV['min_key_rate'], True), fmt(statsCV['min_key_rate'], True)],
        ["Max rate, connected (bit/s)", fmt(statsDV['max_key_rate'], True), fmt(statsCV['max_key_rate'], True)],
        ["Total key rate (bit/s)", fmt(statsDV['total_key_rate'], True), fmt(statsCV['total_key_rate'], True)],
        ["Trusted relays", "0", f"{statsCV.get('total_relays', 0)}"],
    ]
    tbl = axT.table(cellText=cell_rows, colLabels=col_labels,
                    cellLoc="center", colLoc="center", loc="center",
                    colWidths=[0.26, 0.20, 0.20])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(11)
    tbl.scale(1, 1.5)
    # style header
    for c in range(3):
        tbl[0, c].set_facecolor("#dfe6ee")
        tbl[0, c].set_text_props(weight="bold")
    # highlight the winning protocol's column
    win_col = 1 if winner == "DV" else 2
    for r in range(1, len(cell_rows) + 1):
        tbl[r, win_col].set_facecolor("#e7f4e4")

    geo = (f"N={args.n} users, {args.area:.0f}x{args.area:.0f} km   |   "
           f"avg pair dist {avg_dist:.1f} km, total fibre {total_fibre:.0f} km   |   "
           f"Better here: {winner} (by {basis})")
    axT.set_title(geo, fontsize=11, pad=8)

    fig.suptitle(f"Metro network: DV vs CV   (seed={args.seed}, "
                 f"{fmt_floor(REACHABLE_BPS)} floor)",
                 fontsize=15, weight="bold")

    print(f"\n{'metric':<16}{'DV':>14}{'CV':>14}")
    for label, key, e in [("reachability", "success_rate", False),
                          ("avg rate bit/s", "avg_key_rate", True),
                          ("total bit/s", "total_key_rate", True)]:
        dv = statsDV[key] * (100 if key == "success_rate" else 1)
        cv = statsCV[key] * (100 if key == "success_rate" else 1)
        fs = "{:>14.2e}" if e else "{:>13.0f}%"
        print(f"{label:<16}" + fs.format(dv) + fs.format(cv))
    print(f"Better here: {winner} (by {basis})")

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    # Floor in the filename: a 256 bit/s run and a 1 kbit/s run are different
    # figures on the same layout and must not overwrite each other.
    fname = args.save or (
        f"network_N{args.n}_area{args.area:.0f}_seed{args.seed}"
        f"_floor{REACHABLE_BPS:.0f}bps"
        f"{'_relays' if args.relays else ''}.png")
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"Saved: {path}")
    plt.show()


if __name__ == "__main__":
    main()