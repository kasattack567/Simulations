"""
Draw the 6 selected real optical topologies (TopologyBench) as maps.

A verification/overview figure: each network drawn from its real node
coordinates and edges, so you can confirm they loaded correctly (USA100 should
look like the USA, HIBERNIAUK like the UK/Ireland, etc.) and see their structure
before the QKD rate analysis.

Usage:
    python topo_maps.py                 # real geographic scale
    python topo_maps.py --scale metro   # metro-rescaled (40 km span)
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from topo_loader import TOPOLOGY_FILES, load_topology, scaled_xy, edge_lengths


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scale", choices=["real", "metro"], default="real")
    p.add_argument("--metro-span", type=float, default=40.0)
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    plt.rcParams.update({"font.family": "serif", "font.size": 10})
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    axes = axes.ravel()

    for ax, name in zip(axes, TOPOLOGY_FILES):
        topo = load_topology(name)
        xy = scaled_xy(topo, args.scale, args.metro_span)
        elen = edge_lengths(topo, args.scale, args.metro_span)

        # edges coloured by length (longer = redder, to hint at QKD difficulty)
        lmax = elen.max() if len(elen) else 1.0
        for (i, j), L in zip(topo["edges"], elen):
            frac = L / lmax
            ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]],
                    color=plt.cm.plasma(0.15 + 0.7 * frac), lw=1.2, alpha=0.75,
                    zorder=1)
        ax.scatter(xy[:, 0], xy[:, 1], s=35, c="#222", zorder=3)

        unit = "km"
        ax.set_title(f"{name}\n{len(topo['node_ids'])} nodes, "
                     f"{len(topo['edges'])} edges, "
                     f"links {elen.min():.0f}-{elen.max():.0f} {unit}",
                     fontsize=11)
        ax.set_xlabel(f"x ({unit})"); ax.set_ylabel(f"y ({unit})")
        ax.set_aspect("equal"); ax.grid(True, alpha=0.25, lw=0.5)

    scale_txt = ("real geographic scale" if args.scale == "real"
                 else f"metro-rescaled to {args.metro_span:.0f} km span")
    fig.suptitle(f"Six real optical topologies (TopologyBench) — {scale_txt}",
                 fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or f"topo_maps_{args.scale}.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"Saved: {path}")
    plt.show()


if __name__ == "__main__":
    main()
