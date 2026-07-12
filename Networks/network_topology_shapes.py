"""
Topology structures — visual check.

Draws the four canonical topologies (ring, star, tree, mesh) as node/edge
diagrams so the geometry can be verified before any rate analysis. No key-rate
computation and no engine dependency here — this is purely to confirm each
topology is constructed correctly.

Usage:
    python network_topology_shapes.py
    python network_topology_shapes.py --n 12 --area 40
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations


# ---- topology builders (same as network_topology.py) ----
def make_ring(n, area):
    r = area / 2.0
    c = area / 2.0
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pos = np.column_stack([c + r * np.cos(ang), c + r * np.sin(ang)])
    edges = [(i, (i + 1) % n) for i in range(n)]
    return pos, edges


def make_star(n, area):
    r = area / 2.0
    c = area / 2.0
    hub = np.array([[c, c]])
    ang = np.linspace(0, 2 * np.pi, n - 1, endpoint=False)
    spokes = np.column_stack([c + r * np.cos(ang), c + r * np.sin(ang)])
    pos = np.vstack([hub, spokes])
    edges = [(0, k) for k in range(1, n)]
    return pos, edges


def make_tree(n, area):
    c = area / 2.0
    n_branch = max(2, int(round(np.sqrt(n - 1))))
    pos = [[c, c]]
    edges = []
    remaining = n - 1
    ang_b = np.linspace(0, 2 * np.pi, n_branch, endpoint=False)
    r1, r2 = area / 4.0, area / 2.0
    branch_ids = []
    for a in ang_b:
        if remaining <= 0:
            break
        bid = len(pos)
        pos.append([c + r1 * np.cos(a), c + r1 * np.sin(a)])
        edges.append((0, bid)); branch_ids.append((bid, a)); remaining -= 1
    li = 0
    while remaining > 0 and branch_ids:
        bid, a = branch_ids[li % len(branch_ids)]
        jitter = 0.4 * (li // len(branch_ids) + 1) * (-1) ** li
        lid = len(pos)
        pos.append([c + r2 * np.cos(a + 0.2 * jitter),
                    c + r2 * np.sin(a + 0.2 * jitter)])
        edges.append((bid, lid)); remaining -= 1; li += 1
    return np.array(pos, float), edges


def make_mesh(n, area, seed=0):
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, area, size=(n, 2))
    edges = list(combinations(range(n), 2))
    return pos, edges


TOPOS = [("Ring", make_ring), ("Star", make_star),
         ("Tree", make_tree), ("Mesh", make_mesh)]


def draw(ax, pos, edges, name, area):
    # edges
    for a, b in edges:
        ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]],
                color="#377eb8", lw=1.0, alpha=0.55, zorder=1)
    # nodes (hub highlighted for star/tree)
    ax.scatter(pos[:, 0], pos[:, 1], s=70, c="#333333", zorder=3)
    if name in ("Star", "Tree"):
        ax.scatter(pos[0, 0], pos[0, 1], s=140, c="#e41a1c", zorder=4,
                   label="hub")
        ax.legend(loc="upper right", fontsize=9)
    # label nodes with their index (helps verify structure)
    for i, (x, y) in enumerate(pos):
        ax.annotate(str(i), (x, y), textcoords="offset points",
                    xytext=(5, 5), fontsize=8, color="#555")
    ax.set_title(f"{name}  ({len(pos)} nodes, {len(edges)} edges)", fontsize=12)
    ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25, lw=0.5)
    m = area * 0.08
    ax.set_xlim(-m, area + m); ax.set_ylim(-m, area + m)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=12, help="Nodes per topology")
    p.add_argument("--area", type=float, default=40.0, help="Area side (km)")
    p.add_argument("--seed", type=int, default=0, help="Mesh random seed")
    p.add_argument("--save", type=str, default=None)
    args = p.parse_args()

    plt.rcParams.update({"font.family": "serif", "font.size": 11})
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    axes = axes.ravel()

    for ax, (name, builder) in zip(axes, TOPOS):
        if name == "Mesh":
            pos, edges = builder(args.n, args.area, seed=args.seed)
        else:
            pos, edges = builder(args.n, args.area)
        draw(ax, pos, edges, name, args.area)
        print(f"{name}: {len(pos)} nodes, {len(edges)} edges")

    fig.suptitle(f"Topology structures   (N={args.n} nodes, "
                 f"{args.area:.0f}x{args.area:.0f} km)",
                 fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(outdir, exist_ok=True)
    fname = args.save or f"topology_shapes_N{args.n}_area{args.area:.0f}.png"
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=200)
    print(f"\nSaved: {path}")
    plt.show()


if __name__ == "__main__":
    main()
