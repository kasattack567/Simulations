"""
topo_bom.py -- Bill-of-materials cost model: DV vs CV QKD networks.

WHY NOT UNIT COSTS
Hedges-style unit costing (transmitter = 1, receiver = 1, trusted node = +2)
works when the two protocols use EQUIVALENT hardware -- a detector is a
detector. DV and CV do not: DV needs single-photon detectors with cryogenic
cooling, CV needs room-temperature coherent receivers. Collapsing that into one
number requires a price ratio, and no such ratio is published. So this model
does not produce a cost number as its primary output.

WHAT IT PRODUCES INSTEAD
  Stage 1 -- BILL OF MATERIALS (money-free). What must actually be bought,
             counted, per protocol per topology. Fully determined by the
             topology and the relay analysis; no price assumptions whatsoever.
             An operator applies their own procurement costs to these counts.

  Stage 2 -- SINGLE-PARAMETER BREAK-EVEN. The two bills of materials differ in
             exactly one qualitative respect: DV requires cryogenic receiver
             chains, CV requires more of everything else. So the comparison
             reduces to ONE dimensionless number:

                 rho = (cost of one DV cryogenic receiver chain)
                       / (cost of one CV coherent receiver chain)

             We solve for the rho at which the two networks cost the same.
             Nothing else is assumed. This mirrors the methodology of
             Karavias et al. (ONDM 2025), who sweep relative detector cost
             rather than asserting a value.

COUNTING RULES (each is a stated modelling choice)
  * A QKD link is a span terminated by a transmitter at one end and a receiver
    at the other. Inserting r relays into an edge creates r+1 sub-links, so
      links = E + K      (E = topology edges, K = relays for full coverage)
    Bidirectional operation doubles transmitters and receivers equally and
    therefore cancels in the comparison; unidirectional counts are reported.
  * Transmitters = receivers = links (one of each per link).
  * CRYOSTAT SHARING: one cryocooler serves all SNSPD channels co-located at a
    site, so cryocoolers = sites hosting receivers = N + K, not one per
    receiver. This is the Cambridge/BT co-location economy (arXiv:2110.15005)
    and it favours DV -- included so the model is not biased toward CV.
  * Secure relay sites = K. Site hardening is protocol-independent, so it is
    reported as a count, not folded into a cost.
  * Fibre is identical for both protocols on the same topology and cancels.
  * OPEX: cryocoolers draw ~1.5 kW continuously (arXiv:2006.00411, verified in
    arXiv:1303.6381); room-temperature receivers draw single-digit watts.
    Reported as installed cooling power, not converted to money.

USAGE
  python Networks/topo_bom.py
  python Networks/topo_bom.py --rho 20     # evaluate a specific ratio
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Topology data. N nodes, E edges (from topo_loader), and relays needed
# for full pair coverage (from topo_realscale_relays). EDIT to match your
# script output.
# ------------------------------------------------------------------
TOPOLOGIES = {
    #              N,    E,   K_DV, K_CV
    "TATANID":   (142,  180,   45,  816),
    "SAGO":      (18,    17,    0,   38),
    "GERMANY50": (50,    88,   12,  381),
    "CESNET":    (12,    19,    0,   54),
    "LAYER42":   (6,      7,   49,  356),
    "RNPBRAZIL": (10,    12,   39,  309),
}

CRYO_POWER_KW = 1.5      # per cryocooler, continuous
CV_RX_POWER_W = 10.0     # per room-temperature coherent receiver


def bom(n, e, k):
    """Bill of materials for one protocol on one topology."""
    links = e + k
    return dict(links=links, tx=links, rx=links,
                sites=n + k, relay_sites=k)


def break_even_rho(n, e, k_dv, k_cv, share_cryostats=True):
    """Ratio rho = (DV cryo receiver chain) / (CV coherent receiver chain)
    at which the two networks cost the same.

    Cost in units of one CV receiver chain, counting only what differs:
        C_DV = (E + K_DV) * rho          [DV receivers, cryo]
             + (K_DV - K_CV) * s_site    [relative site hardening]
        C_CV = (E + K_CV) * 1

    Everything symmetric (transmitters, key management, fibre) is counted
    identically in both and omitted -- it cannot affect the ratio.
    Site hardening is protocol-independent per relay, so it enters only
    through the DIFFERENCE in relay counts, and is expressed in CV-receiver
    units as s_site (swept separately below; set 0 here for the pure
    hardware ratio)."""
    n_rx_dv = e + k_dv
    n_rx_cv = e + k_cv
    if n_rx_dv == 0:
        return np.nan
    return n_rx_cv / n_rx_dv


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rho", type=float, default=None,
                    help="evaluate at a specific DV/CV receiver cost ratio")
    args = ap.parse_args()

    names = list(TOPOLOGIES)

    # ---------------- Stage 1: bill of materials ----------------
    print("=" * 78)
    print("STAGE 1 -- BILL OF MATERIALS (no price assumptions)")
    print("=" * 78)
    print(f"{'topology':>11} {'proto':>5} {'links':>6} {'Tx':>6} {'Rx':>6} "
          f"{'cryostats':>10} {'secure sites':>13}")
    boms = {}
    for name in names:
        n, e, k_dv, k_cv = TOPOLOGIES[name]
        b_dv = bom(n, e, k_dv)
        b_cv = bom(n, e, k_cv)
        boms[name] = (b_dv, b_cv)
        print(f"{name:>11} {'DV':>5} {b_dv['links']:6d} {b_dv['tx']:6d} "
              f"{b_dv['rx']:6d} {b_dv['sites']:10d} "
              f"{b_dv['relay_sites']:13d}")
        print(f"{'':>11} {'CV':>5} {b_cv['links']:6d} {b_cv['tx']:6d} "
              f"{b_cv['rx']:6d} {0:10d} {b_cv['relay_sites']:13d}")

    print("\nThe two bills differ in exactly two ways:")
    print("  (a) DV requires cryogenic cooling at every receiver site;")
    print("      CV requires none.")
    print("  (b) CV requires more links, receivers and secure relay sites,")
    print("      because its reach is shorter.")
    print("Everything else is counted identically and cannot affect the "
          "comparison.")

    # ---------------- the two countable resources ----------------
    print("\n" + "=" * 78)
    print("THE TRADE, IN COUNTABLE RESOURCES ONLY")
    print("=" * 78)
    print(f"{'topology':>11} {'cryostats (DV)':>15} {'extra secure sites (CV)':>24}"
          f" {'extra receivers (CV)':>21}")
    for name in names:
        n, e, k_dv, k_cv = TOPOLOGIES[name]
        print(f"{name:>11} {n + k_dv:15d} {k_cv - k_dv:24d} "
              f"{(e + k_cv) - (e + k_dv):21d}")
    print("\nThis is the whole cost question with zero price assumptions:")
    print("  DV buys cryogenics.  CV buys secure real estate and hardware "
          "count.")

    # ---------------- Stage 2: break-even ----------------
    print("\n" + "=" * 78)
    print("STAGE 2 -- SINGLE-PARAMETER BREAK-EVEN")
    print("=" * 78)
    print("rho = (cost of one DV cryogenic receiver chain)")
    print("      / (cost of one CV coherent receiver chain)")
    print("DV is cheaper if the true rho is BELOW the break-even value.\n")
    print(f"{'topology':>11} {'Rx_DV':>7} {'Rx_CV':>7} {'rho*':>8}")
    rhos = []
    for name in names:
        n, e, k_dv, k_cv = TOPOLOGIES[name]
        r = break_even_rho(n, e, k_dv, k_cv)
        rhos.append(r)
        print(f"{name:>11} {e + k_dv:7d} {e + k_cv:7d} {r:8.2f}")

    print(f"\nBreak-even range across topologies: "
          f"{np.nanmin(rhos):.1f}x to {np.nanmax(rhos):.1f}x")
    print("Reading: a DV cryogenic receiver chain may cost up to this many")
    print("times a CV coherent receiver chain and DV is still cheaper overall.")

    if args.rho is not None:
        print(f"\nAt rho = {args.rho}:")
        for name, r in zip(names, rhos):
            print(f"  {name:>11}: {'DV' if args.rho < r else 'CV'} cheaper")

    # ---------------- OPEX ----------------
    print("\n" + "=" * 78)
    print("INSTALLED COOLING POWER (OPEX, reported not costed)")
    print("=" * 78)
    print(f"{'topology':>11} {'DV cooling (kW)':>16} {'CV receivers (kW)':>18}"
          f" {'ratio':>7}")
    for name in names:
        n, e, k_dv, k_cv = TOPOLOGIES[name]
        p_dv = (n + k_dv) * CRYO_POWER_KW
        p_cv = (e + k_cv) * CV_RX_POWER_W / 1000.0
        print(f"{name:>11} {p_dv:16.1f} {p_cv:18.2f} "
              f"{p_dv / p_cv:7.0f}x")

    # ---------------- figure ----------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(len(names))
    w = 0.38

    # panel 1: the two countable resources
    cryos = [TOPOLOGIES[k][0] + TOPOLOGIES[k][2] for k in names]
    extra = [TOPOLOGIES[k][3] - TOPOLOGIES[k][2] for k in names]
    ax1.bar(x - w/2, cryos, w, color="tab:blue", label="Cryostats (DV only)")
    ax1.bar(x + w/2, extra, w, color="tab:red",
            label="Extra secure sites (CV)")
    ax1.set_yscale("log")
    ax1.set_xticks(x); ax1.set_xticklabels(names, rotation=25, ha="right",
                                           fontsize=9)
    ax1.set_ylabel("Count")
    ax1.set_title("What each protocol must buy\n(no price assumptions)")
    ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3, axis="y", which="both")

    # panel 2: break-even ratio
    ax2.bar(x, rhos, color="tab:purple", alpha=0.85)
    ax2.axhline(1.0, color="k", linestyle=":", linewidth=1)
    ax2.set_xticks(x); ax2.set_xticklabels(names, rotation=25, ha="right",
                                           fontsize=9)
    ax2.set_ylabel("Break-even ratio $\\rho^*$")
    ax2.set_title("DV cheaper if its cryogenic receiver costs\n"
                  "less than $\\rho^*$ x a CV coherent receiver")
    ax2.grid(True, alpha=0.3, axis="y")
    for xi, v in zip(x, rhos):
        ax2.text(xi, v, f"{v:.1f}x", ha="center", va="bottom", fontsize=9)

    fig.suptitle("Bill-of-materials cost comparison: DV vs CV on real "
                 "optical topologies", fontsize=13)
    plt.tight_layout()
    plt.savefig("cost_bom.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("\nSaved: cost_bom.png")

    print("""
Caveats for the write-up:
 1. The bill of materials is an OUTPUT of the simulation, not an assumption:
    link and relay counts follow from the topology and the reach limits.
 2. Cryostat sharing is modelled (one cooler per site, not per receiver),
    which favours DV. Without sharing, DV's cryogenic count rises to one per
    receiver and the break-even ratio falls accordingly.
 3. Transmitters are assumed comparable across protocols (both are modulated
    telecom lasers with a QRNG). DV uses intensity modulation for decoy
    states, CV uses IQ modulation; no published cost ratio distinguishes them.
 4. Site hardening is protocol-independent per relay and is therefore reported
    as a count. Folding it in would penalise CV further, since CV needs more
    relay sites.
 5. Not all DV is cryogenic: room-temperature InGaAs receivers exist
    (Toshiba). This model prices the SNSPD tier assumed by the parameter
    baseline; a warm-DV variant would remove the cryogenic column entirely.
 6. CV requires high-speed ADC/DSP whose cost is not publicly quantified;
    it is inside the 'CV coherent receiver chain' the ratio is defined over.
""")


if __name__ == "__main__":
    main()