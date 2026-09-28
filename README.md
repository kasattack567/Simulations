# CV vs DV QKD for a realistic network deployment regime

Simulation code for an MSc Physics dissertation (University College London) comparing
discrete-variable (DV) decoy-state BB84 with continuous-variable (CV) GG02 quantum key
distribution. The comparison runs point to point, across synthetic networks, and across
six real optical backbone topologies from TopologyBench at their native geographic
scale, finishing with a trusted-node cost model and a Monte Carlo uncertainty analysis.

The full dissertation is included as `Thesis.pdf`.

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [How the code is organised](#how-the-code-is-organised)
- [Running each stage](#running-each-stage)
- [Configuration](#configuration)
- [Outputs](#outputs)
- [Data](#data)
- [Documentation files](#documentation-files)
- [Citation](#citation)

## Installation

Requires **Python 3.9 to 3.12**. The CV engine (`qosst-skr`) does not yet support
Python 3.13 or NumPy 2, so `requirements.txt` pins NumPy below 2. Tested with Python 3.12.

1. Clone the repository.

   ```bash
   git clone https://github.com/<your-username>/<repo-name>.git
   cd <repo-name>
   ```

2. Create and activate a virtual environment (recommended).

   ```bash
   python -m venv .venv
   source .venv/bin/activate        # Windows: .venv\Scripts\activate
   ```

3. Install the dependencies.

   ```bash
   pip install -r requirements.txt
   ```

   This installs the two key-rate engines, [`qosst-skr`](https://pypi.org/project/qosst-skr/)
   for CV and [`tno.quantum.communication.qkd_key_rate`](https://pypi.org/project/tno.quantum.communication.qkd_key_rate/)
   for DV, together with NumPy, SciPy, Matplotlib, pandas, NetworkX and openpyxl.

## Quick start

Run all commands from the repository root. To reproduce the main point-to-point result:

```bash
python cv_hetro.py
```

This takes about 30 seconds and writes `figure2_bits_per_channel_use.png`,
`figure3_equal_clock_1GHz.png` and `figure4_bits_per_second.png`.

To reproduce the headline real-topology results:

```bash
python topo_realscale_relays.py     # trusted relay counts per topology
python topo_cost.py                 # break-even cost ratios
```

## How the code is organised

Every script draws its parameters and engine calls from two shared modules, so the
baseline is defined once and cannot drift between stages.

| Module | Role |
|---|---|
| `sens_common.py` | Locked parameter baseline, the DV and CV engine wrappers (`dv_rate`, `cv_rate`), clock rates and plotting helpers. Not run directly. |
| `net_common.py` | Network layer built on `sens_common`: link rates, reach distances, the reachability floor and user placement. Not run directly. |
| `topo_loader.py` | Loads the six TopologyBench `.xlsx` files. Not run directly. |
| `cost_common.py` | Cost model machinery (Karavias et al., ONDM 2025). Not run directly. |

The locked baseline and its sources are documented in `param.md`.

## Running each stage

The stages below follow the order of the dissertation. Runtimes are approximate, on a
standard laptop.

### 1. Point-to-point comparison

| Command | What it does | Runtime |
|---|---|---|
| `python cv_hetro.py` | Key rate vs distance for DV, CV heterodyne and CV homodyne, in bits per channel use and bits/s | ~30 s |

### 2. Parameter sensitivity

Six one-at-a-time sweeps, evaluated at a fixed test distance (default 50 km, near the
crossover).

```bash
python run_all_sensitivity.py                  # all six, saved to figures_sensitivity/
python run_all_sensitivity.py --distance 20    # CV-favourable test distance
python run_all_sensitivity.py --distance 80    # DV-favourable test distance
python run_all_sensitivity.py --show           # open interactive windows instead
```

Each sweep can also be run on its own with `python s1_detector_eff.py` and so on.

| Script | Parameter swept |
|---|---|
| `s1_detector_eff.py` | Detector efficiency |
| `s2_reconciliation.py` | CV β and DV f_EC (twin axes, as they are different quantities) |
| `s3_fibre_atten.py` | Fibre attenuation, including hollow-core fibre values |
| `s4_detector_noise.py` | DV dark counts and CV electronic noise |
| `s5_channel_noise.py` | DV QBER and CV receiver-side excess noise ε_b |
| `s6_noise_floor.py` | CV excess-noise floor ε_a + ε_l |

Runtime is about 25 seconds for all six.

### 3. Synthetic networks

Users are placed at random in a square area. These scripts write to `figures/`.

| Command | What it does | Runtime |
|---|---|---|
| `python network_topology_shapes.py` | Draws ring, star, tree and mesh layouts as a visual check | ~1 s |
| `python network.py --n 20 --area 30 --seed 61715` | DV and CV link maps on the same user layout | ~10 s |
| `python network_sweep.py` | Coverage and total key rate as the area grows | several minutes |
| `python network_users.py` | Average key rate as the number of users grows | several minutes |
| `python network_topology.py` | Ring, star, tree and mesh compared across sizes | several minutes |
| `python network_tradeoff.py` | Rate vs trusted relays once both protocols reach full coverage | several minutes |

Most accept `--n`, `--runs`, `--areas` and `--save`. Run any script with `--help` for its
full options. For a fast test run, reduce the work, for example
`python network_sweep.py --runs 1 --areas 10 40`.

### 4. Real topologies

| Command | What it does | Runtime |
|---|---|---|
| `python topo_maps.py` | Maps of the six topologies from their node coordinates | ~2 s |
| `python topo_direct.py` | Link coverage, pair coverage and average rate with no relays | ~10 s |
| `python topo_realscale_relays.py` | Trusted relays needed per topology, and relayed key rate | ~15 s |

`topo_realscale_relays.py` takes `--cv heterodyne|homodyne` and
`--criterion reach|rate`. Use `reach` (the default) for the reach comparison and `rate`
only for costing; the docstring explains why.

### 5. Cost model

```bash
python topo_cost.py                  # break-even analysis, both clock regimes
python topo_cost_summary.py          # break-even ratio mapped over demand and site cost
python topo_cost_summary.py --slices # also write cost_slices.png
```

Both write to `figures_cost/`. `topo_cost.py` takes about 30 seconds and
`topo_cost_summary.py --slices` about 2 minutes. Key options for `topo_cost.py` are
`--cmin` (key-rate demand per link, default 10 Mbit/s), `--site` (trusted-node site
cost, default 5) and `--all-figures`.

### 6. Uncertainty analysis

Latin hypercube sampling over the deployed parameter bands, with standardised rank
regression coefficients for sensitivity. Each sample takes around 1 to 2 seconds, and
the run is resumable: stopping and restarting continues from the last completed row.

```bash
python uncertainty_analysis.py --samples 128                     # sample, then analyse
python uncertainty_analysis.py --samples 1024                    # continue to 1024
python uncertainty_analysis.py --analyse-only                    # re-analyse existing samples
python uncertainty_analysis.py --floor-bps 1000                  # different usability floor
```

Output goes to `uncertainty/` by default (change with `--output-dir`). The dissertation
results use 1024 samples, and those samples are included as `uncertainty_samples.csv`.
To re-analyse them without re-sampling:

```bash
mkdir -p uncertainty
cp uncertainty_samples.csv uncertainty/
python uncertainty_analysis.py --analyse-only
```

**Important:** delete `uncertainty/uncertainty_samples.csv` after any change to the
parameter baseline. The resume logic counts rows only, so it will otherwise extend a
file generated under the old baseline.

## Configuration

Parameters are set in `sens_common.py` and should be changed there only. A few
settings can be overridden per run with environment variables:

| Variable | Default | Effect |
|---|---|---|
| `QKD_REACHABLE_BPS` | `256` | Minimum key rate (bit/s) for a link to count as usable. 256 bit/s is one AES-256 key per second. |
| `QKD_DV_CLOCK_HZ` | `1e9` | DV clock rate in the network scripts |
| `QKD_COST_DV_CLOCK_HZ` | `1e9` | DV clock rate in the cost model |
| `QKD_COST_CV_CLOCK_HZ` | `100e6` | CV deployed clock rate in the cost model |
| `QKD_ENGINE_DIR` | repository root | Folder containing `sens_common.py`, only needed if you move files into subfolders |

Example:

```bash
QKD_REACHABLE_BPS=20 python network_sweep.py     # Jouguet CV field-test floor
QKD_REACHABLE_BPS=1e4 python network_sweep.py    # carrier-grade floor
```

The floor is written into output filenames, so runs at different floors do not
overwrite each other.

## Outputs

| Location | Written by |
|---|---|
| Repository root | `cv_hetro.py`, `topo_direct.py` |
| `figures/` | `network*.py`, `topo_maps.py`, `topo_realscale_relays.py` |
| `figures_sensitivity/` | `run_all_sensitivity.py` |
| `figures_cost/` | `topo_cost.py`, `topo_cost_summary.py` (figures, CSV tables and `tables_cost.md`) |
| `uncertainty/` | `uncertainty_analysis.py` |

The figures committed to this repository are the versions used in the dissertation.
Running the scripts regenerates them.

## Data

The six topologies (TATANID, SAGO, GERMANY50, CESNET, LAYER42, RNPBRAZIL) are taken
from TopologyBench (Matzner et al., *JOCN* 17, 7–27, 2024; Zenodo record 13921775).
They were selected by PCA and k-means clustering over all 105 TopologyBench networks.
Each `.xlsx` file holds a node sheet (ID, latitude, longitude, location, country) and
an edge sheet (source, destination, length in km).

## Documentation files

| File | Contents |
|---|---|
| `param.md` | Locked parameter baseline and the source of each value |
| `engines_summary.md` | How the two engines are used and the corrections applied to each |
| `methodology.md` | Methodological decisions and their justification |
| `AUDIT.md` | Record of a full audit of the code and results |
| `roadmap.md` | Project overview and status |

## Modelling scope

Both engines are asymptotic; finite-key effects are excluded from both arms. Link
distances in the synthetic networks are straight-line, which is a lower bound on real
fibre length. The dissertation discusses these and other limitations in full.

## Citation

If you use this code, please cite the dissertation:

```
Bedford, B. (2026). CV vs DV QKD for a realistic network deployment regime.
MSc dissertation, University College London.
```

## Licence

<!-- Add a licence, e.g. MIT, and a LICENSE file. Without one, others have no legal right to reuse the code. -->
