"""
Run all parameter-sensitivity scripts. By default saves PNGs into
./figures_sensitivity/ ; pass --show to open interactive windows instead.

    python run_all_sensitivity.py                    # save to ./figures_sensitivity
    python run_all_sensitivity.py --output-dir X     # save to X
    python run_all_sensitivity.py --show            # interactive
    python run_all_sensitivity.py --distance 20     # sweep at 20 km instead of 50

The default test distance (50 km) sits at the DV/CV crossover under the Wang
Eq. 12 noise model. Re-run with --distance to check that conclusions are not an
artefact of that choice; 20 km (CV-favourable) and 80 km (DV-favourable) are the
useful contrast points.
"""
import argparse, os, subprocess, sys

SCRIPTS = ["s1_detector_eff.py", "s2_reconciliation.py", "s3_fibre_atten.py",
           "s4_detector_noise.py", "s5_channel_noise.py", "s6_noise_floor.py"]

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=str, default=None)
    p.add_argument("--distance", type=float, default=None,
                   help="Test distance in km (default: sens_common.L_KM = 50)")
    p.add_argument("--show", action="store_true",
                   help="Open interactive windows instead of saving")
    args = p.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out = None if args.show else (args.output_dir or
                                  os.path.join(here, "figures_sensitivity"))
    if out:
        os.makedirs(out, exist_ok=True)
        print(f"Figures -> {out}")

    for i, s in enumerate(SCRIPTS, 1):
        cmd = [sys.executable, os.path.join(here, s)]
        if out:
            cmd += ["--output-dir", out]
        if args.distance is not None:
            cmd += ["--distance", str(args.distance)]
        print(f"[{i}/{len(SCRIPTS)}] {s}")
        subprocess.run(cmd, check=False)
    print("Done.")
