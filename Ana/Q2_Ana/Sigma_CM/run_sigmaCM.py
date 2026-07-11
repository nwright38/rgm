#!/usr/bin/env python3
"""One-command Sigma_CM runner.

Default mode is intentionally quick: nominal stat-only extraction plus plots.
Use --full for cut/GCF/combined toys, profile scans, budget assembly, and plots.
"""

import argparse
import subprocess
from pathlib import Path


def run(cmd):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def exe(build_dir, name):
    top_level = Path(build_dir) / "Ana" / "Q2_Ana" / "Sigma_CM" / f"sigmacm_{name}"
    if top_level.exists():
        return top_level
    return Path(build_dir) / f"sigmacm_{name}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data")
    ap.add_argument("mc")
    ap.add_argument("out_prefix")
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--seed", default="17")
    ap.add_argument("--full", action="store_true", help="Run systematic toys/profile scans too")
    ap.add_argument("--n-toys", default="100")
    ap.add_argument("--n-bootstrap", default="200")
    args = ap.parse_args()

    prefix = Path(args.out_prefix)
    plot_dir = prefix.parent / f"{prefix.name}_plots"
    common = ["--seed", args.seed]

    nominal = prefix.with_suffix(".nominal.root")
    run([exe(args.build_dir, "run_nominal"), args.data, args.mc, nominal, *common])

    roots = [nominal]
    budget_json = None
    if args.full:
        cut = prefix.with_suffix(".cut_toys.root")
        gcf = prefix.with_suffix(".gcf_toys.root")
        combined = prefix.with_suffix(".combined_toys.root")
        ranges = prefix.with_suffix(".fit_ranges.root")
        closure = prefix.with_suffix(".closure.root")
        profiles = []
        run([exe(args.build_dir, "run_cut_toys"), args.data, args.mc, cut, *common,
             f"--n-cut-toys={args.n_toys}", f"--n-bootstrap={args.n_bootstrap}"])
        run([exe(args.build_dir, "run_gcf_toys"), args.data, args.mc, gcf, *common,
             f"--n-toys={args.n_toys}"])
        run([exe(args.build_dir, "run_combined_toys"), args.data, args.mc, combined, *common,
             f"--n-toys={args.n_toys}"])
        run([exe(args.build_dir, "run_fit_range_scan"), args.data, args.mc, ranges, *common])
        run([exe(args.build_dir, "run_closure"), args.mc, closure, *common])
        for axis in range(3):
            prof = prefix.with_suffix(f".profile_axis{axis}.root")
            run([exe(args.build_dir, "run_profile_scan"), args.data, args.mc, prof, *common,
                 f"--axis={axis}", "--scan-min=0.08", "--scan-max=0.26", "--n-points=61"])
            profiles.append(prof)
        budget = prefix.with_suffix(".budget")
        run([Path(__file__).with_name("budget_assembler.py"), "--nominal", nominal,
             "--cut-toys", cut, "--gcf-toys", gcf, "--fit-range", ranges,
             "--closure", closure, "--out-prefix", budget])
        budget_json = budget.with_suffix(".json")
        roots.extend([cut, gcf, combined, ranges, closure, *profiles])

    plot_cmd = [Path(__file__).with_name("make_plots.py"), *roots, "--out-dir", plot_dir]
    if budget_json:
        plot_cmd.extend(["--budget-json", budget_json])
    run(plot_cmd)
    print(f"Wrote {nominal}")
    print(f"Wrote plots in {plot_dir}")


if __name__ == "__main__":
    main()
