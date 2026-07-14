#!/usr/bin/env python3
"""One-command Sigma_CM runner.

Default mode is intentionally quick: nominal stat-only extraction plus plots.
Use --from-hipo for hipo inputs.
Use --full for systematic toys/scans and plots.
Use --profiles when you explicitly want the slower profile-scan ROOT files.
Use --export-budget when you explicitly want JSON/CSV/TeX budget sidecars.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run(cmd):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def exe(build_dir, name):
    exe_name = f"sigmacm_{name}"
    build = Path(build_dir)
    here = Path.cwd()
    script = Path(__file__).resolve()
    source_root = script.parents[3]
    candidates = [
        here / exe_name,
        build / "Ana" / "Q2_Ana" / "Sigma_CM" / exe_name,
        build / exe_name,
        source_root / "build" / "Ana" / "Q2_Ana" / "Sigma_CM" / exe_name,
        source_root / "build" / exe_name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    tried = "\n  ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        f"Could not find {exe_name}. Tried:\n  {tried}\n"
        "Run from the build/Ana/Q2_Ana/Sigma_CM directory, or pass --build-dir."
    )


def helper(name):
    return Path(__file__).resolve().with_name(name)


def run_python(script, args):
    run([sys.executable, helper(script), *args])


def remove_stale(path):
    path = Path(path)
    if path.exists():
        path.unlink()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data")
    ap.add_argument("mc")
    ap.add_argument("out_prefix")
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--from-hipo", action="store_true",
                    help="Read data/mc as hipo files instead of skim ROOT files")
    ap.add_argument("--A", default="4", help="Target mass number for --from-hipo")
    ap.add_argument("--beam-energy", default="5.98636")
    ap.add_argument("--max-events",
                    help="Quick hipo test run over the first N events in each input")
    ap.add_argument("--cache-dir",
                    help="Directory for hipo-to-skim caches; default is OUT_PREFIX_cache")
    ap.add_argument("--seed", default="17")
    ap.add_argument("--cut-range-xy", help="Nominal X/Y fit half-window, e.g. 0.55")
    ap.add_argument("--fit-z-min", help="Nominal Z fit lower edge")
    ap.add_argument("--fit-z-max", help="Nominal Z fit upper edge")
    ap.add_argument("--xy-ranges", default="0.45,0.50,0.55",
                    help="Comma-separated X/Y half-windows for fit-range scan")
    ap.add_argument("--full", action="store_true", help="Run systematic toys/scans too")
    ap.add_argument("--profiles", action="store_true",
                    help="Also run diagnostic profile scans for X/Y/Z")
    ap.add_argument("--export-budget", action="store_true",
                    help="Write budget JSON/CSV/TeX sidecars during --full runs")
    ap.add_argument("--skip-python", action="store_true",
                    help="Run only the C++ executables; skip budget and plot scripts")
    ap.add_argument("--n-toys", default="100")
    ap.add_argument("--n-bootstrap", default="200")
    ap.add_argument("--n-gcf-toys", type=int, default=100,
                    help="Number of GCF toy weight branches to write in hipo MC caches during --full")
    args = ap.parse_args()

    prefix = Path(args.out_prefix)
    plot_dir = prefix.parent / f"{prefix.name}_plots"
    common = ["--seed", args.seed]
    if args.cut_range_xy:
        common.append(f"--cut-range-xy={args.cut_range_xy}")
    if args.fit_z_min:
        common.append(f"--fit-z-min={args.fit_z_min}")
    if args.fit_z_max:
        common.append(f"--fit-z-max={args.fit_z_max}")
    hipo_common = []
    if args.from_hipo:
        hipo_common.extend(["--beam-energy", args.beam_energy])
        if args.max_events:
            hipo_common.extend(["--max-events", args.max_events])

    data_input = Path(args.data)
    mc_input = Path(args.mc)
    if args.from_hipo:
        cache_dir = Path(args.cache_dir) if args.cache_dir else prefix.parent / f"{prefix.name}_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        data_input = cache_dir / "data_skim.root"
        mc_input = cache_dir / "mc_skim.root"
        run([exe(args.build_dir, "make_skim"), args.A, "data", data_input, args.data, *hipo_common])
        mc_skim_cmd = [exe(args.build_dir, "make_skim"), args.A, "mc", mc_input, args.mc, *hipo_common]
        if args.full and args.n_gcf_toys > 0:
            mc_skim_cmd.append(f"--gcf-toys={args.n_gcf_toys}")
        run(mc_skim_cmd)

    nominal = prefix.with_suffix(".nominal.root")
    run([exe(args.build_dir, "extract"), "--from-skim", nominal, data_input, mc_input, *common])

    roots = [nominal]
    budget_json = None
    if args.full:
        cut = prefix.with_suffix(".cut_toys.root")
        gcf = prefix.with_suffix(".gcf_toys.root")
        combined = prefix.with_suffix(".combined_toys.root")
        ranges = prefix.with_suffix(".fit_ranges.root")
        closure = prefix.with_suffix(".closure.root")
        profiles = []
        run([exe(args.build_dir, "run_cut_toys"), data_input, mc_input, cut, *common,
             f"--n-cut-toys={args.n_toys}", f"--n-bootstrap={args.n_bootstrap}"])
        remove_stale(gcf)
        run([exe(args.build_dir, "run_gcf_toys"), data_input, mc_input, gcf, *common])
        run([exe(args.build_dir, "run_combined_toys"), data_input, mc_input, combined, *common,
             f"--n-toys={args.n_toys}"])
        run([exe(args.build_dir, "run_fit_range_scan"), data_input, mc_input, ranges,
             *common, f"--xy-ranges={args.xy_ranges}"])
        run([exe(args.build_dir, "run_closure"), mc_input, closure, *common])
        if args.profiles:
            for axis in range(3):
                prof = prefix.with_suffix(f".profile_axis{axis}.root")
                run([exe(args.build_dir, "run_profile_scan"), data_input, mc_input, prof,
                     *common, f"--axis={axis}", "--scan-min=0.08", "--scan-max=0.26",
                     "--n-points=61"])
                profiles.append(prof)
        budget = prefix.with_suffix(".budget")
        if args.export_budget and not args.skip_python:
            budget_cmd = ["--nominal", nominal, "--cut-toys", cut,
                          "--fit-range", ranges, "--closure", closure,
                          "--out-prefix", budget]
            if gcf.exists():
                budget_cmd.extend(["--gcf-toys", gcf])
            run_python("budget_assembler.py", budget_cmd)
            budget_json = budget.with_suffix(".json")
        roots.extend([cut, combined, ranges, closure, *profiles])
        if gcf.exists():
            roots.append(gcf)

    if not args.skip_python:
        plot_cmd = [*roots, "--out-dir", plot_dir]
        if budget_json:
            plot_cmd.extend(["--budget-json", budget_json])
        run_python("make_plots.py", plot_cmd)
    print(f"Wrote {nominal}")
    if args.skip_python:
        print("Skipped Python budget/plot helpers")
    else:
        print(f"Wrote plots in {plot_dir}")


if __name__ == "__main__":
    main()
