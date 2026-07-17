#!/usr/bin/env python3
"""Simple plotting for Sigma_CM ROOT result files."""

import argparse
import json
from pathlib import Path

np = None
plt = None
uproot = None


def require_plot_modules():
    global np, plt, uproot
    if np is not None and plt is not None and uproot is not None:
        return
    try:
        import numpy as numpy_module
        import matplotlib.pyplot as pyplot
        import uproot as uproot_module
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "make_plots.py needs numpy, matplotlib, and uproot. Install them in this "
            "Python environment, or run from an environment where they are available."
        ) from exc
    np = numpy_module
    plt = pyplot
    uproot = uproot_module


DIRECTIONS = ("X", "Y", "Z")
MIN_TOY_ENTRIES = 25


def read_tree(path, tree):
    with uproot.open(path) as f:
        if tree not in f:
            return None
        return f[tree].arrays(library="np")


def systematic_lookup(path):
    if not path:
        return None
    rows = json.loads(Path(path).read_text())
    return {row["direction"]: float(row.get("total_systematic", 0.0)) for row in rows}


def read_budget(path):
    if not path:
        return None
    return json.loads(Path(path).read_text())


def savefig(out, name):
    plt.tight_layout()
    plt.savefig(out / f"{name}.pdf")
    plt.close()


def remove_stale_scan_plots(out, stem):
    for suffix in ("sigma_vs_fit_window", "chi2ndf_vs_fit_window", "z_window_vs_fit_window"):
        path = out / f"{stem}_{suffix}.pdf"
        if path.exists():
            path.unlink()


def _as_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def q2_bin_centers(arr, mask):
    if "q2BinLow" in arr and "q2BinHigh" in arr:
        low = np.asarray(arr["q2BinLow"])[mask]
        high = np.asarray(arr["q2BinHigh"])[mask]
        return 0.5 * (low + high)

    if "configJson" in arr:
        centers = []
        for cfg_text, bin_index in zip(np.asarray(arr["configJson"])[mask],
                                       np.asarray(arr["q2BinIndex"])[mask]):
            cfg = json.loads(_as_text(cfg_text))
            edges = cfg.get("q2Edges", [])
            idx = int(bin_index)
            if idx < 0 or idx + 1 >= len(edges):
                centers.append(float("nan"))
            else:
                centers.append(0.5 * (float(edges[idx]) + float(edges[idx + 1])))
        return np.asarray(centers, dtype=float)

    return 0.5 * (np.asarray(arr["q2Lower"])[mask] + np.asarray(arr["q2Upper"])[mask])


def integrated_sigmas(arr):
    if arr is None or len(arr.get("sigmaX", [])) == 0:
        return {}
    use = 0
    if "q2BinIndex" in arr:
        integrated = np.where(np.asarray(arr["q2BinIndex"]) < 0)[0]
        if integrated.size:
            use = int(integrated[0])
    out = {}
    for d in DIRECTIONS:
        value = float(arr[f"sigma{d}"][use])
        if np.isfinite(value) and value != 0.0:
            out[d] = value
    return out


def nominal_sigmas_from_roots(root_files):
    preferred = [p for p in root_files if Path(p).stem == "nominal"]
    for root in preferred + list(root_files):
        arr = read_tree(root, "sigmaCM")
        values = integrated_sigmas(arr)
        if values:
            return values
    return {}


def plot_toy_distributions(arr, stem, out):
    for d in DIRECTIONS:
        values = np.asarray(arr[f"sigma{d}"], dtype=float)
        if values.size < MIN_TOY_ENTRIES:
            continue
        plt.figure(figsize=(6.2, 4.4))
        plt.hist(values, bins=30, histtype="stepfilled", alpha=0.65, color="#4477aa")
        plt.axvline(np.mean(values), color="#cc6677", lw=1.8)
        plt.xlabel(rf"$\hat{{\sigma}}_{{CM,{d}}}$ [GeV/c]")
        plt.ylabel("Toy fits")
        plt.title(f"{stem}: sigma-hat toys, {d}")
        savefig(out, f"{stem}_sigma{d}_toy_distribution")


def plot_budget_sources(rows, out, nominal_sigmas):
    if not rows:
        return
    by_direction = {row["direction"]: row for row in rows}
    directions = [d for d in DIRECTIONS if d in by_direction]
    if nominal_sigmas:
        directions = [d for d in directions if d in nominal_sigmas]
    if not directions:
        return

    source_keys = [
        (("statistical_data", "data_bootstrap"), "data stat"),
        (("statistical_mc_estimated",), "MC stat"),
        (("cut_toys_used_in_total", "cut_toys_stat_subtracted"), "cuts"),
        (("gcf_used_in_total", "gcf_toys"), "GCF"),
        (("fit_range_used_in_total", "fit_range_envelope"), "fit range"),
        (("closure_used_in_total",), "closure"),
        (("total_systematic",), "sys total"),
    ]
    def source_value(row, keys):
        for key in keys:
            if key in row:
                return float(row.get(key, 0.0))
        return 0.0

    values = np.array([
        [source_value(by_direction[d], keys) for keys, _ in source_keys]
        for d in directions
    ])
    if nominal_sigmas:
        denominators = np.array([abs(nominal_sigmas[d]) for d in directions], dtype=float)
        values = 100.0 * values / denominators[:, np.newaxis]

    x = np.arange(len(directions))
    width = min(0.12, 0.75 / max(len(source_keys), 1))
    offsets = (np.arange(len(source_keys)) - 0.5 * (len(source_keys) - 1)) * width
    colors = ["#222222", "#777777", "#66c2a5", "#fc8d62", "#8da0cb", "#a6d854", "#999999"]

    plt.figure(figsize=(7.2, 4.5))
    for i, (_, label) in enumerate(source_keys):
        plt.bar(x + offsets[i], values[:, i], width=width, label=label,
                color=colors[i], alpha=0.88)

    if any("fit_range_envelope_raw" in by_direction[d] for d in directions):
        raw = np.array([
            float(by_direction[d].get("fit_range_envelope_raw", np.nan))
            for d in directions
        ])
        if nominal_sigmas:
            raw = 100.0 * raw / np.array([abs(nominal_sigmas[d]) for d in directions], dtype=float)
        plt.scatter(x + offsets[4], raw, marker="_", s=160, linewidths=2.0,
                    color="#4c4c4c", label="fit range raw envelope")

    plt.xticks(x, directions)
    if nominal_sigmas:
        plt.ylabel(r"Relative uncertainty on $\sigma_{CM}$ [%]")
    else:
        plt.ylabel(r"Uncertainty on $\sigma_{CM}$ [GeV/c]")
    plt.title("uncertainty source budget")
    plt.legend(frameon=False, ncol=3)
    savefig(out, "budget_uncertainty_sources")


def plot_sigma_vs_q2(arr, stem, out, sys, show_integrated_sys):
    if "q2BinIndex" not in arr:
        return
    mask = np.asarray(arr["q2BinIndex"]) >= 0
    if not np.any(mask):
        return
    x = q2_bin_centers(arr, mask)
    finite = np.isfinite(x)
    if not np.any(finite):
        return
    mask_indices = np.where(mask)[0][finite]
    x = x[finite]
    order = np.argsort(x)
    x = x[order]
    plt.figure(figsize=(7.2, 4.8))
    colors = {"X": "#4477aa", "Y": "#228833", "Z": "#cc6677"}
    for d in DIRECTIONS:
        y = np.asarray(arr[f"sigma{d}"])[mask_indices][order]
        stat = np.asarray(arr[f"sigma{d}ErrHigh"])[mask_indices][order]
        if sys is not None and show_integrated_sys:
            total = np.sqrt(stat * stat + sys.get(d, 0.0) ** 2)
            plt.fill_between(x, y - total, y + total, color=colors[d], alpha=0.13, linewidth=0)
        plt.errorbar(x, y, yerr=stat, marker="o", linestyle="-", color=colors[d], label=f"{d} stat")
    plt.xlabel(r"$Q^2$ [GeV$^2$]")
    plt.ylabel(r"$\sigma_{CM}$ [GeV/c]")
    plt.title(f"{stem}: sigma vs Q2")
    plt.legend(frameon=False, ncol=3)
    suffix = "stat_and_integrated_total" if sys is not None and show_integrated_sys else "stat"
    savefig(out, f"{stem}_sigma_vs_q2_{suffix}")


def plot_integrated_summary(arr, stem, out, sys):
    if arr is None or len(arr.get("sigmaX", [])) == 0:
        return
    use = 0
    if "q2BinIndex" in arr:
        integrated = np.where(np.asarray(arr["q2BinIndex"]) < 0)[0]
        if integrated.size:
            use = int(integrated[0])
    labels = list(DIRECTIONS)
    y = np.array([arr[f"sigma{d}"][use] for d in labels], dtype=float)
    stat = np.array([arr[f"sigma{d}ErrHigh"][use] for d in labels], dtype=float)
    x = np.arange(len(labels))
    plt.figure(figsize=(5.8, 4.2))
    if sys is not None:
        total = np.sqrt(stat * stat + np.array([sys.get(d, 0.0) for d in labels]) ** 2)
        plt.errorbar(x - 0.04, y, yerr=total, fmt="none", ecolor="#999999",
                     elinewidth=6, alpha=0.45, label="stat+sys")
    plt.errorbar(x, y, yerr=stat, fmt="o", color="#222222", label="stat")
    plt.xticks(x, labels)
    plt.ylabel(r"$\sigma_{CM}$ [GeV/c]")
    plt.title(f"{stem}: integrated widths")
    plt.legend(frameon=False)
    suffix = "stat_and_total" if sys is not None else "stat"
    savefig(out, f"{stem}_integrated_sigma_{suffix}")


def remove_integrated_total(out, stem):
    path = out / f"{stem}_integrated_sigma_stat_and_total.pdf"
    if path.exists():
        path.unlink()


def plot_closure_response(arr, stem, out):
    if arr is None or "closureInjectedSigma" not in arr:
        return
    injected = np.asarray(arr["closureInjectedSigma"], dtype=float)
    mask = np.isfinite(injected)
    if "q2BinIndex" in arr:
        mask &= np.asarray(arr["q2BinIndex"]) < 0
    if "converged" in arr:
        mask &= np.asarray(arr["converged"], dtype=bool)
    if np.count_nonzero(mask) < 2:
        return

    x = injected[mask]
    order = np.argsort(x)
    x = x[order]
    colors = {"X": "#4477aa", "Y": "#228833", "Z": "#cc6677"}

    plt.figure(figsize=(6.2, 5.0))
    lo = float(np.nanmin(x))
    hi = float(np.nanmax(x))
    for d in DIRECTIONS:
        y = np.asarray(arr[f"sigma{d}"], dtype=float)[mask][order]
        err = np.asarray(arr[f"sigma{d}ErrHigh"], dtype=float)[mask][order]
        plt.errorbar(x, y, yerr=err, marker="o", linestyle="-", color=colors[d], label=d)
        lo = min(lo, float(np.nanmin(y)))
        hi = max(hi, float(np.nanmax(y)))
    pad = 0.05 * max(hi - lo, 1.0e-6)
    plt.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#555555", lw=1.2, ls="--", label="ideal")
    plt.xlabel(r"Injected $\sigma_{CM}$ [GeV/c]")
    plt.ylabel(r"Extracted $\sigma_{CM}$ [GeV/c]")
    plt.title(f"{stem}: closure response")
    plt.legend(frameon=False, ncol=2)
    savefig(out, f"{stem}_closure_injected_vs_extracted")

    plt.figure(figsize=(6.2, 4.4))
    for d in DIRECTIONS:
        y = np.asarray(arr[f"sigma{d}"], dtype=float)[mask][order]
        plt.plot(x, y - x, marker="o", linestyle="-", color=colors[d], label=d)
    plt.axhline(0.0, color="#555555", lw=1.2, ls="--")
    plt.xlabel(r"Injected $\sigma_{CM}$ [GeV/c]")
    plt.ylabel(r"Extracted - injected [GeV/c]")
    plt.title(f"{stem}: closure bias")
    plt.legend(frameon=False, ncol=3)
    savefig(out, f"{stem}_closure_bias_vs_injected")


def plot_fit_range_scan(arr, stem, out):
    if "cutRangeXY" not in arr or "chi2" not in arr or "ndf" not in arr:
        remove_stale_scan_plots(out, stem)
        return
    if len(arr.get("sigmaX", [])) < 2:
        remove_stale_scan_plots(out, stem)
        return
    if "q2BinIndex" in arr:
        mask = np.asarray(arr["q2BinIndex"]) < 0
    else:
        mask = np.ones(len(arr["sigmaX"]), dtype=bool)
    if np.count_nonzero(mask) < 2:
        remove_stale_scan_plots(out, stem)
        return

    x = np.asarray(arr["cutRangeXY"], dtype=float)[mask]
    finite = np.isfinite(x)
    if np.count_nonzero(finite) < 2:
        remove_stale_scan_plots(out, stem)
        return
    indices = np.where(mask)[0][finite]
    x = x[finite]
    if np.unique(np.round(x, decimals=8)).size < 2:
        remove_stale_scan_plots(out, stem)
        return
    order = np.argsort(x)
    indices = indices[order]
    x = x[order]

    colors = {"X": "#4477aa", "Y": "#228833", "Z": "#cc6677"}
    plt.figure(figsize=(7.0, 4.6))
    for d in DIRECTIONS:
        y = np.asarray(arr[f"sigma{d}"], dtype=float)[indices]
        err = np.asarray(arr[f"sigma{d}ErrHigh"], dtype=float)[indices]
        plt.errorbar(x, y, yerr=err, marker="o", linestyle="-", color=colors[d], label=d)
    plt.xlabel(r"X/Y fit half-window [GeV/c]")
    plt.ylabel(r"$\sigma_{CM}$ [GeV/c]")
    plt.title(f"{stem}: sigma vs fit window")
    plt.legend(frameon=False, ncol=3)
    savefig(out, f"{stem}_sigma_vs_fit_window")

    ndf = np.asarray(arr["ndf"], dtype=float)[indices]
    chi2 = np.asarray(arr["chi2"], dtype=float)[indices]
    chi2ndf = np.divide(chi2, ndf, out=np.full_like(chi2, np.nan), where=ndf > 0)
    plt.figure(figsize=(6.6, 4.3))
    plt.plot(x, chi2ndf, marker="o", linestyle="-", color="#444444")
    plt.xlabel(r"X/Y fit half-window [GeV/c]")
    plt.ylabel(r"$\chi^2 / ndf$")
    plt.title(f"{stem}: fit quality vs fit window")
    savefig(out, f"{stem}_chi2ndf_vs_fit_window")

    if "fitZMin" in arr and "fitZMax" in arr:
        zmin = np.asarray(arr["fitZMin"], dtype=float)[indices]
        zmax = np.asarray(arr["fitZMax"], dtype=float)[indices]
        plt.figure(figsize=(6.6, 4.3))
        plt.plot(x, zmin, marker="o", linestyle="-", label="Z min", color="#4477aa")
        plt.plot(x, zmax, marker="o", linestyle="-", label="Z max", color="#cc6677")
        plt.xlabel(r"X/Y fit half-window [GeV/c]")
        plt.ylabel(r"Z fit edge [GeV/c]")
        plt.title(f"{stem}: Z window used in scan")
        plt.legend(frameon=False)
        savefig(out, f"{stem}_z_window_vs_fit_window")


def plot_profiles(path, stem, out):
    prof = read_tree(path, "profile")
    if prof is None or len(prof.get("sigma", [])) == 0:
        return
    axes = sorted(set(int(a) for a in prof["axis"]))
    names = {0: "X", 1: "Y", 2: "Z"}
    for axis in axes:
        mask = np.asarray(prof["axis"]) == axis
        x = np.asarray(prof["sigma"])[mask]
        y = np.asarray(prof["chi2"])[mask]
        order = np.argsort(x)
        y = y[order]
        x = x[order]
        plt.figure(figsize=(6.2, 4.4))
        plt.plot(x, y - np.min(y), marker="o", ms=3, color="#4477aa")
        plt.axhline(1.0, color="#cc6677", lw=1.4, ls="--")
        plt.xlabel(r"$\sigma_{CM}$ [GeV/c]")
        plt.ylabel(r"$\Delta\chi^2$")
        plt.title(f"{stem}: profile chi2, {names.get(axis, axis)}")
        savefig(out, f"{stem}_profile_chi2_{names.get(axis, axis)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root_files", nargs="+", help="Nominal, toy, profile, or scan result ROOT files")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--budget-json", help="JSON written by budget_assembler.py")
    ap.add_argument("--show-integrated-sys-on-q2", action="store_true",
                    help=("Draw integrated budget systematics on Q2-binned plots. "
                          "By default Q2 plots stay stat-only because the budget is integrated-only."))
    args = ap.parse_args()
    require_plot_modules()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    budget_rows = read_budget(args.budget_json)
    sys = systematic_lookup(args.budget_json)
    nominal_sigmas = nominal_sigmas_from_roots(args.root_files) if budget_rows else {}
    plot_budget_sources(budget_rows, out, nominal_sigmas)
    if args.budget_json and not args.show_integrated_sys_on_q2:
        print("Q2 plots are stat-only: the supplied budget is integrated-only.")

    stems = [Path(root).stem for root in args.root_files]
    nominal_total_stem = next((stem for stem in stems if "nominal" in stem), stems[0] if stems else "")
    for root in args.root_files:
        arr = read_tree(root, "sigmaCM")
        if arr is None:
            continue
        stem = Path(root).stem
        if sys is None or stem != nominal_total_stem:
            remove_integrated_total(out, stem)
            plot_integrated_summary(arr, stem, out, None)
        else:
            plot_integrated_summary(arr, stem, out, sys)
        plot_sigma_vs_q2(arr, stem, out, sys, args.show_integrated_sys_on_q2)
        plot_fit_range_scan(arr, stem, out)
        plot_toy_distributions(arr, stem, out)
        plot_closure_response(arr, stem, out)
        plot_profiles(root, stem, out)


if __name__ == "__main__":
    main()
