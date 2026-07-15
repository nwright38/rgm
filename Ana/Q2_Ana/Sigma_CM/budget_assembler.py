#!/usr/bin/env python3
"""Assemble sigma_CM systematic budgets from driver result TTrees."""

import argparse
import json
import re
from pathlib import Path

np = None
uproot = None


def require_budget_modules():
    global np, uproot
    if np is not None and uproot is not None:
        return
    try:
        import numpy as numpy_module
        import uproot as uproot_module
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "budget_assembler.py needs numpy and uproot. Install them in this "
            "Python environment, or run from an environment where they are available."
        ) from exc
    np = numpy_module
    uproot = uproot_module


DIRECTIONS = ["X", "Y", "Z"]


def read(path):
    with uproot.open(path) as f:
        return f["sigmaCM"].arrays(library="np")


def width(values):
    values = np.asarray(values, dtype=float)
    return float(np.std(values, ddof=1)) if values.size > 1 else 0.0


def weighted_spread(values, errors):
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    value_mask = np.isfinite(values)
    if np.count_nonzero(value_mask) < 2:
        return 0.0
    vals = values[value_mask]
    if errors.size != values.size:
        return width(vals)
    errs_all = errors[value_mask]
    err_mask = np.isfinite(errs_all) & (errs_all > 0.0)
    if np.count_nonzero(err_mask) < 2:
        return width(vals)
    vals = vals[err_mask]
    errs = errs_all[err_mask]
    weights = 1.0 / np.square(errs)
    mean = np.average(vals, weights=weights)
    return float(np.sqrt(np.average(np.square(vals - mean), weights=weights)))


def text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def config_value(arr, key, default=None):
    if "configJson" not in arr:
        return np.full(len(next(iter(arr.values()))), default)
    out = []
    for item in arr["configJson"]:
        try:
            out.append(json.loads(text(item)).get(key, default))
        except Exception:
            out.append(default)
    return np.asarray(out)


def integrated_mask(arr):
    n = len(next(iter(arr.values())))
    mask = np.ones(n, dtype=bool)
    if "q2BinIndex" in arr:
        mask &= np.asarray(arr["q2BinIndex"]) < 0
    if "converged" in arr:
        mask &= np.asarray(arr["converged"], dtype=bool)
    return mask


def first_stat(nominal):
    mask = integrated_mask(nominal)
    idx = int(np.where(mask)[0][0]) if np.any(mask) else 0
    return {d: float(nominal[f"sigma{d}ErrHigh"][idx]) for d in DIRECTIONS}


def envelope(arr):
    mask = integrated_mask(arr)
    return {
        d: 0.5 * float(np.max(np.asarray(arr[f"sigma{d}"])[mask]) -
                       np.min(np.asarray(arr[f"sigma{d}"])[mask]))
        if np.any(mask) else 0.0
        for d in DIRECTIONS
    }


def fit_range_spreads(arr):
    mask = integrated_mask(arr)
    return {
        d: weighted_spread(
            np.asarray(arr[f"sigma{d}"], dtype=float)[mask],
            np.asarray(arr[f"sigma{d}ErrHigh"], dtype=float)[mask],
        ) if np.any(mask) else 0.0
        for d in DIRECTIONS
    }


def bound_warnings(arr, label):
    mask = integrated_mask(arr)
    if not np.any(mask):
        return []
    warnings = []
    for d in DIRECTIONS:
        vals = np.asarray(arr[f"sigma{d}"], dtype=float)[mask]
        mins = np.asarray(arr.get("sigmaMin", np.full(len(mask), np.nan)), dtype=float)[mask]
        maxs = np.asarray(arr.get("sigmaMax", np.full(len(mask), np.nan)), dtype=float)[mask]
        if vals.size and np.any(np.isfinite(mins)) and np.any(np.isfinite(maxs)):
            near_min = np.isclose(vals, mins, rtol=0.0, atol=1e-5)
            near_max = np.isclose(vals, maxs, rtol=0.0, atol=1e-5)
            n_bound = int(np.count_nonzero(near_min | near_max))
            if n_bound:
                warnings.append(f"{label}: {n_bound} integrated sigma{d} fits are at/near Minuit bounds")
    return warnings


def closure_injected(arr):
    if "closureInjectedSigma" in arr:
        return np.asarray(arr["closureInjectedSigma"], dtype=float)
    vals = []
    for item in arr.get("status", []):
        m = re.search(r"closure injected sigma=([-+]?[0-9]*\.?[0-9]+)", text(item))
        vals.append(float(m.group(1)) if m else np.nan)
    return np.asarray(vals, dtype=float)


def closure_bias(closure):
    mask = integrated_mask(closure)
    injected = closure_injected(closure)
    if injected.size:
        mask &= np.isfinite(injected)
    if not np.any(mask):
        return {d: 0.0 for d in DIRECTIONS}
    return {
        d: float(np.max(np.abs(np.asarray(closure[f"sigma{d}"], dtype=float)[mask] - injected[mask])))
        for d in DIRECTIONS
    }


def cut_and_bootstrap_widths(cut):
    mask = integrated_mask(cut)
    use_boot = config_value(cut, "useBootstrapWeights", False).astype(bool)
    cut_mask = mask & ~use_boot
    boot_mask = mask & use_boot
    raw = {d: width(np.asarray(cut[f"sigma{d}"])[cut_mask]) for d in DIRECTIONS}
    boot = {d: width(np.asarray(cut[f"sigma{d}"])[boot_mask]) for d in DIRECTIONS}
    sub = {d: float(np.sqrt(max(raw[d] * raw[d] - boot[d] * boot[d], 0.0))) for d in DIRECTIONS}
    return raw, boot, sub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nominal", required=True)
    ap.add_argument("--cut-toys", required=True)
    ap.add_argument("--gcf-toys")
    ap.add_argument("--fit-range", required=True)
    ap.add_argument("--closure", required=True)
    ap.add_argument("--exclude-cuts", action="store_true",
                    help="Do not include stat-subtracted cut toys in total systematic; still report them")
    ap.add_argument("--exclude-gcf", action="store_true",
                    help="Do not include GCF toys in total systematic; still report them")
    ap.add_argument("--exclude-fit-range", action="store_true",
                    help="Do not include fit-range spread in total systematic; still report it")
    ap.add_argument("--exclude-closure", action="store_true",
                    help="Do not include closure bias in total systematic; still report it")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()
    require_budget_modules()

    nominal = read(args.nominal)
    cut = read(args.cut_toys)
    gcf = read(args.gcf_toys) if args.gcf_toys else None
    fit_range = read(args.fit_range)
    closure_raw = read(args.closure)
    cut_raw, bootstrap, cut_sub = cut_and_bootstrap_widths(cut)
    closure = closure_bias(closure_raw)
    fit_range_env = envelope(fit_range)
    fit_range_weighted = fit_range_spreads(fit_range)

    rows = []
    for d in DIRECTIONS:
        gcf_width = width(np.asarray(gcf[f"sigma{d}"])[integrated_mask(gcf)]) if gcf is not None else 0.0
        sources = {
            "statistical": first_stat(nominal)[d],
            "cut_toys_raw": cut_raw[d],
            "data_bootstrap": bootstrap[d],
            "cut_toys_stat_subtracted": cut_sub[d],
            "cut_toys_used_in_total": 0.0 if args.exclude_cuts else cut_sub[d],
            "gcf_toys": gcf_width,
            "gcf_used_in_total": 0.0 if args.exclude_gcf else gcf_width,
            "fit_range_envelope_raw": fit_range_env[d],
            "fit_range_weighted_spread": fit_range_weighted[d],
            "fit_range_used_in_total": 0.0 if args.exclude_fit_range else fit_range_weighted[d],
            "fit_range_envelope": fit_range_weighted[d],
            "closure_bias_uncorrected": closure[d],
            "closure_used_in_total": 0.0 if args.exclude_closure else closure[d],
        }
        syst = np.sqrt(
            sources["cut_toys_used_in_total"] ** 2
            + sources["gcf_used_in_total"] ** 2
            + sources["fit_range_used_in_total"] ** 2
            + sources["closure_used_in_total"] ** 2
        )
        systematics_only = {
            "cut_toys_used_in_total": sources["cut_toys_used_in_total"],
            "gcf_used_in_total": sources["gcf_used_in_total"],
            "fit_range_used_in_total": sources["fit_range_used_in_total"],
            "closure_used_in_total": sources["closure_used_in_total"],
        }
        dominant = max(systematics_only, key=systematics_only.get)
        rows.append({"direction": d, **sources,
                     "dominant_systematic": dominant,
                     "total_systematic": float(syst)})

    out = Path(args.out_prefix)
    out.with_suffix(".json").write_text(json.dumps(rows, indent=2))
    with out.with_suffix(".csv").open("w") as fp:
        keys = list(rows[0])
        fp.write(",".join(keys) + "\n")
        for row in rows:
            fp.write(",".join(str(row[k]) for k in keys) + "\n")
    with out.with_suffix(".tex").open("w") as fp:
        fp.write("\\begin{tabular}{lrrrrr}\\hline\n")
        fp.write("Direction & Stat. & Cuts & GCF & Range & Total sys.\\\\\\hline\n")
        for row in rows:
            fp.write(
                f"{row['direction']} & {row['statistical']:.4f} & {row['cut_toys_used_in_total']:.4f} & "
                f"{row['gcf_used_in_total']:.4f} & {row['fit_range_used_in_total']:.4f} & "
                f"{row['total_systematic']:.4f}\\\\\n"
            )
        fp.write("\\hline\\end{tabular}\n")

    print("Sigma_CM budget summary")
    for row in rows:
        print(
            f"  {row['direction']}: stat={row['statistical']:.5f}, "
            f"cuts(raw/sub/used)={row['cut_toys_raw']:.5f}/{row['cut_toys_stat_subtracted']:.5f}/"
            f"{row['cut_toys_used_in_total']:.5f}, "
            f"boot={row['data_bootstrap']:.5f}, gcf(raw/used)={row['gcf_toys']:.5f}/"
            f"{row['gcf_used_in_total']:.5f}, "
            f"range(raw/used)={row['fit_range_envelope_raw']:.5f}/{row['fit_range_used_in_total']:.5f}, "
            f"closure={row['closure_bias_uncorrected']:.5f}, "
            f"closure_used={row['closure_used_in_total']:.5f}, "
            f"sys={row['total_systematic']:.5f}, dominant={row['dominant_systematic']}"
        )
    warnings = []
    warnings.extend(bound_warnings(fit_range, "fit_range"))
    warnings.extend(bound_warnings(closure_raw, "closure"))
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"  {warning}")


if __name__ == "__main__":
    main()
