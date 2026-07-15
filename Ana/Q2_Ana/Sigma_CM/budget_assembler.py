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
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()
    require_budget_modules()

    nominal = read(args.nominal)
    cut = read(args.cut_toys)
    gcf = read(args.gcf_toys) if args.gcf_toys else None
    fit_range = read(args.fit_range)
    closure = read(args.closure)
    cut_raw, bootstrap, cut_sub = cut_and_bootstrap_widths(cut)
    closure = closure_bias(closure)
    fit_range_env = envelope(fit_range)

    rows = []
    for d in DIRECTIONS:
        sources = {
            "statistical": first_stat(nominal)[d],
            "cut_toys_raw": cut_raw[d],
            "data_bootstrap": bootstrap[d],
            "cut_toys_stat_subtracted": cut_sub[d],
            "gcf_toys": width(np.asarray(gcf[f"sigma{d}"])[integrated_mask(gcf)]) if gcf is not None else 0.0,
            "fit_range_envelope": fit_range_env[d],
            "closure_bias_uncorrected": closure[d],
        }
        syst = np.sqrt(
            sources["cut_toys_stat_subtracted"] ** 2
            + sources["gcf_toys"] ** 2
            + sources["fit_range_envelope"] ** 2
            + sources["closure_bias_uncorrected"] ** 2
        )
        rows.append({"direction": d, **sources, "total_systematic": float(syst)})

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
                f"{row['direction']} & {row['statistical']:.4f} & {row['cut_toys_raw']:.4f} & "
                f"{row['gcf_toys']:.4f} & {row['fit_range_envelope']:.4f} & "
                f"{row['total_systematic']:.4f}\\\\\n"
            )
        fp.write("\\hline\\end{tabular}\n")


if __name__ == "__main__":
    main()
