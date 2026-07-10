#!/usr/bin/env python3
"""Assemble sigma_CM systematic budgets from driver result TTrees."""

import argparse
import json
from pathlib import Path

import numpy as np
import uproot


DIRECTIONS = ["X", "Y", "Z"]


def read(path):
    with uproot.open(path) as f:
        return f["sigmaCM"].arrays(library="np")


def width(values):
    values = np.asarray(values, dtype=float)
    return float(np.std(values, ddof=1)) if values.size > 1 else 0.0


def first_stat(nominal):
    return {d: float(nominal[f"sigma{d}ErrHigh"][0]) for d in DIRECTIONS}


def envelope(arr):
    return {d: 0.5 * float(np.max(arr[f"sigma{d}"]) - np.min(arr[f"sigma{d}"])) for d in DIRECTIONS}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nominal", required=True)
    ap.add_argument("--cut-toys", required=True)
    ap.add_argument("--gcf-toys", required=True)
    ap.add_argument("--fit-range", required=True)
    ap.add_argument("--closure", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    nominal = read(args.nominal)
    cut = read(args.cut_toys)
    gcf = read(args.gcf_toys)
    fit_range = read(args.fit_range)
    closure = read(args.closure)

    rows = []
    for d in DIRECTIONS:
        sources = {
            "statistical": first_stat(nominal)[d],
            "cut_toys_raw": width(cut[f"sigma{d}"]),
            "gcf_toys": width(gcf[f"sigma{d}"]),
            "fit_range_envelope": envelope(fit_range)[d],
            "closure_bias_uncorrected": float(np.max(np.abs(closure[f"sigma{d}"] - np.mean(closure[f"sigma{d}"])))),
        }
        syst = np.sqrt(
            sources["cut_toys_raw"] ** 2
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
