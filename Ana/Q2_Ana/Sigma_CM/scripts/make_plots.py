#!/usr/bin/env python3
"""Basic sigma_CM result plots from driver TTrees."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


def arrays(path, tree="sigmaCM"):
    with uproot.open(path) as f:
        return f[tree].arrays(library="np")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root_files", nargs="+")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    for root in args.root_files:
        arr = arrays(root)
        stem = Path(root).stem
        for d in "XYZ":
            plt.figure()
            plt.hist(arr[f"sigma{d}"], bins=30, histtype="stepfilled", alpha=0.6)
            plt.xlabel(rf"$\sigma_{{CM,{d}}}$ [GeV/c]")
            plt.ylabel("Fits")
            plt.tight_layout()
            plt.savefig(out / f"{stem}_sigma{d}_distribution.png", dpi=160)
            plt.close()

        if "q2BinIndex" in arr and len(arr["q2BinIndex"]) > 1:
            mask = arr["q2BinIndex"] >= 0
            if np.any(mask):
                plt.figure()
                for d in "XYZ":
                    plt.errorbar(
                        arr["q2BinIndex"][mask],
                        arr[f"sigma{d}"][mask],
                        yerr=arr[f"sigma{d}ErrHigh"][mask],
                        marker="o",
                        linestyle="-",
                        label=d,
                    )
                plt.xlabel(r"$Q^2$ bin index")
                plt.ylabel(r"$\sigma_{CM}$ [GeV/c]")
                plt.legend()
                plt.tight_layout()
                plt.savefig(out / f"{stem}_sigma_vs_q2.png", dpi=160)
                plt.close()


if __name__ == "__main__":
    main()
