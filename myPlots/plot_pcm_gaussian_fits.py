"""
plot_pcm_gaussian_fits.py

Integrated pcmx/pcmy distributions with Gaussian fits overlaid.

Reads the same Main_Figs_Binned-style ROOT files used by plot_1D.py and
ExtractFitQuantities.C. The preferred input is the nominal TH1D under
`hists/nominal/nominal_pcmx_epp` and `hists/nominal/nominal_pcmy_epp`;
older top-level histogram names are also tried.

Usage:
    python plot_pcm_gaussian_fits.py Data_He.root
    python plot_pcm_gaussian_fits.py Data_He.root --out-dir pdf/pcm_fits
    python plot_pcm_gaussian_fits.py Data_He.root --method stddev
"""

from __future__ import division

import argparse
import math
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())

import matplotlib.pyplot as plt
import numpy as np
import uproot
from scipy.optimize import curve_fit


PLOTS = [
    {
        "key": "pcmx",
        "task": "pcmx_epp",
        "xlabel": r"$\vec{p}_{C.M.} \cdot \hat{v}_{x} [GeV]$",
        "fit_range": (-0.25, 0.25),
    },
    {
        "key": "pcmy",
        "task": "pcmy_epp",
        "xlabel": r"$\vec{p}_{C.M.} \cdot \hat{v}_{y} [GeV]$",
        "fit_range": (-0.25, 0.25),
    },
]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("data_file", help="Main_Figs_Binned/BuildGraphs ROOT output")
    p.add_argument("--out-dir", default="pdf/pcm_fits")
    p.add_argument("--format", default="pdf", choices=["pdf", "png"])
    p.add_argument("--method", default="gaussian", choices=["gaussian", "stddev"],
                   help="Use a Gaussian fit sigma or histogram std dev in the range.")
    p.add_argument("--xlim", nargs=2, type=float, default=(-0.75, 0.75),
                   metavar=("MIN", "MAX"))
    p.add_argument("--fit-range-pcmx", nargs=2, type=float, default=None,
                   metavar=("MIN", "MAX"))
    p.add_argument("--fit-range-pcmy", nargs=2, type=float, default=None,
                   metavar=("MIN", "MAX"))
    return p.parse_args()


def gaussian(x, norm, mu, sigma):
    z = (x - mu) / sigma
    return norm / (sigma * math.sqrt(2.0 * math.pi)) * np.exp(-0.5 * z * z)


def candidate_hist_names(task):
    return [
        "hists/nominal/nominal_" + task,
        "nominal_" + task,
        task,
        "h_" + task,
    ]


def read_hist(root_file, task):
    tried = candidate_hist_names(task)
    for name in tried:
        try:
            hist = root_file[name]
        except Exception:
            continue

        values, edges = hist.to_numpy(flow=False)
        variances = hist.variances(flow=False)
        errors = np.sqrt(np.clip(variances, 0.0, None))
        return np.asarray(values, dtype=float), np.asarray(edges, dtype=float), errors, name

    raise KeyError("Could not find histogram for %s. Tried: %s"
                   % (task, ", ".join(tried)))


def initial_guess(x, y, fit_mask):
    xf = x[fit_mask]
    yf = y[fit_mask]
    positive = yf > 0.0
    if np.any(positive):
        weights = yf[positive]
        xp = xf[positive]
        mu = np.average(xp, weights=weights)
        var = np.average((xp - mu) ** 2, weights=weights)
        sigma = math.sqrt(max(var, 1.0e-6))
        norm = np.sum(weights) * np.mean(np.diff(x)) if len(x) > 1 else np.sum(weights)
    else:
        mu = 0.5 * (xf[0] + xf[-1])
        sigma = 0.25 * (xf[-1] - xf[0])
        norm = 1.0
    return norm, mu, sigma


def fit_gaussian(centers, counts, errors, fit_min, fit_max):
    fit_mask = ((centers >= fit_min) & (centers <= fit_max) &
                np.isfinite(counts) & np.isfinite(errors) & (counts > 0.0))
    if np.count_nonzero(fit_mask) < 3:
        raise RuntimeError("Need at least 3 populated bins in the fit range")

    x = centers[fit_mask]
    y = counts[fit_mask]
    yerr = errors[fit_mask]
    yerr = np.where(yerr > 0.0, yerr, 1.0)

    p0 = initial_guess(centers, counts, fit_mask)
    lower = [0.0, fit_min, 1.0e-6]
    upper = [np.inf, fit_max, fit_max - fit_min]
    popt, pcov = curve_fit(
        gaussian, x, y, p0=p0, sigma=yerr, absolute_sigma=True,
        bounds=(lower, upper), maxfev=20000)

    perr = np.sqrt(np.diag(pcov)) if pcov is not None else np.full(3, np.nan)
    return popt, perr


def range_stddev(centers, counts, fit_min, fit_max):
    mask = ((centers >= fit_min) & (centers <= fit_max) &
            np.isfinite(counts) & (counts > 0.0))
    if np.count_nonzero(mask) == 0:
        raise RuntimeError("Need at least 1 populated bin in the std-dev range")

    x = centers[mask]
    w = counts[mask]
    mean = np.average(x, weights=w)
    variance = np.average((x - mean) ** 2, weights=w)
    stddev = math.sqrt(max(variance, 0.0))

    # ROOT's TH1::GetStdDevError approximation for a weighted distribution.
    # This is useful as a visual guide; systematic treatment still belongs
    # in the dedicated extraction macros.
    n_eff = (np.sum(w) ** 2) / np.sum(w ** 2) if np.sum(w ** 2) > 0.0 else 0.0
    stddev_err = stddev / math.sqrt(2.0 * max(n_eff - 1.0, 1.0))
    return mean, stddev, stddev_err


def draw_plot(counts, edges, errors, hist_name, spec, fit_range, args):
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    fit_min, fit_max = fit_range

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.errorbar(centers, counts, yerr=errors, xerr=0.5 * widths,
                color="black", linestyle="", marker="o", markersize=4,
                linewidth=1.0, label="Data")
    ax.stairs(counts, edges, color="black", linewidth=1.2)

    ax.axvspan(fit_min, fit_max, color="firebrick", alpha=0.08, linewidth=0)

    peak = np.max(counts + errors)
    if args.method == "gaussian":
        popt, perr = fit_gaussian(centers, counts, errors, fit_min, fit_max)
        xfit = np.linspace(fit_min, fit_max, 400)
        yfit = gaussian(xfit, *popt)
        ax.plot(xfit, yfit, color="firebrick", linewidth=2.0,
                label="Gaussian fit")
        peak = max(peak, np.max(yfit))
        text = (r"$\mu = %.4f \pm %.4f$ GeV" "\n"
                r"$\sigma = %.4f \pm %.4f$ GeV"
                % (popt[1], perr[1], popt[2], perr[2]))
        result = (popt[1], perr[1], popt[2], perr[2])
    else:
        mean, stddev, stddev_err = range_stddev(centers, counts, fit_min, fit_max)
        ax.axvline(mean, color="firebrick", linewidth=1.6, linestyle="--",
                   label="Range mean")
        text = (r"$\langle x\rangle_{range} = %.4f$ GeV" "\n"
                r"$StdDev_{range} = %.4f \pm %.4f$ GeV"
                % (mean, stddev, stddev_err))
        result = (mean, 0.0, stddev, stddev_err)

    ax.set_xlim(args.xlim)
    ax.set_ylim(0.0, 1.2 * peak if peak > 0.0 else 1.0)
    ax.set_xlabel(spec["xlabel"], fontsize=14)
    ax.set_ylabel("Counts", fontsize=14)
    ax.text(0.04, 0.92, r"$(e,e^{\prime}pp)$", transform=ax.transAxes,
            fontsize=18, ha="left", va="top")
    ax.text(
        0.04, 0.80,
        text,
        transform=ax.transAxes, fontsize=11, ha="left", va="top")
    ax.legend(loc="upper right", frameon=False)
    ax.set_title(hist_name, fontsize=10)
    fig.tight_layout()

    out_path = os.path.join(args.out_dir, "%s_epp_gaussian_fit.%s"
                            % (spec["key"], args.format))
    if args.method == "stddev":
        out_path = os.path.join(args.out_dir, "%s_epp_range_stddev.%s"
                                % (spec["key"], args.format))
    fig.savefig(out_path)
    plt.close(fig)
    return out_path, result


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    fit_range_overrides = {
        "pcmx": args.fit_range_pcmx,
        "pcmy": args.fit_range_pcmy,
    }

    root_file = uproot.open(args.data_file)
    for spec in PLOTS:
        counts, edges, errors, hist_name = read_hist(root_file, spec["task"])
        fit_range = fit_range_overrides[spec["key"]] or spec["fit_range"]
        out_path, result = draw_plot(
            counts, edges, errors, hist_name, spec, fit_range, args)
        if args.method == "gaussian":
            print(
                "Wrote %s  mu=%.6g +/- %.3g, sigma=%.6g +/- %.3g"
                % (out_path, result[0], result[1], result[2], result[3]))
        else:
            print(
                "Wrote %s  range mean=%.6g, stddev=%.6g +/- %.3g"
                % (out_path, result[0], result[2], result[3]))


if __name__ == "__main__":
    main()
