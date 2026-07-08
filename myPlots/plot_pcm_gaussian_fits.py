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
    python plot_pcm_gaussian_fits.py Data_He.root --overlay-pcmxy
    python plot_pcm_gaussian_fits.py Data_He.root --rebin 2
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


LEGEND_FONTSIZE = 12
STATS_FONTSIZE = 11
FIT_LINEWIDTH = 1.6
FIT_ALPHA = 0.65


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
    p.add_argument("--rebin", type=positive_int, default=1,
                   help="Combine this many adjacent bins before fitting/drawing.")
    p.add_argument("--fit-range-pcmx", nargs=2, type=float, default=None,
                   metavar=("MIN", "MAX"))
    p.add_argument("--fit-range-pcmy", nargs=2, type=float, default=None,
                   metavar=("MIN", "MAX"))
    p.add_argument("--overlay-pcmxy", action="store_true",
                   help="Draw pcm_x and pcm_y on one figure.")
    return p.parse_args()


def positive_int(value):
    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive integer")
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


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


def rebin_hist(counts, edges, errors, factor):
    if factor == 1:
        return counts, edges, errors

    counts = np.asarray(counts, dtype=float)
    edges = np.asarray(edges, dtype=float)
    errors = np.asarray(errors, dtype=float)
    variances = np.clip(errors, 0.0, None) ** 2

    rebinned_counts = []
    rebinned_variances = []
    rebinned_edges = [edges[0]]

    for start in range(0, len(counts), factor):
        stop = min(start + factor, len(counts))
        rebinned_counts.append(np.sum(counts[start:stop]))
        rebinned_variances.append(np.sum(variances[start:stop]))
        rebinned_edges.append(edges[stop])

    return (
        np.asarray(rebinned_counts, dtype=float),
        np.asarray(rebinned_edges, dtype=float),
        np.sqrt(np.asarray(rebinned_variances, dtype=float)),
    )


def read_hist_for_plot(root_file, task, rebin):
    counts, edges, errors, hist_name = read_hist(root_file, task)
    counts, edges, errors = rebin_hist(counts, edges, errors, rebin)
    return counts, edges, errors, hist_name


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


def stats_text(mean, mean_err, sigma, sigma_err):
    return "\n".join([
        r"$\mu = %.4f \pm %.4f$ GeV" % (mean, mean_err),
        r"$\sigma = %.4f \pm %.4f$ GeV" % (sigma, sigma_err),
    ])


def stats_box_kwargs():
    return {
        "boxstyle": "round,pad=0.25",
        "facecolor": "white",
        "edgecolor": "none",
        "alpha": 0.82,
    }


def draw_plot(counts, edges, errors, hist_name, spec, fit_range, args):
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    fit_min, fit_max = fit_range

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.errorbar(centers, counts, yerr=errors, xerr=0.5 * widths,
                color="black", linestyle="", marker="o", markersize=4,
                linewidth=1.0, label="Data")

    peak = np.max(counts + errors)
    if args.method == "gaussian":
        popt, perr = fit_gaussian(centers, counts, errors, fit_min, fit_max)
        xfit = np.linspace(fit_min, fit_max, 400)
        yfit = gaussian(xfit, *popt)
        ax.plot(xfit, yfit, color="firebrick", linewidth=FIT_LINEWIDTH,
                alpha=FIT_ALPHA, label="Gaussian fit")
        peak = max(peak, np.max(yfit))
        text = stats_text(popt[1], perr[1], popt[2], perr[2])
        result = (popt[1], perr[1], popt[2], perr[2])
    else:
        mean, stddev, stddev_err = range_stddev(centers, counts, fit_min, fit_max)
        ax.axvline(mean, color="firebrick", linewidth=1.6, linestyle="--",
                   label="Range mean")
        text = stats_text(mean, 0.0, stddev, stddev_err)
        result = (mean, 0.0, stddev, stddev_err)

    ax.set_xlim(args.xlim)
    ax.set_ylim(0.0, 1.2 * peak if peak > 0.0 else 1.0)
    ax.set_xlabel(spec["xlabel"], fontsize=14)
    ax.set_ylabel("Counts", fontsize=14)
    ax.text(
        0.04, 0.92,
        text,
        transform=ax.transAxes, fontsize=STATS_FONTSIZE, ha="left", va="top",
        bbox=stats_box_kwargs())
    ax.legend(loc="upper right", frameon=False, fontsize=LEGEND_FONTSIZE)
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


def fit_result(centers, counts, errors, fit_range, method):
    fit_min, fit_max = fit_range
    if method == "gaussian":
        popt, perr = fit_gaussian(centers, counts, errors, fit_min, fit_max)
        return {
            "mean": popt[1],
            "mean_err": perr[1],
            "sigma": popt[2],
            "sigma_err": perr[2],
            "fit_params": popt,
            "fit_errors": perr,
        }

    mean, stddev, stddev_err = range_stddev(centers, counts, fit_min, fit_max)
    return {
        "mean": mean,
        "mean_err": 0.0,
        "sigma": stddev,
        "sigma_err": stddev_err,
        "fit_params": None,
        "fit_errors": None,
    }


def draw_overlay_pcmxy(root_file, fit_range_overrides, args):
    colors = {"pcmx": "black", "pcmy": "royalblue"}
    markers = {"pcmx": "o", "pcmy": "s"}
    labels = {
        "pcmx": r"$p_{C.M.,x}$",
        "pcmy": r"$p_{C.M.,y}$",
    }

    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    peak = 0.0
    results = []

    for spec in PLOTS:
        counts, edges, errors, hist_name = read_hist_for_plot(
            root_file, spec["task"], args.rebin)
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = np.diff(edges)
        fit_range = fit_range_overrides[spec["key"]] or spec["fit_range"]
        fit_min, fit_max = fit_range
        result = fit_result(centers, counts, errors, fit_range, args.method)
        results.append((spec, hist_name, result))

        color = colors[spec["key"]]
        ax.errorbar(centers, counts, yerr=errors, xerr=0.5 * widths,
                    color=color, linestyle="", marker=markers[spec["key"]],
                    markersize=4, linewidth=1.0, label=labels[spec["key"]])

        if args.method == "gaussian":
            xfit = np.linspace(fit_min, fit_max, 400)
            yfit = gaussian(xfit, *result["fit_params"])
            ax.plot(xfit, yfit, color=color, linewidth=FIT_LINEWIDTH,
                    alpha=FIT_ALPHA)
            peak = max(peak, np.max(yfit))
        else:
            ax.axvline(result["mean"], color=color, linewidth=1.5,
                       linestyle="--", alpha=0.9)

        peak = max(peak, np.max(counts + errors))

    ax.set_xlim(args.xlim)
    ax.set_ylim(0.0, 1.2 * peak if peak > 0.0 else 1.0)
    ax.set_xlabel(r"$p_{C.M.}$ component [GeV]", fontsize=14)
    ax.set_ylabel("Counts", fontsize=14)

    text_lines = []
    for spec, _, result in results:
        short_label = r"$p_x$" if spec["key"] == "pcmx" else r"$p_y$"
        text_lines.append(r"%s $\mu = %.4f \pm %.4f$ GeV"
                          % (short_label, result["mean"],
                             result["mean_err"]))
        text_lines.append(r"%s $\sigma = %.4f \pm %.4f$ GeV"
                          % (short_label, result["sigma"],
                             result["sigma_err"]))
    ax.text(0.04, 0.92, "\n".join(text_lines),
            transform=ax.transAxes, fontsize=STATS_FONTSIZE, ha="left", va="top",
            bbox=stats_box_kwargs())
    ax.legend(loc="upper right", frameon=False, fontsize=LEGEND_FONTSIZE)

    suffix = "gaussian_fit" if args.method == "gaussian" else "range_stddev"
    out_path = os.path.join(args.out_dir, "pcmxy_epp_overlay_%s.%s"
                            % (suffix, args.format))
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path, results


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    fit_range_overrides = {
        "pcmx": args.fit_range_pcmx,
        "pcmy": args.fit_range_pcmy,
    }

    root_file = uproot.open(args.data_file)
    if args.overlay_pcmxy:
        out_path, results = draw_overlay_pcmxy(root_file, fit_range_overrides, args)
        print("Wrote %s" % out_path)
        for spec, hist_name, result in results:
            if args.method == "gaussian":
                print(
                    "%s (%s): mu=%.6g +/- %.3g, sigma=%.6g +/- %.3g"
                    % (spec["key"], hist_name, result["mean"], result["mean_err"],
                       result["sigma"], result["sigma_err"]))
            else:
                print(
                    "%s (%s): range mean=%.6g, stddev=%.6g +/- %.3g"
                    % (spec["key"], hist_name, result["mean"],
                       result["sigma"], result["sigma_err"]))
        return

    for spec in PLOTS:
        counts, edges, errors, hist_name = read_hist_for_plot(
            root_file, spec["task"], args.rebin)
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
