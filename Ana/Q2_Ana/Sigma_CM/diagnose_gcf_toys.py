#!/usr/bin/env python3
"""Diagnose GCF toy weight pathologies against fitted sigma_CM results."""

import argparse
import json
from pathlib import Path

np = None
plt = None
PdfPages = None
uproot = None


def require_modules():
    global np, plt, PdfPages, uproot
    if np is not None:
        return
    try:
        import numpy as numpy_module
        import matplotlib.pyplot as pyplot
        from matplotlib.backends.backend_pdf import PdfPages as pdf_pages
        import uproot as uproot_module
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "diagnose_gcf_toys.py needs numpy, matplotlib, and uproot. "
            "Run it in the same Python environment used for plot_sigmaCM.py."
        ) from exc
    np = numpy_module
    plt = pyplot
    PdfPages = pdf_pages
    uproot = uproot_module


def text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def read_tree(path, tree_name):
    with uproot.open(path) as f:
        if tree_name not in f:
            return None
        return f[tree_name].arrays(library="np")


def read_skim(path):
    with uproot.open(path) as f:
        tree_name = "srcTree" if "srcTree" in f else "skim"
        if tree_name not in f:
            raise SystemExit(f"{path} has neither srcTree nor skim")
        tree = f[tree_name]
        branches = set(tree.keys())
        required = [
            "Q2", "xB", "mMiss", "kMiss", "pLead", "thetaLead", "leadRegion",
            "hasRecoil", "pRec", "pcmY", "pcmY_truth", "genWeight",
        ]
        missing = [b for b in required if b not in branches]
        if missing:
            raise SystemExit(f"{path} is missing required branches: {', '.join(missing)}")
        gcf = sorted(b for b in branches if b.startswith("w_gcf_toy_"))
        if not gcf:
            raise SystemExit(f"{path} contains no w_gcf_toy_* branches")
        params = None
        if "gcfToyParams" in f:
            params = f["gcfToyParams"].arrays(library="np")
        return tree.arrays(required + gcf, library="np"), gcf, params


def toy_param_lookup(params):
    if params is None or "branchName" not in params:
        return {}
    out = {}
    n = len(params["branchName"])
    for i in range(n):
        branch = text(params["branchName"][i])
        out[branch] = {
            key: float(params[key][i])
            for key in params
            if key != "branchName"
        }
    return out


def config_for_row(gcf_rows, index):
    if "configJson" not in gcf_rows:
        return {}
    try:
        return json.loads(text(gcf_rows["configJson"][index]))
    except Exception:
        return {}


def selected_mask(mc, cfg):
    n = len(mc["Q2"])
    mask = np.ones(n, dtype=bool)
    mask &= np.asarray(mc["xB"]) > float(cfg.get("xBLower", 1.2))
    mask &= np.asarray(mc["Q2"]) > float(cfg.get("q2Lower", 1.5))
    mask &= np.asarray(mc["Q2"]) < float(cfg.get("q2Upper", 5.0))
    mask &= np.asarray(mc["mMiss"]) > float(cfg.get("mMissLower", 0.65))
    mask &= np.asarray(mc["mMiss"]) < float(cfg.get("mMissUpper", 1.10))
    mask &= np.asarray(mc["kMiss"]) > float(cfg.get("kMissLower", 0.3))
    mask &= np.asarray(mc["pLead"]) > float(cfg.get("pLeadLower", 1.0))
    mask &= np.asarray(mc["hasRecoil"]).astype(bool)
    mask &= np.asarray(mc["pRec"]) > float(cfg.get("pRecLower", 0.3))
    if bool(cfg.get("requirePcmLtPrel", False)) and "pCM" in mc and "pRel" in mc:
        mask &= np.asarray(mc["pCM"]) < np.asarray(mc["pRel"])

    lead_mode = str(cfg.get("leadMode", "FD")).upper()
    lead_region = np.asarray(mc["leadRegion"])
    theta = np.asarray(mc["thetaLead"])
    fd_value = int(cfg.get("fdLeadRegionValue", 2000))
    cd_value = int(cfg.get("cdLeadRegionValue", 4000))
    fd = (lead_region == fd_value) & (theta < float(cfg.get("thetaFDUpper", 37.0)))
    cd = (lead_region == cd_value) & (theta > float(cfg.get("thetaCDLower", 45.0)))
    if lead_mode == "FD":
        mask &= fd
    elif lead_mode == "CD":
        mask &= cd
    elif lead_mode == "BOTH":
        mask &= fd | cd

    if not bool(cfg.get("integratedQ2", True)):
        edges = cfg.get("q2Edges", [])
        idx = int(cfg.get("q2BinIndex", -1))
        if 0 <= idx and idx + 1 < len(edges):
            mask &= np.asarray(mc["Q2"]) >= float(edges[idx])
            mask &= np.asarray(mc["Q2"]) < float(edges[idx + 1])
    return mask


def weighted_rms(values, weights):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights != 0.0)
    if not np.any(mask):
        return float("nan")
    vals = values[mask]
    w = np.abs(weights[mask])
    mean = np.average(vals, weights=w)
    return float(np.sqrt(np.average((vals - mean) ** 2, weights=w)))


def effective_entries(weights):
    weights = np.asarray(weights, dtype=float)
    weights = weights[np.isfinite(weights)]
    if weights.size == 0:
        return 0.0
    sw = np.sum(np.abs(weights))
    sw2 = np.sum(weights * weights)
    return float(sw * sw / sw2) if sw2 > 0.0 else 0.0


def row_mask(gcf_rows):
    n = len(gcf_rows["sigmaY"])
    mask = np.ones(n, dtype=bool)
    if "q2BinIndex" in gcf_rows:
        mask &= np.asarray(gcf_rows["q2BinIndex"]) < 0
    if "converged" in gcf_rows:
        mask &= np.asarray(gcf_rows["converged"], dtype=bool)
    return mask


def build_rows(mc, branches, gcf_rows, toy_params=None):
    mask_rows = row_mask(gcf_rows)
    branch_set = set(branches)
    toy_params = toy_params or {}
    rows = []
    nominal_cache = {}
    for i in np.where(mask_rows)[0]:
        branch = text(gcf_rows["auxWeightBranch"][i]) if "auxWeightBranch" in gcf_rows else ""
        if branch not in branch_set:
            continue
        cfg = config_for_row(gcf_rows, i)
        key = json.dumps(cfg, sort_keys=True)
        if key not in nominal_cache:
            sel = selected_mask(mc, cfg)
            base = np.asarray(mc["genWeight"], dtype=float)[sel]
            nominal_cache[key] = {
                "sel": sel,
                "base": base,
                "n_selected": int(np.count_nonzero(sel)),
                "pcmY_reco_nominal_rms": weighted_rms(np.asarray(mc["pcmY"])[sel], base),
                "pcmY_truth_nominal_rms": weighted_rms(np.asarray(mc["pcmY_truth"])[sel], base),
            }
        cached = nominal_cache[key]
        sel = cached["sel"]
        aux = np.asarray(mc[branch], dtype=float)[sel]
        base = cached["base"]
        total = base * aux
        finite_aux = aux[np.isfinite(aux)]
        abs_aux = np.abs(finite_aux)
        if abs_aux.size == 0:
            continue
        reco_rms = weighted_rms(np.asarray(mc["pcmY"])[sel], total)
        truth_rms = weighted_rms(np.asarray(mc["pcmY_truth"])[sel], total)
        row = {
            "branch": branch,
            "sigmaY": float(gcf_rows["sigmaY"][i]),
            "sigmaYErr": float(gcf_rows["sigmaYErrHigh"][i]) if "sigmaYErrHigh" in gcf_rows else float("nan"),
            "chi2ndf": float(gcf_rows["chi2"][i] / gcf_rows["ndf"][i]) if "chi2" in gcf_rows and "ndf" in gcf_rows and gcf_rows["ndf"][i] > 0 else float("nan"),
            "n_selected": cached["n_selected"],
            "aux_mean": float(np.mean(finite_aux)),
            "aux_std": float(np.std(finite_aux)),
            "aux_abs_max": float(np.max(abs_aux)),
            "aux_abs_p95": float(np.percentile(abs_aux, 95)),
            "aux_abs_p99": float(np.percentile(abs_aux, 99)),
            "aux_abs_p999": float(np.percentile(abs_aux, 99.9)),
            "aux_frac_gt2": float(np.mean(abs_aux > 2.0)),
            "aux_frac_gt5": float(np.mean(abs_aux > 5.0)),
            "aux_frac_gt10": float(np.mean(abs_aux > 10.0)),
            "neff_total": effective_entries(total),
            "neff_frac": effective_entries(total) / cached["n_selected"] if cached["n_selected"] else 0.0,
            "pcmY_reco_rms_ratio": reco_rms / cached["pcmY_reco_nominal_rms"] if cached["pcmY_reco_nominal_rms"] > 0 else float("nan"),
            "pcmY_truth_rms_ratio": truth_rms / cached["pcmY_truth_nominal_rms"] if cached["pcmY_truth_nominal_rms"] > 0 else float("nan"),
        }
        for key, value in toy_params.get(branch, {}).items():
            row[f"param_{key}"] = value
        rows.append(row)
    return rows


def scatter_page(pdf, rows, xkey, xlabel, logx=False):
    x = np.asarray([r[xkey] for r in rows], dtype=float)
    y = np.asarray([r["sigmaY"] for r in rows], dtype=float)
    c = np.asarray([r["neff_frac"] for r in rows], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        return
    plt.figure(figsize=(7.0, 5.0))
    sc = plt.scatter(x[mask], y[mask], c=c[mask], s=36, alpha=0.85, cmap="viridis")
    plt.colorbar(sc, label="effective MC fraction")
    if logx:
        positive = x[mask] > 0.0
        if np.any(positive):
            plt.xscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(r"$\hat{\sigma}_{CM,Y}$ [GeV/c]")
    plt.title(f"sigmaY vs {xlabel}")
    plt.tight_layout()
    pdf.savefig()
    plt.close()


def summary_page(pdf, rows, out_csv):
    sigma = np.asarray([r["sigmaY"] for r in rows], dtype=float)
    order = np.argsort(sigma)[::-1]
    lines = [
        "GCF toy sigmaY diagnostic",
        "",
        f"toy rows used: {len(rows)}",
        f"sigmaY mean/std: {np.mean(sigma):.5f} / {np.std(sigma, ddof=1):.5f}",
        f"sigmaY min/median/max: {np.min(sigma):.5f} / {np.median(sigma):.5f} / {np.max(sigma):.5f}",
        "",
        "largest sigmaY toys:",
        "branch                 sigmaY   p99|w|  max|w|  Neff/N  truthY rms ratio",
    ]
    for idx in order[:10]:
        r = rows[int(idx)]
        lines.append(
            f"{r['branch']:<20} {r['sigmaY']:.5f}  {r['aux_abs_p99']:.3g}  "
            f"{r['aux_abs_max']:.3g}  {r['neff_frac']:.3g}  {r['pcmY_truth_rms_ratio']:.3g}"
        )
    param_keys = sorted(k for k in rows[0] if k.startswith("param_"))
    if param_keys:
        lines += ["", "GCF toy parameters found:", ", ".join(k[6:] for k in param_keys)]
    else:
        lines += ["", "No gcfToyParams tree found in MC skim."]
    lines += ["", f"CSV written to: {out_csv}"]
    plt.figure(figsize=(8.3, 6.2))
    plt.axis("off")
    plt.text(0.02, 0.98, "\n".join(lines), va="top", family="monospace", fontsize=10)
    pdf.savefig()
    plt.close()


def write_csv(rows, path):
    keys = list(rows[0])
    with Path(path).open("w") as fp:
        fp.write(",".join(keys) + "\n")
        for row in rows:
            fp.write(",".join(str(row[k]) for k in keys) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mc_skim", help="MC skim ROOT file containing w_gcf_toy_* branches")
    ap.add_argument("gcf_toys", help="gcf_toys.root result file")
    ap.add_argument("out_pdf", help="multi-page output PDF")
    ap.add_argument("--out-csv", help="optional CSV path; default is PDF stem + .csv")
    args = ap.parse_args()
    require_modules()

    mc, branches, params = read_skim(args.mc_skim)
    gcf_rows = read_tree(args.gcf_toys, "sigmaCM")
    if gcf_rows is None:
        raise SystemExit(f"{args.gcf_toys} is missing sigmaCM tree")
    rows = build_rows(mc, branches, gcf_rows, toy_param_lookup(params))
    if not rows:
        raise SystemExit("No matching integrated GCF toy rows were found")

    out_pdf = Path(args.out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    out_csv = Path(args.out_csv) if args.out_csv else out_pdf.with_suffix(".csv")
    write_csv(rows, out_csv)

    with PdfPages(out_pdf) as pdf:
        summary_page(pdf, rows, out_csv)
        scatter_page(pdf, rows, "aux_abs_p99", "99th percentile |toy/nominal weight|", logx=True)
        scatter_page(pdf, rows, "aux_abs_max", "max |toy/nominal weight|", logx=True)
        scatter_page(pdf, rows, "neff_frac", "effective MC entries / selected MC entries")
        scatter_page(pdf, rows, "pcmY_truth_rms_ratio", "truth pcmY weighted RMS / nominal")
        scatter_page(pdf, rows, "pcmY_reco_rms_ratio", "reco pcmY weighted RMS / nominal")
        scatter_page(pdf, rows, "aux_frac_gt5", "fraction of events with |toy/nominal weight| > 5")
        param_labels = {
            "param_sigmaCM": "randomized GCF sigma_CM",
            "param_Cpp0": "Cpp0",
            "param_Cpn0": "Cpn0",
            "param_Cnn0": "Cnn0",
            "param_Cpn1": "Cpn1",
            "param_TN": "TN",
            "param_TNN": "TNN",
        }
        p_keys = ["param_sigmaCM", "param_Cpp0", "param_Cpn0", "param_Cnn0", "param_Cpn1",
                  "param_P00", "param_P01", "param_P10", "param_P11",
                  "param_P20", "param_P21", "param_P30", "param_P31",
                  "param_TN", "param_TNN"]
        for key in p_keys:
            if key in rows[0]:
                scatter_page(pdf, rows, key, param_labels.get(key, key[6:]))

    print(f"Wrote {out_pdf}")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
