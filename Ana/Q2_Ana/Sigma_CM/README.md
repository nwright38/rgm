# Sigma_CM extraction

This directory is intentionally flat. The C++ source, executable drivers, tests,
and optional Python helpers all live directly in `Ana/Q2_Ana/Sigma_CM`.

The fit basis is unchanged from the previous Sigma_CM code: MC templates are
reweighted from `sigma_gen`, the default fit uses one shared scale, and the chi2
is

```text
sum (data - scale * mc)^2 / (data_err^2 + scale^2 * mc_err^2)
```

## 1. Build

From the repository root:

```bash
cmake -S . -B build
cmake --build build --target sigmacm_tests
```

The executables will be under:

```bash
build/Ana/Q2_Ana/Sigma_CM/
```

## 2. Run With C++ Only

This is the simplest path and does not require any Python packages.

Daily hipo-to-result path, from a full repo build with CLAS12/hipo available:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_extract 4 nominal.root data.hipo sim.hipo
```

Fast/debug skim path:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_extract --from-skim nominal.root data.root mc.root
```

The standalone Sigma_CM build can still compile without CLAS12/hipo; in that
case hipo mode prints a clear "build does not include hipo support" error and
`--from-skim` remains available.

Nominal/stat-only extraction, integrated plus Q2-binned:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_nominal data.root mc.root nominal.root
```

`nominal.root` contains both:

```text
sigmaCM/profile TTrees for the new workflow
standard ROOT hists, graphs, best-fit overlays, and canvases
```

The standard plotting objects are written by the new Sigma_CM code. You do not need
to run `Main_sigmaCM_Hists` to get the ROOT plotting surface from a nominal
Sigma_CM extraction.

Profile-chi2 scan for one projection:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_profile_scan data.root mc.root profile_x.root --axis=0
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_profile_scan data.root mc.root profile_y.root --axis=1
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_profile_scan data.root mc.root profile_z.root --axis=2
```

Toy and systematic inputs:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_cut_toys data.root mc.root cut_toys.root --n-cut-toys=100 --n-bootstrap=200
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_gcf_toys data.root mc.root gcf_toys.root --n-toys=100
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_combined_toys data.root mc.root combined_toys.root --n-toys=100
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_fit_range_scan data.root mc.root fit_ranges.root --xy-ranges=0.45,0.50,0.55
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_closure mc.root closure.root
```

Useful common options:

```text
--seed N
--q2-bin N
--aux-weight branchName
--independent-scales
--pcm-lt-prel
--cut-range-xy=0.55
--fit-z-min=-0.5
--fit-z-max=1.0
--beam-energy=5.98636
--max-events=N
```

Histogram comparison ranges:

```text
X chi2 compares bins from -cutRangeXY to +cutRangeXY
Y chi2 compares bins from -cutRangeXY to +cutRangeXY
Z chi2 compares bins from fitZMin to fitZMax
T chi2 compares bins from 0 to sqrt(2)*cutRangeXY
```

With defaults, that is X/Y `[-0.55, 0.55]`, Z `[-0.5, 1.0]`,
and T `[0, 0.778]`.

## 3. Plotting And Budget Scripts

The plotting and budget scripts are in this same directory:

```text
plot_sigmaCM.py       obvious-name plotting script
make_plots.py         same plotting code; kept for compatibility
budget_assembler.py   systematic budget tables
run_sigmaCM.py        convenience wrapper around executables plus scripts
```

`plot_sigmaCM.py` reads the ROOT files made by the C++ executables and writes:

```text
sigma vs Q2 with stat and stat+sys bands
profile-chi2 curves
toy sigma-hat distributions
integrated sigma summaries
```

You can run it with only the nominal ROOT file. In that case it makes the plots
that are possible from nominal/stat-only results. The script uses `uproot` to
read the ROOT outputs directly. `--budget-json` is optional; only add it when
you have explicitly exported a budget sidecar and want stat+sys bands from that
file.

`budget_assembler.py` is an optional exporter. It reads the
nominal/toy/scan/closure ROOT files and writes the systematic budget as JSON,
CSV, and TeX only when you run it.

The data+best-fit overlay canvases with the familiar ROOT object names are
also made by `Main_sigmaCM_Hists`, which remains useful for comparisons to the
historical stage-2 workflow.

They need these Python packages:

```text
numpy
uproot
matplotlib
```

In the current shell I checked, `uproot` is not installed for `python3`, so the
C++ executables can run but the Python plotting scripts will fail until the
Python environment has those packages.

One clean setup is:

```bash
python3 -m venv .venv_sigmacm
source .venv_sigmacm/bin/activate
python -m pip install numpy uproot matplotlib
```

Then run the helpers with that environment active:

```bash
Ana/Q2_Ana/Sigma_CM/plot_sigmaCM.py \
  nominal.root \
  --out-dir plots_nominal_only

Ana/Q2_Ana/Sigma_CM/plot_sigmaCM.py \
  nominal.root cut_toys.root gcf_toys.root combined_toys.root profile_x.root profile_y.root profile_z.root \
  --out-dir plots

Ana/Q2_Ana/Sigma_CM/budget_assembler.py \
  --nominal nominal.root \
  --cut-toys cut_toys.root \
  --gcf-toys gcf_toys.root \
  --fit-range fit_ranges.root \
  --closure closure.root \
  --out-prefix budget
```

The plot command writes:

```text
plots/*.pdf
plots/*.png
```

The optional budget export writes `budget.json`, `budget.csv`, and
`budget.tex`.

## 4. One-Command Wrapper

`run_sigmaCM.py` is just a convenience wrapper around the C++ executables plus
the Python helpers. It can be run from the repository root or from the build
output directory that contains `sigmacm_run_nominal`.

If the Python environment is not ready yet, add `--skip-python`. That runs the
C++ executable steps and skips budget/plot scripts.

Quick nominal run plus plots, from the repository root:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out
```

Quick nominal run with no Python helpers, from `build/Ana/Q2_Ana/Sigma_CM`:

```bash
python3 ../../../../Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out --skip-python
```

Full run with toys, profiles, and plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out --full
```

Add table sidecars only when you want them:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out --full --export-budget
```

Full run with only the C++ ROOT outputs:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out --full --skip-python
```

If you built Sigma_CM directly instead of through the full repo, pass that build
directory:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out --build-dir /tmp/sigmacm_build
```

## 5. Standard ROOT Plotting Output

`sigmacm_run_nominal` now writes the standard ROOT plotting surface into its output
ROOT file, alongside the compact `sigmaCM` result TTree. `Main_sigmaCM_Hists`
is still untouched and remains available for historical comparisons, but the new
Sigma_CM nominal output has the same ROOT objects available for everyday
plotting and inspection.

Important object names include:

```text
sigmacmx_int, sigmacmy_int, sigmacmz_int, sigmacmT_int
sigmacmx_Q2, sigmacmy_Q2, sigmacmz_Q2, sigmacmT_Q2
pcmx_epp, pcmy_epp, pcmz_epp, pcmT_epp
pcmx_epp_fit, pcmy_epp_fit, pcmz_epp_fit, pcmT_epp_fit
c_chi2_*, c_scale_*, c_overlay_*
```

For Q2-binned overlays it also writes names like:

```text
h_pcmx_epp_SRC_Q2_0, h_pcmx_epp_SRC_Q2_0_fit
g_chi2_pcmx_epp_SRC_Q2_0, g_scale_pcmx_epp_SRC_Q2_0
c_overlay_pcmx_Q2_0
```

`Main_sigmaCM_Hists` should remain useful for comparing against the historical
output. The new `Sigma_CM` directory is now the place to get the same style of
ROOT plotting objects from the current skim-based fit.

## 6. Input Requirements

The skim loader expects a TTree named `srcTree` or `skim`.

Data and MC need:

```text
run, event, Q2, xB, mMiss, kMiss, pMiss, pLead, thetaLead,
leadRegion, hasRecoil, pRec, thetaRec, pcmX, pcmY, pcmZ, pRel, pCM
```

MC also needs:

```text
isMC, genWeight, pcmX_truth, pcmY_truth, pcmZ_truth
```

Metadata must include `sigma_gen` and FD/CD lead-region enum values under one of
these name pairs:

```text
leadRegionFD / leadRegionCD
lead_region_fd / lead_region_cd
fdLeadRegionValue / cdLeadRegionValue
```
