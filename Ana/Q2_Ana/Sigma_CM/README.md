# Sigma_CM Extraction

This directory is the current place to extract `sigma_CM` and to write the
standard ROOT plotting objects used by `myPlots/plot_sigmaCM.C` and the older
`Main_sigmaCM_Hists` inspection style.

The intended daily workflow is:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_extract 4 sigmaCM.root data.hipo sim.hipo
```

That one command reads hipo data and hipo simulation, runs the nominal
integrated and Q2-binned extraction, and writes one ROOT file.

## What Gets Written

`sigmaCM.root` contains:

```text
sigmaCM                 compact result TTree
profile                 profile TTree, if profile scans were run
sigmacmx_int            integrated sigma_CMx graph
sigmacmy_int            integrated sigma_CMy graph
sigmacmz_int            integrated sigma_CMz graph
sigmacmT_int            integrated transverse sigma_CM graph
sigmacmx_Q2             sigma_CMx vs Q2 graph
sigmacmy_Q2             sigma_CMy vs Q2 graph
sigmacmz_Q2             sigma_CMz vs Q2 graph
sigmacmT_Q2             transverse sigma_CM vs Q2 graph
pcmx_epp                data histogram
pcmx_epp_fit            best-fit MC histogram
g_chi2_*                chi2 scan graphs
g_chi2_joint_*          total X+Y+Z chi2 scans for the simultaneous fit
g_scale_*               scale scan graphs
c_chi2_*                chi2 canvases
c_scale_*               scale canvases
c_overlay_*             data + best-fit MC overlay canvases
```

For Q2 bins, the file also contains objects such as:

```text
h_pcmx_epp_SRC_Q2_0
h_pcmx_epp_SRC_Q2_0_fit
g_chi2_pcmx_epp_SRC_Q2_0
g_scale_pcmx_epp_SRC_Q2_0
c_overlay_pcmx_Q2_0
```

In other words: you do **not** need to run `Main_sigmaCM` and then
`Main_sigmaCM_Hists` just to see the chi2 curves, scale curves, overlays, and
sigma-vs-Q2 graphs. `sigmacm_extract` writes that plotting surface directly.

The `g_chi2_pcmx_*`, `g_chi2_pcmy_*`, and `g_chi2_pcmz_*` graphs are
per-projection diagnostics. The simultaneous-fit diagnostic graphs are named
`g_chi2_joint_sigmaX_*`, `g_chi2_joint_sigmaY_*`, and
`g_chi2_joint_sigmaZ_*`; those use the total X+Y+Z chi2.

## Build

From the repository root on a machine with the full CLAS12/hipo environment:

```bash
cmake -S . -B build
cmake --build build --target sigmacm_extract
```

The executable will be:

```bash
build/Ana/Q2_Ana/Sigma_CM/sigmacm_extract
```

For local development without CLAS12/hipo, this directory can still be built
standalone. In that case hipo input is unavailable, but `--from-skim` still
works:

```bash
cmake -S Ana/Q2_Ana/Sigma_CM -B /tmp/sigmacm_build
cmake --build /tmp/sigmacm_build
```

## Run Nominal Extraction

Hipo input:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_extract 4 sigmaCM.root data.hipo sim.hipo
```

Skim ROOT input, useful for fast reruns or debugging:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_extract --from-skim sigmaCM.root data_skim.root sim_skim.root
```

The skim files must contain a `srcTree` or `skim` tree with the branch contract
listed at the end of this README.

## Common Options

```text
--lead-mode FD|CD|BOTH       default: FD
--independent-scales         fit separate X/Y/Z normalization scales
--pcm-lt-prel                require p_CM < p_rel
--cut-range-xy=0.55          chi2 range for X/Y and transverse projection
--fit-z-min=-0.5             lower chi2 range for Z
--fit-z-max=1.0              upper chi2 range for Z
--xy-ranges=0.45,0.50,0.55  wrapper: fit-window scan values for --full
--couple-z-to-xy             wrapper/scan: set Z window to [-w,2w] for each X/Y scan point
--beam-energy=5.98636        hipo mode beam energy
--max-events=N               quick test run over first N hipo events
--n-gcf-toys=N               wrapper: GCF toy branches in hipo MC cache during --full
--seed=N                     seed used by toy/expert drivers
```

Fit ranges:

```text
X chi2 compares bins from -cutRangeXY to +cutRangeXY
Y chi2 compares bins from -cutRangeXY to +cutRangeXY
Z chi2 compares bins from fitZMin to fitZMax
T chi2 compares bins from 0 to sqrt(2)*cutRangeXY
```

For the nominal fit, set:

```bash
--cut-range-xy=0.55 --fit-z-min=-0.5 --fit-z-max=1.0
```

For the optional fit-range scan, set:

```bash
--xy-ranges=0.40,0.45,0.50,0.55,0.60
```

The scan varies the X/Y half-window and keeps the configured Z window fixed.
Pass `--couple-z-to-xy` only if you intentionally want the older behavior where
each X/Y half-window `w` also sets `fitZMin=-w` and `fitZMax=2w`.

The fit-range scan writes `fit_ranges.root`. Python plotting turns that into
diagnostic PDFs such as:

```text
fit_ranges_sigma_vs_fit_window.pdf
fit_ranges_chi2ndf_vs_fit_window.pdf
fit_ranges_z_window_vs_fit_window.pdf
```

## Plotting

The ROOT file already contains canvases. If you want PDF files from Python:

```bash
Ana/Q2_Ana/Sigma_CM/plot_sigmaCM.py sigmaCM.root --out-dir plots
```

The Python plotter:

- reads ROOT files with `uproot`;
- writes PDF files only;
- does not write PNG files;
- shows stat-only errors unless you explicitly pass a budget JSON.

Python dependencies:

```text
numpy
uproot
matplotlib
```

## Systematics And Sidecar Tables

Systematic/toy studies are intentionally optional because they are slow.

The easiest hipo workflow is the wrapper:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py \
  --from-hipo --A 4 \
  data.hipo sim.hipo sigmacm_out \
  --full
```

That writes:

```text
sigmacm_out_cache/data_skim.root
sigmacm_out_cache/mc_skim.root
sigmacm_out.nominal.root
sigmacm_out.cut_toys.root
sigmacm_out.combined_toys.root
sigmacm_out.fit_ranges.root
sigmacm_out.closure.root
sigmacm_out_plots/*.pdf
```

In hipo mode, the wrapper reads the hipo files once and writes the two cached
skim ROOT files first. Nominal extraction, toys, fit-range scans, closure, and
profile scans then run from those cached ROOT files. This mirrors the old
`Main_sigmaCM` then `Main_sigmaCM_Hists` split: hipo/event reading is stage 1;
ROOT-based fitting/plotting/systematics are stage 2.

To also export `sigmacm_out.budget.json`, `.csv`, and `.tex`, add:

```bash
--export-budget
```

Profile scans are extra diagnostic outputs. Add `--profiles` if you want:

```text
sigmacm_out.profile_axis0.root
sigmacm_out.profile_axis1.root
sigmacm_out.profile_axis2.root
```

For hipo `--full` runs, the MC cache includes `w_gcf_toy_*` auxiliary weight
branches by default, so GCF toys can be included in the budget. Control the
number of branches with `--n-gcf-toys=N`.

If you only want to make the reusable hipo cache:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_make_skim 4 data data_skim.root data.hipo
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_make_skim 4 mc mc_skim.root sim.hipo --gcf-toys=100
```

The systematic/profile drivers operate on skim ROOT files:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_cut_toys data.root mc.root cut_toys.root --n-cut-toys=100 --n-bootstrap=200
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_gcf_toys data.root mc.root gcf_toys.root
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_combined_toys data.root mc.root combined_toys.root --n-toys=100
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_fit_range_scan data.root mc.root fit_ranges.root
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_closure mc.root closure.root
```

Budget sidecars are opt-in. They are not written by default:

```bash
Ana/Q2_Ana/Sigma_CM/budget_assembler.py \
  --nominal sigmaCM.root \
  --cut-toys cut_toys.root \
  --gcf-toys gcf_toys.root \
  --fit-range fit_ranges.root \
  --closure closure.root \
  --out-prefix budget
```

For a hipo budget without GCF toys, omit `--gcf-toys`.

This writes:

```text
budget.json
budget.csv
budget.tex
```

The headline `total_systematic` uses a conservative correlated sum of the
included systematic components. The old quadrature combination is still written
as `total_systematic_quadrature_independence_approx`, and the bootstrap-
subtracted cut spread is still reported separately.

To include stat+sys bands in Python plots, pass the exported JSON:

```bash
Ana/Q2_Ana/Sigma_CM/plot_sigmaCM.py sigmaCM.root --budget-json budget.json --out-dir plots
```

With a budget JSON, integrated plots show stat+sys bands. Q2 plots remain
stat-only unless you also pass `--show-integrated-sys-on-q2`.

## Convenience Wrapper

`run_sigmaCM.py` is a convenience wrapper. It can run either skim ROOT inputs
or hipo inputs.

Quick skim nominal plus plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data_skim.root sim_skim.root sigmacm_out
```

Full skim-based systematics plus plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data_skim.root sim_skim.root sigmacm_out --full
```

Full hipo-based systematics plus plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py --from-hipo --A 4 data.hipo sim.hipo sigmacm_out --full
```

Add profile scans only when needed:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py --from-hipo --A 4 data.hipo sim.hipo sigmacm_out --full --profiles
```

Export budget JSON/CSV/TeX only when wanted:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data_skim.root sim_skim.root sigmacm_out --full --export-budget
```

## Manual Troubleshooting Order

If you run the steps individually, use this order:

```text
1. make skims from hipo
2. nominal extraction
3. optional systematic/diagnostic runners
4. optional budget assembly
5. optional Python plots
```

Starting from hipo:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_make_skim 4 data data_skim.root data.hipo
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_make_skim 4 mc mc_skim.root sim.hipo --gcf-toys=100
```

Then nominal:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_extract --from-skim nominal.root data_skim.root mc_skim.root
```

Then optional systematics:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_cut_toys data_skim.root mc_skim.root cut_toys.root --n-cut-toys=100 --n-bootstrap=200
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_combined_toys data_skim.root mc_skim.root combined_toys.root --n-toys=100
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_fit_range_scan data_skim.root mc_skim.root fit_ranges.root
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_closure mc_skim.root closure.root
```

Optional profiles:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_profile_scan data_skim.root mc_skim.root profile_axis0.root --axis=0
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_profile_scan data_skim.root mc_skim.root profile_axis1.root --axis=1
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_profile_scan data_skim.root mc_skim.root profile_axis2.root --axis=2
```

Run GCF toys when the MC skim has `w_gcf_toy_*` branches. Hipo-made MC skims
get these branches when made with `--gcf-toys=N`, and the wrapper does that
automatically for hipo `--full` runs.

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_gcf_toys data_skim.root mc_skim.root gcf_toys.root
```

Budget without GCF:

```bash
Ana/Q2_Ana/Sigma_CM/budget_assembler.py \
  --nominal nominal.root \
  --cut-toys cut_toys.root \
  --fit-range fit_ranges.root \
  --closure closure.root \
  --out-prefix budget
```

Budget with GCF:

```bash
Ana/Q2_Ana/Sigma_CM/budget_assembler.py \
  --nominal nominal.root \
  --cut-toys cut_toys.root \
  --gcf-toys gcf_toys.root \
  --fit-range fit_ranges.root \
  --closure closure.root \
  --out-prefix budget
```

Plots:

```bash
Ana/Q2_Ana/Sigma_CM/plot_sigmaCM.py nominal.root cut_toys.root combined_toys.root fit_ranges.root closure.root --out-dir plots
```

Plots with stat+sys bands:

```bash
Ana/Q2_Ana/Sigma_CM/plot_sigmaCM.py nominal.root --budget-json budget.json --out-dir plots
```

Dependency map:

```text
make_skim -> everything else if starting from hipo
nominal.root -> budget, main plots
cut_toys.root -> budget
fit_ranges.root -> budget
closure.root -> budget
gcf_toys.root -> budget only if available
combined_toys.root -> plots/diagnostics, not currently used by budget
profile_axis*.root -> profile plots only
```

## Code Map

```text
SigmaCMConfig.*       options and cut configuration
SigmaCMEvent.h        in-memory event record used by the fitter
SigmaCMInput.*        hipo input and optional skim ROOT input
SigmaCMExtractor.*    chi2 model, minimization, profile scan
SigmaCMPlotOutput.*   ROOT histograms, graphs, chi2/scale curves, canvases
SigmaCMResultIO.*     compact sigmaCM/profile TTrees
run_extract.cpp       main daily executable
run_*toys/scans.cpp   expert systematic/profile drivers
make_plots.py         optional PDF plotting from ROOT via uproot
budget_assembler.py   optional JSON/CSV/TeX budget exporter
```

## Skim Branch Contract

This only matters when using `--from-skim` or the expert systematic drivers.

Data and MC need:

```text
run, event, Q2, xB, mMiss, kMiss, pMiss, pLead, thetaLead,
leadRegion, hasRecoil, pRec, thetaRec, pcmX, pcmY, pcmZ, pRel, pCM
```

MC also needs:

```text
isMC, genWeight, pcmX_truth, pcmY_truth, pcmZ_truth
```

Metadata must include `sigma_gen` and FD/CD lead-region enum values under one
of these name pairs:

```text
leadRegionFD / leadRegionCD
lead_region_fd / lead_region_cd
fdLeadRegionValue / cdLeadRegionValue
```
