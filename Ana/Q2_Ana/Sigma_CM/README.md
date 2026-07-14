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
--beam-energy=5.98636        hipo mode beam energy
--max-events=N               quick test run over first N hipo events
--seed=N                     seed used by toy/expert drivers
```

Fit ranges:

```text
X chi2 compares bins from -cutRangeXY to +cutRangeXY
Y chi2 compares bins from -cutRangeXY to +cutRangeXY
Z chi2 compares bins from fitZMin to fitZMax
T chi2 compares bins from 0 to sqrt(2)*cutRangeXY
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
sigmacm_out.profile_axis0.root
sigmacm_out.profile_axis1.root
sigmacm_out.profile_axis2.root
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

The wrapper skips GCF toys for hipo input because those toys require
`w_gcf_toy_*` auxiliary weight branches. Use skim ROOT input if you need that
source included in the budget.

If you only want to make the reusable hipo cache:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_make_skim 4 data data_skim.root data.hipo
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_make_skim 4 mc mc_skim.root sim.hipo
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

To include stat+sys bands in Python plots, pass the exported JSON:

```bash
Ana/Q2_Ana/Sigma_CM/plot_sigmaCM.py sigmaCM.root --budget-json budget.json --out-dir plots
```

## Convenience Wrapper

`run_sigmaCM.py` is a convenience wrapper. It can run either skim ROOT inputs
or hipo inputs.

Quick skim nominal plus plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data_skim.root sim_skim.root sigmacm_out
```

Full skim-based toys/profiles plus plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data_skim.root sim_skim.root sigmacm_out --full
```

Full hipo-based toys/profiles plus plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py --from-hipo --A 4 data.hipo sim.hipo sigmacm_out --full
```

Export budget JSON/CSV/TeX only when wanted:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data_skim.root sim_skim.root sigmacm_out --full --export-budget
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
