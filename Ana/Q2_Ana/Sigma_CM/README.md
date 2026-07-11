# Sigma_CM extraction

This is the flat, event-skim based sigma_CM workflow. All source, drivers,
tests, and plotting helpers live in this directory on purpose.

The quantitative basis is unchanged: MC templates are reweighted from
`sigma_gen`, fits minimize the shared-scale chi2
`sum (d - s e)^2 / (sigma_d^2 + s^2 sigma_e^2)`, and template bin errors are
recomputed from `sum w^2` at each parameter point. The default is nominal/stat
only for quick iteration. Use `--full` in the wrapper, or `--full-errors` in the
C++ drivers, when assembling systematic bands.

## Build

From the repository root:

```bash
cmake -S . -B build
cmake --build build --target sigmacm_tests
```

## One-command running

Quick nominal/stat-only extraction and PDF plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out
```

Full pass with cut toys, GCF toys, combined toys, fit-range scan, closure,
profile-chi2 scans, budget tables, and plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out --full
```

Outputs are written next to the chosen prefix:

- `sigmacm_out.nominal.root`: integrated and Q2-binned nominal results
- `sigmacm_out.*.root`: toy, scan, closure, and profile outputs in full mode
- `sigmacm_out.budget.{json,csv,tex}`: systematic budget in full mode
- `sigmacm_out_plots/*.pdf`: sigma-vs-Q2 bands, profile-chi2 curves, and toy
  sigma-hat distributions

## Legacy ROOT compatibility

The old plotting contract is still owned by
`Ana/Q2_Ana/Main_sigmaCM_Hists.cpp`. Its output ROOT object names are unchanged,
including:

- `sigmacmx_int`, `sigmacmy_int`, `sigmacmz_int`, `sigmacmT_int`
- `sigmacmx_Q2`, `sigmacmy_Q2`, `sigmacmz_Q2`, `sigmacmT_Q2`
- `pcmx_epp`, `pcmy_epp`, `pcmz_epp`, `pcmT_epp`
- `pcmx_epp_fit`, `pcmy_epp_fit`, `pcmz_epp_fit`, `pcmT_epp_fit`
- `c_chi2_*`, `c_scale_*`, and `c_overlay_*` canvases

Use that macro when you need a ROOT file that is byte-for-byte shaped for
legacy plotting. Use this directory when you want the cleaner skim-based
nominal/profile/toy/systematic workflow and additional PDF plots.

## Individual drivers

All drivers accept common options from `SigmaCMConfig.cpp`, including
`--seed`, `--lead-mode`, `--q2-bin`, `--aux-weight`,
`--legacy-independent-scales`, `--full-errors`, `--cut-range-xy`,
`--fit-z-min`, and `--fit-z-max`.

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_nominal data.root mc.root nominal.root
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_profile_scan data.root mc.root profile.root --axis=0
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_cut_toys data.root mc.root cut.root --n-cut-toys=100 --n-bootstrap=200
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_gcf_toys data.root mc.root gcf.root --n-toys=100
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_combined_toys data.root mc.root combined.root --n-toys=100
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_fit_range_scan data.root mc.root ranges.root --xy-ranges=0.45,0.50,0.55
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_closure mc.root closure.root
```

Budget and plotting helpers:

```bash
Ana/Q2_Ana/Sigma_CM/budget_assembler.py \
  --nominal nominal.root --cut-toys cut.root --gcf-toys gcf.root \
  --fit-range ranges.root --closure closure.root --out-prefix budget

Ana/Q2_Ana/Sigma_CM/make_plots.py nominal.root cut.root profile.root \
  --budget-json budget.json --out-dir plots
```

## Input requirements

The loader requires fixed skim branches. MC files must also contain `isMC`,
`genWeight`, and `pcmX/Y/Z_truth`. Metadata must contain `sigma_gen` and FD/CD
lead-region enum values under one of these names:
`leadRegionFD`/`leadRegionCD`, `lead_region_fd`/`lead_region_cd`, or
`fdLeadRegionValue`/`cdLeadRegionValue`. Missing schema elements are fatal.
