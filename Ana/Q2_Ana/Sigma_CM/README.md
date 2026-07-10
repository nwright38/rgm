# Sigma_CM extraction

This directory contains the C++17 + ROOT pipeline for the SRC `(e,e'pp)`
sigma_CM extraction from flat ROOT skims. It is built as part of the repository
through `Ana/Q2_Ana/CMakeLists.txt`, while keeping the code in its own directory
so this analysis can evolve separately from the older Q2 programs. It
intentionally does not depend on clas12root, HIPO, or experiment framework
libraries.

## Build

From the repository root, using the normal repository build:

```bash
cmake -S . -B build
cmake --build build --target sigmacm_tests
```

For quick local iteration on only this package, the subdirectory can also be
configured directly:

```bash
cmake -S Ana/Q2_Ana/Sigma_CM -B /tmp/sigmacm_build
cmake --build /tmp/sigmacm_build --target sigmacm_tests
```

Driver executables are named `sigmacm_run_nominal`, `sigmacm_run_cut_toys`,
`sigmacm_run_gcf_toys`, `sigmacm_run_combined_toys`, `sigmacm_run_closure`,
`sigmacm_run_fit_range_scan`, `sigmacm_run_diagnostics`, and
`sigmacm_run_profile_scan`.

## Input requirements

The loader requires the fixed skim branches from the task spec. MC files must
also contain `isMC`, `genWeight`, and `pcmX/Y/Z_truth`. Metadata must contain
`sigma_gen` and FD/CD lead-region enum values under one of these names:
`leadRegionFD`/`leadRegionCD`, `lead_region_fd`/`lead_region_cd`, or
`fdLeadRegionValue`/`cdLeadRegionValue`. Missing schema elements are fatal.

## Examples

Nominal integrated plus Q2-binned extraction:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_nominal data.root mc.root nominal.root --seed 17
```

Cut toys plus data bootstrap:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_cut_toys data.root mc.root cut_toys.root --seed 22 --n-cut-toys=100 --n-bootstrap=200
```

GCF model-variation toys:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_gcf_toys data.root mc.root gcf_toys.root --seed 23
```

Closure:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_closure mc.root closure.root --seed 24
```

Fit-range scan:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_fit_range_scan data.root mc.root ranges.root --xy-ranges=0.45,0.50,0.55,0.60
```

Profile scan:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_profile_scan data.root mc.root profile_x.root --axis=0 --scan-min=0.08 --scan-max=0.26 --n-points=61
```

Budget assembly:

```bash
python3 Ana/Q2_Ana/Sigma_CM/scripts/budget_assembler.py \
  --nominal nominal.root --cut-toys cut_toys.root --gcf-toys gcf_toys.root \
  --fit-range ranges.root --closure closure.root --out-prefix budget
```

Plots:

```bash
python3 Ana/Q2_Ana/Sigma_CM/scripts/make_plots.py nominal.root cut_toys.root closure.root --out-dir plots
```

## Notes

The fit minimizes the shared-scale chi2
`sum (d - s e)^2 / (sigma_d^2 + s^2 sigma_e^2)` with Minuit2. Template bin
errors are recomputed from `sum w^2` at every parameter point. The legacy
independent-scale mode is parsed for validation compatibility, but the default
and recommended extraction is the shared-scale fit.
