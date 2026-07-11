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

Nominal/stat-only extraction, integrated plus Q2-binned:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_run_nominal data.root mc.root nominal.root
```

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
--legacy-independent-scales
--cut-range-xy=0.55
--fit-z-min=-0.5
--fit-z-max=1.0
```

## 3. Python Helpers

The Python scripts are optional. They are only for assembling systematic budget
tables and making PDF/PNG plots from the ROOT files produced by the C++
executables.

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
Ana/Q2_Ana/Sigma_CM/budget_assembler.py \
  --nominal nominal.root \
  --cut-toys cut_toys.root \
  --gcf-toys gcf_toys.root \
  --fit-range fit_ranges.root \
  --closure closure.root \
  --out-prefix budget

Ana/Q2_Ana/Sigma_CM/make_plots.py \
  nominal.root cut_toys.root gcf_toys.root combined_toys.root profile_x.root profile_y.root profile_z.root \
  --budget-json budget.json \
  --out-dir plots
```

This writes:

```text
budget.json
budget.csv
budget.tex
plots/*.pdf
plots/*.png
```

## 4. One-Command Wrapper

`run_sigmaCM.py` is just a convenience wrapper around the C++ executables plus
the Python helpers.

Use it only after the Python environment above is working.

Quick nominal run plus plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out
```

Full run with toys, profiles, budget, and plots:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out --full
```

If you built Sigma_CM directly instead of through the full repo, pass that build
directory:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data.root mc.root sigmacm_out --build-dir /tmp/sigmacm_build
```

## 5. Legacy ROOT Output

The old ROOT-object contract is still produced by
`Ana/Q2_Ana/Main_sigmaCM_Hists.cpp`, not by the new skim-based executables.

Use `Main_sigmaCM_Hists` when legacy plotting needs the exact old object names:

```text
sigmacmx_int, sigmacmy_int, sigmacmz_int, sigmacmT_int
sigmacmx_Q2, sigmacmy_Q2, sigmacmz_Q2, sigmacmT_Q2
pcmx_epp, pcmy_epp, pcmz_epp, pcmT_epp
pcmx_epp_fit, pcmy_epp_fit, pcmz_epp_fit, pcmT_epp_fit
c_chi2_*, c_scale_*, c_overlay_*
```

Use this `Sigma_CM` directory when you want the newer skim-based nominal,
profile, toy, systematic-budget, and plotting workflow.

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
