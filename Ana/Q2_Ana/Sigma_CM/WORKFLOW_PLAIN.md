# Sigma_CM Plain-Language Workflow

This file explains what this directory is trying to do, without assuming you
want to remember every executable name or implementation detail.

## The Big Picture

The goal is to extract the center-of-mass momentum widths:

```text
sigma_CM,x
sigma_CM,y
sigma_CM,z
```

The extraction compares data to simulation. The simulation is reweighted to
different trial `sigma_CM` values. For each trial, the code asks:

```text
How well does this reweighted MC describe the data?
```

The preferred `sigma_CM` value is the one with the best chi2.

This directory is meant to reproduce the useful parts of the older workflow:

```text
Main_sigmaCM       reads hipo and makes ROOT objects
Main_sigmaCM_Hists reads ROOT objects and makes chi2/sigma plots
```

The newer workflow keeps that same idea:

```text
hipo files
  -> optional cached skim ROOT files
  -> nominal sigma_CM extraction
  -> optional systematic studies
  -> optional PDF plots / JSON budget
```

## The Daily Nominal Workflow

For the usual stat-only extraction from hipo:

```bash
./build/Ana/Q2_Ana/Sigma_CM/sigmacm_extract 4 sigmaCM.root data.hipo sim.hipo
```

This reads the data and MC hipo files, applies the analysis cuts, fits the
integrated `sigma_CM`, fits each Q2 bin, and writes one ROOT file.

That ROOT file is the main output. It contains:

```text
sigmaCM tree              fitted sigma values, errors, chi2, settings
sigma-vs-Q2 graphs        sigma_CM in each Q2 bin
chi2 scan graphs          chi2 as a function of trial sigma
joint chi2 scan graphs    total X+Y+Z chi2 as a function of trial sigma
scale scan graphs         best MC normalization scale vs trial sigma
overlay histograms        data and best-fit MC projections
ROOT canvases             quick ROOT-viewable diagnostic canvases
```

In plain terms:

```text
sigmacm_extract answers:
"What sigma_CM best describes my data, and what do the chi2 curves look like?"
```

There are two kinds of chi2 graphs:

```text
g_chi2_pcmx_*, g_chi2_pcmy_*, g_chi2_pcmz_*
  per-projection diagnostic chi2 curves

g_chi2_joint_sigmaX_*, g_chi2_joint_sigmaY_*, g_chi2_joint_sigmaZ_*
  total X+Y+Z chi2 curves connected to the simultaneous fit
```

If you are diagnosing the simultaneous fit, start with the `joint` graphs.
If you are asking which projection is causing tension, then look at the
projection-specific graphs.

## Why There Is A Skim Cache

For slow studies, rereading hipo over and over is painful. So the wrapper can
make a compact ROOT cache first:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py \
  --from-hipo --A 4 \
  data.hipo sim.hipo sigmacm_out \
  --full
```

In hipo mode, the wrapper first writes:

```text
sigmacm_out_cache/data_skim.root
sigmacm_out_cache/mc_skim.root
```

Those cached files contain the event quantities needed by the fitter:

```text
Q2, xB, missing mass, missing momentum, lead/recoil kinematics,
pcmX, pcmY, pcmZ, pRel, pCM, MC truth pcm values, MC weights
```

Then the slow steps run from the cache instead of rereading hipo. This is the
new version of the old two-stage idea:

```text
stage 1: hipo -> reusable ROOT cache
stage 2: ROOT cache -> fits, toys, plots, budgets
```

If you already have skim ROOT files, you can skip hipo entirely:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py data_skim.root sim_skim.root sigmacm_out --full
```

## What The Main Files Mean

### `sigmacm_extract`

This is the main nominal extractor.

It can run directly from hipo:

```bash
sigmacm_extract 4 sigmaCM.root data.hipo sim.hipo
```

or from cached skim ROOT files:

```bash
sigmacm_extract --from-skim sigmaCM.root data_skim.root sim_skim.root
```

It writes the main result ROOT file.

Use this when you want the central result, Q2-binned result, chi2 curves, and
ROOT plotting objects.

### `sigmacm_make_skim`

This only makes the cache.

```bash
sigmacm_make_skim 4 data data_skim.root data.hipo
sigmacm_make_skim 4 mc mc_skim.root sim.hipo --gcf-toys=100
```

It does not fit `sigma_CM`. It just converts hipo into the compact ROOT event
format used by the rest of the directory.

Use this when you want to read hipo once, then rerun fits/systematics quickly.

### `run_sigmaCM.py`

This is the convenience wrapper.

Nominal only:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py --from-hipo --A 4 data.hipo sim.hipo sigmacm_out
```

Nominal plus slow systematic studies:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py --from-hipo --A 4 data.hipo sim.hipo sigmacm_out --full
```

Nominal plus systematics plus JSON/CSV/TeX budget sidecars:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py --from-hipo --A 4 data.hipo sim.hipo sigmacm_out --full --export-budget
```

Add profile scans only when you really want that diagnostic:

```bash
Ana/Q2_Ana/Sigma_CM/run_sigmaCM.py --from-hipo --A 4 data.hipo sim.hipo sigmacm_out --full --profiles
```

In plain terms:

```text
run_sigmaCM.py is the "please do the usual sequence for me" script.
```

## What The Systematic Studies Mean

Systematic studies are optional because they take longer and produce extra
files. They all follow the same basic idea:

```text
rerun the extraction under a controlled variation
see how much sigma_CM moves
```

That movement is used as an uncertainty estimate or diagnostic.

### Cut Toys: `sigmacm_run_cut_toys`

Cut toys randomly vary reasonable analysis-cut choices, then rerun the fit.

Examples of varied choices:

```text
xB lower cut
Q2 lower cut
missing-mass window
kMiss cut
pLead cut
theta cut
pRec cut
```

Physics/statistics meaning:

```text
If I had made slightly different but still reasonable analysis choices,
how much would sigma_CM change?
```

This is a systematic-uncertainty proxy for cut sensitivity.

### Bootstrapping: also in `sigmacm_run_cut_toys`

The same executable also does bootstrap replicas.

Bootstrapping gives each data event a random Poisson weight. Some events count
zero times, some once, some more than once. Then the fit is rerun.

Physics/statistics meaning:

```text
If I had another statistically similar data sample,
how much might sigma_CM fluctuate?
```

This is mostly a statistical robustness check. Be careful not to blindly add it
on top of the nominal fit error, because both can describe finite-statistics
uncertainty.

The clean interpretation is:

```text
nominal fit error       official statistical error
bootstrap spread        check of statistical stability
cut-toy spread          systematic sensitivity to cut choices
```

### Combined Toys: `sigmacm_run_combined_toys`

Combined toys vary several things together. They can include cut variations
and, when available in the skim, auxiliary MC weight variations.

Physics/statistics meaning:

```text
If several allowed analysis/model choices move at the same time,
what spread in sigma_CM do I get?
```

This is useful as a stress test. It is not always the cleanest final budget
component because it can mix several sources at once.

### GCF Toys: `sigmacm_run_gcf_toys`

GCF toys use MC auxiliary weight branches named like:

```text
w_gcf_toy_*
```

Each branch represents a different model-weight variation. The fit is rerun
for each one.

Physics/statistics meaning:

```text
How much does sigma_CM depend on this model-weight uncertainty?
```

Important practical note:

```text
Hipo-made MC caches contain w_gcf_toy_* branches when made with --gcf-toys=N.
```

The wrapper adds those branches automatically for hipo `--full` runs. If you
make the MC cache by hand and want GCF toys, remember to pass `--gcf-toys=N`.

### Fit-Range Scan: `sigmacm_run_fit_range_scan`

This changes the histogram range included in the chi2 comparison.

Physics/statistics meaning:

```text
Does sigma_CM depend strongly on exactly which part of the distribution
I include in the fit?
```

If it moves a lot, that can mean the model describes the core and tails
differently, or that the chosen fit range matters more than expected.

By default the scan tries X/Y half-windows:

```text
0.45, 0.50, 0.55 GeV/c
```

By default, the Z window stays at the configured `fitZMin`/`fitZMax` while the
X/Y half-window changes. If `--couple-z-to-xy` is passed, each X/Y half-window
`w` also sets the Z window to:

```text
[-w, 2w]
```

The useful diagnostic PDFs are:

```text
fit_ranges_sigma_vs_fit_window.pdf
fit_ranges_chi2ndf_vs_fit_window.pdf
fit_ranges_z_window_vs_fit_window.pdf
```

The first plot answers:

```text
How much does sigma_CM move as I change the fit window?
```

The second answers:

```text
Does chi2/ndf improve or worsen as I change the fit window?
```

### Closure: `sigmacm_run_closure`

Closure uses MC to make pseudo-data with known injected widths, then fits it
back.

Physics/statistics meaning:

```text
If I know the true sigma_CM, does this fitter recover it?
```

Closure is a bias check. If the fitter systematically misses the injected
value, that bias should be understood before trusting the extraction too much.

### Profile Scans: `sigmacm_run_profile_scan`

Profile scans fix one sigma value at a time, refit the other parameters, and
record the chi2.

Physics/statistics meaning:

```text
What does the chi2 curve look like near the best fit when I force sigma_CM
to different values?
```

This is a diagnostic for the uncertainty shape and fit stability. It is useful,
but not needed for every routine run. That is why `run_sigmaCM.py` only runs it
when you pass:

```bash
--profiles
```

## What The Budget Means

`budget_assembler.py` reads the ROOT outputs from the nominal and systematic
studies and makes optional sidecar tables:

```text
budget.json
budget.csv
budget.tex
```

The budget is meant to collect uncertainty components like:

```text
statistical
cut-toy spread
GCF-toy spread, if available
fit-range envelope
closure bias estimate
total systematic
```

The headline `total_systematic` is a conservative fully correlated sum:

```text
total systematic = cut + gcf + fit_range + closure
```

The JSON also keeps `total_systematic_quadrature_independence_approx` for
comparison. By default the cut term uses the raw cut-toy spread; the older
bootstrap-subtracted cut term is still reported and can be selected explicitly
with `--use-stat-subtracted-cuts`.

## What The Python Plots Do

`make_plots.py` and `plot_sigmaCM.py` read ROOT output with `uproot` and write
PDFs. They do not write PNGs.

They make things like:

```text
integrated sigma_CM summary
sigma_CM vs Q2
toy distributions, only when enough toy entries exist
profile chi2 plots, only when profile output exists
```

If you pass a budget JSON, the integrated sigma plot shows stat+systematic
bands. Q2 plots stay stat-only by default because the budget is integrated-only.
Use `--show-integrated-sys-on-q2` only when you explicitly want the old
presentation choice of painting the integrated systematic onto every Q2 point.

## Which Output Should I Look At First?

For a normal run, start with:

```text
sigmacm_out.nominal.root
sigmacm_out_plots/*sigma_vs_q2_stat.pdf
sigmacm_out_plots/*integrated_sigma_stat.pdf
```

Inside the ROOT file, useful objects include:

```text
sigmaCM                 compact numeric result tree
sigmacmx_Q2             sigma_CMx vs Q2
sigmacmy_Q2             sigma_CMy vs Q2
sigmacmz_Q2             sigma_CMz vs Q2
g_chi2_*                chi2 curves
g_chi2_joint_*          total X+Y+Z chi2 curves
c_overlay_*             data/MC overlays
```

For a systematic run with budget export, also inspect:

```text
sigmacm_out.budget.json
sigmacm_out.budget.csv
sigmacm_out_plots/*stat_and_total.pdf
```

## Mental Model

The workflow is easiest to remember this way:

```text
1. Convert or read events.
2. Build data and MC distributions.
3. Reweight MC to trial sigma_CM values.
4. Compute chi2(data, reweighted MC).
5. Pick the sigma_CM with the best chi2.
6. Repeat under variations to estimate systematics.
7. Plot the central values and uncertainties.
```

The nominal extraction answers:

```text
What is sigma_CM?
```

The systematic studies answer:

```text
How much would sigma_CM move if reasonable analysis/model choices changed?
```

The bootstrap answers:

```text
Does the fitted value look stable under finite-data-sample fluctuations?
```

The closure test answers:

```text
Can the fitter recover a known injected answer?
```
