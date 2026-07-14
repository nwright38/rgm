# Sigma_CM workflow refactor plan

## Goal

Make `Ana/Q2_Ana/Sigma_CM` feel like the old `Main_sigmaCM` +
`Main_sigmaCM_Hists` workflow again:

1. one daily command that reads data/sim inputs and writes the familiar ROOT
   plotting surface;
2. no fragile requirement that precomputed skim branches already exist;
3. profile chi2, convergence, fit errors, and unfolding diagnostics written as
   first-class output objects;
4. systematic uncertainties available behind an explicit slow flag;
5. quantitative unfolding machinery kept close to the nominal extraction rather
   than bolted on through scattered helper scripts.

## Current shape

The current code already has a useful core:

- `SigmaCMExtractor` performs the multi-parameter chi2 minimization.
- `run_nominal` runs integrated plus Q2-binned fits.
- `SigmaCMLegacyOutput` writes the standard ROOT plotting objects expected by
  `myPlots/plot_sigmaCM.C` and the old inspection macros.
- `run_sigmaCM.py` can orchestrate nominal, toys, profiles, budgets, and plots.

The pain points are mostly at the workflow boundary:

- `SigmaCMSkimIO` requires a fixed `srcTree`/`skim` branch contract.
- Direct hipo-to-result running is not the main path.
- Nominal, profile scans, toys, budgets, and plotting are split across many
  executables/scripts.
- Profile-chi2 information exists, but is not produced by the nominal command.
- The output is split between compact TTrees, plotting graphs, canvases, JSON,
  CSV, TeX, and Python plots, so the "where do I look?" answer is not obvious.
  The refactor should make ROOT the default analysis artifact and make sidecar
  table files an explicit export step.

## Target daily workflow

Keep the old mental model:

```bash
sigmacm_extract A out.root data.hipo sim.hipo
root -l -q 'myPlots/plot_sigmaCM.C("out.root","sigmaCM.pdf")'
```

Optional slow systematics:

```bash
sigmacm_extract A out.root data.hipo sim.hipo --systematics
```

Optional explicit skim mode remains available for fast iteration:

```bash
sigmacm_extract --from-skim out.root data_skim.root sim_skim.root
```

The output ROOT file should be enough for normal inspection. Python helpers can
read the ROOT output with uproot. JSON/CSV/TeX files should be optional exports,
not files produced by default in every run.

## Proposed architecture

### 1. Input adapters

Introduce a small input-adapter layer:

- `SigmaCMHipoIO`: reads hipo files with `HipoChain`/`clas12ana`, computes
  the kinematics currently expected as branches, and fills `sigmacm::Event`.
- `SigmaCMSkimIO`: keeps the current ROOT tree loader for debugging and
  fast reruns.

This lets the extractor stay independent of whether events came from hipo or a
ROOT skim. It also removes the feeling that branch production is a prerequisite
for the analysis.

### 2. One daily driver

Add `sigmacm_extract` as the main executable:

- default: hipo data + hipo sim to one ROOT file;
- `--from-skim`: use existing skim inputs;
- always run integrated and Q2-binned nominal fits;
- always write standard ROOT plotting graphs/hists/canvases;
- always write compact result/profile TTrees;
- optionally write a small text/JSON config note into the ROOT file.
- write JSON/CSV/TeX sidecars only with an explicit export flag.

Existing executables can remain during migration, but should become expert
tools or implementation details.

### 3. Profile and convergence output by default

The nominal output should include:

- `sigmaCM`: one row per integrated/Q2-bin fit;
- `profile`: chi2 scan points for x/y/z, integrated and Q2-binned where
  requested;
- `fitDiagnostics`: converged flag, status, chi2, ndf, per-axis chi2, scale,
  scale error, effective MC entries, and covariance/correlation summaries;
- standard `g_chi2_*` graphs with visible best-fit and error markers.

This directly supports quick questions like:

- Did Minuit converge?
- Where are the +-1 chi2 crossings?
- Is one Q2 bin failing?
- Is one projection driving the chi2?

### 4. Systematics as an explicit mode

Keep systematics slow and opt-in:

- `--systematics cut,gcf,fit-range,closure`
- `--systematics all`
- `--n-toys`, `--n-bootstrap`, `--seed`

The slow mode should append to the same output file when possible:

- `sigmaCM_systematics`
- `sigmaCM_budget`
- `toy_sigmaCM`
- `closure`

The budget JSON/CSV/TeX outputs are still useful, but should be opt-in exports.
The ROOT file should contain enough budget information for uproot-based plotting
so routine runs do not clog the working directory with sidecars.

### 5. Better unfolding framework

Treat unfolding as a named extraction method rather than a separate universe:

- `--method chi2-template` for the current multi-parameter template fit;
- `--method width-match` for the current width-matching cross-check;
- `--method response-unfold` for a future response-matrix method.

For response unfolding, define the output contract early:

- response matrix per axis and Q2 bin;
- regularization choice and strength;
- closure/pull diagnostics;
- covariance matrix for unfolded sigma values;
- same final `sigmacm*_int` and `sigmacm*_Q2` graphs.

That way plotting code does not care which quantitative method produced the
numbers.

## Migration steps

1. Done: keep current `loadSkim` in `SigmaCMSkimIO`.
2. Done: add `SigmaCMHipoIO` to compute hipo event kinematics directly into
   `sigmacm::Event`.
3. Done: add `sigmacm_extract` and make it reproduce current `run_nominal`
   output in `--from-skim` mode.
4. Done: make `sigmacm_extract A out.root data.hipo sim.hipo` the default
   hipo input path in full CLAS12/hipo builds.
5. Fold profile scans into the nominal output behind a cheap default, e.g.
   integrated profiles always, Q2 profiles with `--profiles q2`.
6. Add ROOT-stored systematics budget and make `--systematics` append those
   outputs. Add `--export-budget-tables` or similar for JSON/CSV/TeX sidecars.
7. Update `myPlots/plot_sigmaCM.C` or add a companion macro that draws:
   sigma vs Q2, integrated bands, profile chi2 with +-1 crossings, convergence
   summaries, and fit/data overlays.
8. Deprecate the pile of individual `sigmacm_run_*` commands in the README once
   the daily driver covers the real workflow.

## First implementation checkpoint

The first checkpoint should be intentionally modest:

```bash
sigmacm_extract --from-skim out.root data.root sim.root --profiles integrated
```

It should produce:

- the same `sigmaCM` tree as `run_nominal`;
- the same standard ROOT plotting objects as `run_nominal`;
- profile chi2 curves in the same file;
- a clear console summary table of integrated and Q2-binned sigma values.

Once that is stable, direct hipo input becomes a contained adapter problem
instead of a rewrite of the fitter.
