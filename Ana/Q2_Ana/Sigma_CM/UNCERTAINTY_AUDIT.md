# Sigma<sub>CM</sub> Extraction and Uncertainty Audit

## Purpose

This note explains the full `Sigma_CM` extraction chain in plain language and audits how the uncertainty budget is currently constructed. It is written as if you were briefing a collaborator before a talk: the emphasis is on what the analysis is doing physically and statistically, what each study is meant to probe, and where the present implementation hides correlations or analysis choices that matter for interpretation.

The code path discussed here is the current workflow in `Ana/Q2_Ana/Sigma_CM`, with cross-checks against the older cut-variation logic in `Ana/Q2_Ana/Main_Figs_Sys_Err.cpp`.

---

## Executive summary

At a high level, the extraction does the following:

1. Build a compact event sample from data and MC.
2. Define the pair center-of-mass momentum components \(p_{CM,x}\), \(p_{CM,y}\), \(p_{CM,z}\).
3. Reweight each MC event from the generator width to trial \(\sigma_{CM,x}\), \(\sigma_{CM,y}\), \(\sigma_{CM,z}\).
4. Compare data and reweighted MC in the three 1D projections with a simultaneous \(\chi^2\) fit.
5. Take the best-fit widths as the nominal result.
6. Re-run the same extraction under cut, model, fit-range, and closure variations.
7. Turn the observed shifts into an uncertainty budget.

The most important audit conclusion is that the final systematic budget is **not a sum of strictly independent ingredients**. The present machinery is sensible as a practical uncertainty study, but several components share the same data sample, the same MC pool, overlapping fit windows, and in some cases the same random-seed lineage. The code then combines those components in quadrature, which is an analysis choice rather than something guaranteed by the implementation.

Two especially important caveats:

- the budget JSON is built from **integrated** results only, while the plotting code applies those integrated systematic widths as uniform bands on the **Q2-binned** plots;
- the cut-toy contribution used in the total is a **bootstrap-subtracted spread**, which implicitly assumes that cut variation and bootstrap variation are orthogonal.

---

## The full extraction flow

### 1. Build skims from hipo

The first stage converts hipo files into compact ROOT skims containing only the quantities needed downstream.

Conceptually this stage:

- reads each event once;
- reconstructs the relevant kinematics;
- applies detector/data-vs-MC corrections and MC smearing at skim time;
- computes the pair CM observables and truth-level MC quantities;
- stores a reusable event table for repeated fitting.

For each accepted event, the skim keeps:

- inclusive kinematics such as \(Q^2\), \(x_B\), missing mass, missing momentum, and \(k_{miss}\);
- leading and recoil proton quantities;
- \(p_{CM,x}\), \(p_{CM,y}\), \(p_{CM,z}\), \(p_{rel}\), \(p_{CM}\);
- for MC, the generator weight and truth-level \(p_{CM}\) components used for reweighting.

If requested, the MC skim also stores many auxiliary `w_gcf_toy_*` branches. Each of those branches is a different model-weight realization, so the later GCF study can be done without rereading hipo.

### What matters for uncertainties here

This stage already fixes some later correlations:

- every later study reuses the same skimmed event sample unless it explicitly bootstraps or splits it;
- all GCF toy branches are generated from the same underlying MC events;
- any detector smearing or correction choice made here is inherited by every later uncertainty study.

---

### 2. Do the nominal extraction

The nominal extraction is the central result: it finds the \(\sigma_{CM}\) values that best describe the data.

### Physics picture

The analysis compares the measured data distribution to a family of reweighted MC hypotheses:

- \(\sigma_{CM,x}\)
- \(\sigma_{CM,y}\)
- \(\sigma_{CM,z}\)

The MC was generated with one width, \(\sigma_{gen}\), and each event is reweighted to a trial target width using a Gaussian ratio in each truth component. In other words, the code does not regenerate MC for each hypothesis; it transforms the existing MC into a trial hypothesis event-by-event.

### Event selection

Both data and MC are passed through the same cut logic:

- \(x_B\) lower bound
- \(Q^2\) window
- missing-mass window
- \(k_{miss}\) lower bound
- leading-proton momentum threshold
- recoil-proton threshold
- lead-region selection (FD, CD, or both)
- optional \(p_{CM} < p_{rel}\)

For Q2-binned results, the same extraction is repeated independently in each Q2 interval.

### Fit model

The nominal fit is a simultaneous three-projection \(\chi^2\) minimization:

- one histogram in \(p_{CM,x}\)
- one in \(p_{CM,y}\)
- one in \(p_{CM,z}\)

The fit varies:

- \(\sigma_{CM,x}\)
- \(\sigma_{CM,y}\)
- \(\sigma_{CM,z}\)
- one shared normalization scale by default

The data and MC statistical variances are both included in the per-bin \(\chi^2\). The Minuit Hessian then supplies the nominal statistical errors and the scale-width correlations.

### What the nominal statistical error really is

The quoted nominal error bars are fit errors from the Hessian around the best minimum. That means they are:

- local;
- Gaussian/parabolic in spirit;
- tied to the shared-scale simultaneous fit.

They are not the same object as a toy spread or a bootstrap width.

---

### 3. Run cut toys and bootstraps

The cut-toy executable actually performs two different studies in one file.

### Cut toys

For each toy, the code jitters the nominal analysis cuts with Gaussian offsets:

- \(x_B\) and \(Q^2\) thresholds by 0.01;
- missing-mass limits by 0.03;
- \(k_{miss}\) and \(p_{lead}\) thresholds by 0.03;
- FD lead-angle cut by 1 degree;
- recoil threshold by 0.045.

Each toy then reruns the full extraction from scratch.

**Interpretation:**  
This is asking how much the extracted width moves if the analyst had made slightly different but still nearby analysis choices.

### Bootstrap replicas

The same driver also creates bootstrap replicas by giving each data event a Poisson(1) weight and rerunning the same fit.

**Interpretation:**  
This is a data-statistics robustness study. It probes how much the answer fluctuates under sample resampling, not how much it moves under a physics-model change.

### How the budget uses them

The budget script computes:

- the raw cut-toy spread;
- the bootstrap spread;
- a **bootstrap-subtracted cut spread**

using

\[
\sigma_{cuts,used} = \sqrt{\max(\sigma_{cuts,raw}^2 - \sigma_{boot}^2, 0)}.
\]

That subtraction is meant to remove the finite-data-statistics part from the cut-toy spread before using it as a systematic.

### Audit comment

This subtraction is only clean if the cut-toy spread and bootstrap spread can be treated as orthogonal pieces. In practice they are produced from the same data sample and the same extraction machinery, so the subtraction is a pragmatic approximation, not a guaranteed decomposition.

---

### 4. Run GCF toys

The GCF study reruns the extraction once per auxiliary MC weight branch `w_gcf_toy_*`.

Each branch represents a different model-weight realization, all evaluated on the same MC events. The spread across those reruns is taken as the model-weight systematic.

**Interpretation:**  
This answers how much the extracted \(\sigma_{CM}\) depends on the model-reweighting uncertainty encoded in those auxiliary branches.

### Audit comment

This is not independent of everything else:

- it uses the same MC pool as the nominal result;
- it uses the same cut logic as the cut toys;
- it probes a model effect on the same accepted event sample rather than on an independent pseudo-experiment.

So it is a distinct study, but not a statistically independent experiment.

---

### 5. Run combined toys

The combined-toy driver varies two things at once:

- the cut thresholds;
- a randomly chosen GCF auxiliary weight branch, if available.

This is best thought of as a stress test rather than a clean budget component. It measures what happens when several allowed analysis/model shifts move together.

### Important implementation detail

The combined-toy output is **not currently used by the budget assembler**. It is diagnostic only.

---

### 6. Run the fit-range scan

The fit-range scan reruns the nominal fit for a list of X/Y half-windows. Depending on how it is invoked, the Z fit window is either:

- tied to the X/Y choice as \([-w, 2w]\), or
- kept fixed if explicit `--fit-z-min` and `--fit-z-max` were supplied.

This distinction matters for the exact shell script you pasted:

- your first `sigmacm_run_fit_range_scan` command uses only `--xy-ranges`;
- your second command writes to the same output file again but also supplies explicit Z bounds;
- therefore the second command overwrites the first, and the surviving scan is the one with **fixed Z window**.

**Interpretation:**  
This study asks whether the answer depends strongly on exactly which part of the data/MC distributions are allowed to drive the fit.

### How the budget uses it

The budget script computes two related quantities:

- a raw envelope \((\max - \min)/2\);
- a weighted spread across scan points

and uses the **weighted spread** in the total systematic.

### Audit comment

Because all scan points are fit on the same data and mostly overlapping bins, those scan results are highly correlated. The weighted spread is therefore a useful sensitivity diagnostic, but not an independent random uncertainty estimate.

---

### 7. Run closure

The closure study uses MC alone:

1. split the MC into even and odd entries;
2. use the even half as pseudo-data;
3. reweight that pseudo-data to injected target widths;
4. fit it back with the odd half used as the template pool.

The injected targets are 0.10, 0.14, 0.18, and 0.22 GeV/c.

**Interpretation:**  
This asks whether the fitter can recover a known answer, and therefore whether there is an extraction bias.

### How the budget uses it

For each direction, the closure uncertainty is the maximum absolute difference between fitted and injected width across the tested closure points.

### Audit comment

This is a conservative bias measure, but it is again not independent of the full rest of the workflow:

- all closure points reuse the same even/odd split;
- the same fit machinery is being interrogated repeatedly, not rebuilt from independent pseudo-data constructions.

---

### 8. Assemble the budget

The budget script works only on the **integrated** rows of the result trees. It extracts the following pieces:

- nominal statistical uncertainty;
- bootstrap spread;
- bootstrap-subtracted cut-toy spread;
- GCF toy spread;
- fit-range weighted spread;
- closure bias estimate.

It then forms the total systematic as

\[
\sigma_{sys,total} =
\sqrt{
\sigma_{cuts}^2 +
\sigma_{gcf}^2 +
\sigma_{range}^2 +
\sigma_{closure}^2 }.
\]

The script supports excluding selected components, which is exactly what your pasted command does:

- `--exclude-closure`
- `--exclude-fit-range`

So for the workflow you showed, the total systematic in `budget.json` is intentionally **not** the full quadrature sum of every study that was run. It only includes the components you leave enabled at assembly time.

### Audit comment

The budget combination is transparent, but the independence assumption is stronger than what the upstream studies justify. The present total should therefore be described as a **working combined budget under an independence approximation**, not as a covariance-complete uncertainty decomposition.

---

### 9. Make plots

The plotting step reads the ROOT outputs and, if a budget JSON is supplied, adds stat+sys bands.

### Important caveat for Q2 plots

The budget JSON has one row per direction for the **integrated** extraction only. The plotting code then uses those integrated systematic values as if they were the uncertainty for **every Q2 bin** in that direction.

That means the Q2 plots with systematic bands are currently showing:

- Q2-dependent central values and statistical errors;
- but a direction-wise systematic width imported from the integrated budget.

This is convenient for visualization, but it is not the same as a true bin-by-bin systematic propagation.

---

## Where hidden correlations enter

The main hidden or partially hidden correlations are:

### 1. Shared event samples

Almost every study reuses the same skimmed data and MC events. That means many uncertainty components are probing different perturbations of the same underlying sample rather than independent samples.

### 2. Shared MC pool across nominal, cuts, GCF, and fit-range studies

The nominal extraction, cut toys, GCF toys, and fit-range scan all use the same accepted MC event pool. Model uncertainty, cut sensitivity, and fit-range sensitivity are therefore entangled through the same finite-MC fluctuations.

### 3. Shared fit model

Every study goes through the same simultaneous three-axis fit with the same normalization strategy. If there is a structural fit bias or a sensitivity to the shared scale parameter, that coupling appears in multiple budget components at once.

### 4. Shared normalization scale

By default the extraction uses one shared normalization parameter across X, Y, and Z. That couples the fitted widths even before any systematic study is performed.

### 5. Overlapping fit windows

The fit-range scan is highly correlated point-to-point because the windows overlap heavily and are applied to the same events.

### 6. Bootstrap subtraction

The cut-study term used in the total is derived by subtracting bootstrap variance from cut variance in quadrature. That assumes the two can be cleanly separated.

### 7. Reused closure split

All closure injections are tested on the same even/odd MC split, so they are correlated probes of bias rather than independent closure experiments.

### 8. Integrated-only budget applied to Q2 plots

This is not a fit correlation inside the extractor, but it is an important correlation assumption in the final presentation layer: one integrated systematic is reused across all Q2 points in a given direction.

### 9. Seed lineage and reproducibility choices

The current `Sigma_CM` workflow uses deterministic seeded toy generation. That is good for reproducibility, but it also means the various toy ensembles are not independent in a philosophical sense; they are coordinated views of one seeded workflow.

---

## Cross-check against `Main_Figs_Sys_Err`

The old `Main_Figs_Sys_Err.cpp` is not the source of truth for `Sigma_CM`, but it is still useful as a comparison point.

## What is similar

Both workflows vary essentially the same cut families with essentially the same Gaussian widths:

- \(x_B\), \(Q^2\)
- missing-mass bounds
- \(k_{miss}\)
- \(p_{lead}\)
- angle cut
- recoil threshold

So at the level of what is considered a “reasonable cut excursion,” the two approaches are aligned.

## What is different

### 1. Old code works at histogram level, new code works at fit-result level

`Main_Figs_Sys_Err` fills 100 altered histograms and turns their per-bin spread into an uncertainty band.  
`Sigma_CM` reruns the full fit and measures the spread of the extracted widths themselves.

That makes the new workflow more directly tied to the observable of interest.

### 2. Old code mixes cut and reweight variations together

The older code varies the cuts and also uses a randomized reweighter configuration in the same systematic loop. So its systematic band is not a pure cut study.

By contrast:

- `run_cut_toys` isolates cut variations;
- `run_gcf_toys` isolates model-weight variations;
- `run_combined_toys` exists precisely for the “mix everything together” stress test.

This separation is cleaner.

### 3. Old code combines bin spread with nominal stat error directly

In `Main_Figs_Sys_Err`, the final error shown per bin is

\[
\sqrt{\sigma_{sys,bin}^2 + \sigma_{stat,bin}^2}.
\]

That is a straightforward display choice, but it carries no covariance information and is not derived from the full fit model.

### 4. Old code does not carry fit-parameter covariance

The old histogram-band method has no analogue of the current fit-level correlations such as the shared normalization coupling to \(\sigma_X\), \(\sigma_Y\), and \(\sigma_Z\).

### 5. RNG handling is much less transparent in the old code

The old code reinitializes `TRandom3(0)` inside the cut-randomization function each time it is called. Even if ROOT gives distinct seeds in that mode, the behavior is opaque and differs sharply from the deterministic, explicitly seeded toy ensemble used in the new workflow. At minimum, this is a reproducibility discrepancy worth flagging.

---

## Overall assessment of the current uncertainty extraction

## What is strong

- The observable is extracted with a dedicated simultaneous fit rather than inferred indirectly from histogram widths.
- The workflow separates several distinct systematic ideas instead of blending them all together.
- The budget script is simple and auditable.
- Closure, bootstrap, and GCF studies are all present and conceptually sensible.
- The result tree stores enough metadata to reconstruct what each toy actually did.

## What should be said carefully in a presentation

### 1. “Total systematic” is modelled, not covariance-complete

It is a practical combined uncertainty under a quadrature assumption, not the output of a joint covariance analysis.

### 2. The cut systematic is a residual after bootstrap subtraction

That makes it more refined than a raw toy spread, but also more assumption-dependent.

### 3. The Q2 systematic bands are not true Q2-specific systematic extractions

They are integrated systematic numbers painted onto Q2 points.

### 4. Combined toys are a stress test, not a budget ingredient

So there is currently no direct check that the quadrature sum reproduces the spread seen when multiple allowed variations move together.

### 5. Closure is treated as a conservative bias bound

That is reasonable, but it should be described as a bias proxy rather than as a purely statistical uncertainty.

---

## Practical script-specific observations from the pasted workflow

Your pasted shell sequence contains a few interpretation-level details worth noting:

1. `data_skim.root` is assumed to already exist, because the data skim command is commented out.
2. `fit_ranges.root` is produced twice; the second command overwrites the first.
3. Because the second fit-range command supplies explicit `--fit-z-min` and `--fit-z-max`, the surviving scan keeps the Z window fixed while varying only the X/Y half-window.
4. The budget command explicitly excludes closure and fit-range from the combined total, even though those studies were run.
5. The `statOnly` flag is stored in the result trees, but in the present implementation it behaves as metadata rather than as a switch that changes the extractor or budget logic.

So the “full chain” and the “uncertainty components actually included in the final total” are not identical in the example workflow.

---

## Bottom line

The current `Sigma_CM` workflow is a coherent extraction framework with a well-defined nominal fit and a useful suite of uncertainty studies. The main caution is not that the uncertainties are missing, but that several of them are **correlated by construction** while the final budget treats them as if they were independent enough for quadrature.

If this is presented to collaborators, the most accurate summary is:

> We extract \(\sigma_{CM}\) by fitting data to reweighted MC in the three CM momentum projections, then estimate robustness against cut choices, model weights, fit windows, and closure bias with dedicated reruns of the same extractor. The current total systematic is a practical quadrature combination of those studies, but some components share the same samples and fit structure, so the independence assumption should be viewed as approximate rather than exact.
