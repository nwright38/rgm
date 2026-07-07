# Ana/Q2_Ana — Script Reference

This directory contains the C++ analysis programs and helper libraries that process raw CLAS12
HIPO event files and produce the ROOT histograms and graphs consumed by the `myPlots/` plotting
scripts.

## Overall Pipeline

```
HIPO data/sim files
        ↓
[Andrew_skim]         ← optional: skims to single-electron + lead-SRC-proton events
        ↓
Main_Figs_Binned      ← main analysis: kinematics, cuts, histogram filling (nominal + 100 toys)
        ↓
output.root  (diffTable, integratedTable, hists/nominal, hists/toy_NNN)
        ↓
BuildGraphs           ← converts flat tables into TGraphAsymmErrors objects in graphs/ directory
        ↓
ExtractFitQuantities.C  ← adds fit-derived graphs (σ_CM, E_miss moments) to graphs/
        ↓
Python plotting scripts (myPlots/)
```

### Q² Reweighting (optional):
```
Q2_Reweight  →  weights.root
Main_Figs_Binned --q2-reweight weights.root  →  Q²-reweighted histograms
```

### σ_CM extraction pipeline:
```
Main_sigmaCM          ← stage 1: fills template bank of σ_CM scan histograms
Main_sigmaCM_WidthMatch.C  ← stage 2: width-matching extraction
plot_sigmaCM.C        ← final plots
```

---

## Standard Kinematic Cuts

Most programs in this directory use the same nominal SRC selection:

| Variable | Cut |
|---|---|
| x_B | > 1.2 |
| Q² | 1.5 – 5.0 GeV² |
| p_lead | > 1.0 GeV |
| θ_lead | < 37° (FD only) |
| k_miss (ZQ) | > 0.3 GeV |
| m_miss | 0.65 – 1.10 GeV |
| p_recoil | > 0.30–0.35 GeV (CD only) |

Beam energy: **5.98636 GeV** (or 5.984792 GeV in some calibration programs).

---

## Headers

### `Q2Reweight.h`

Minimal helper class that loads a 1D ROOT histogram of Q² weights and provides a per-event
weight lookup. Used by `Main_Figs_Binned.cpp` when `--q2-reweight` is passed.

**Class: `Q2Reweight`**

| Method | Description |
|---|---|
| `load(fileName)` | Opens a ROOT file, clones the histogram `q2_reweighting/h_Q2_reweight` (with fallback to `h_Q2_reweight` at the top level), detaches it from the file. Returns `bool` — `false` if the file or histogram is missing. |
| `enabled()` | Returns `true` if a valid weight histogram is loaded. |
| `weight(q2)` | Returns the bin content for the given Q² value. Returns `1.0` if not enabled, `0.0` if Q² is out of range. |

---

### `many_plots.h` / `many_plots.cpp`

Container for matched (e,e'p) and (e,e'pp) histogram sets including per-Q² sub-histograms. Used
by the older **analysis-note** programs (`Analysis_Note_Q2.cpp`, etc.).

**Class: `many_plots`**

Owns: `h_ep`, `h_epp`, `h_ep_Q2bin[10]`, `h_epp_Q2bin[10]`.

Q² bin edges (hard-coded, 10 bins):
`{1.50, 1.65, 1.80, 1.95, 2.10, 2.25, 2.40, 2.70, 3.00, 3.50, 5.00}`

| Method | Description |
|---|---|
| `many_plots(name, title, xmin, xmax)` | Books fresh histograms: `ep` with 100 bins, per-Q² ep/epp with 100/25 bins |
| `many_plots(name, rootFile)` | Retrieves existing histograms from a ROOT file by base name |
| `Fill_hist_set(is_epp, Q2, x, wep, wepp)` | Fills the ep histogram always; fills epp if `is_epp`. Both use the appropriate weight. Also fills the correct Q² bin histogram. |
| `Write_hist_set()` | Writes the ep histogram (and Q² bin ep histograms) to the current ROOT directory |
| `Write_hist_set_epp()` | Writes the epp histograms |
| `Write_ratio_set()` | Clones epp, divides by ep bin-by-bin, writes ratio histograms, draws them on a canvas |
| `Draw_*` helpers | Various methods to overlay or draw individual histogram sets |
| `getScale_ep()` / `getScale_epp()` | Return the ep/epp histogram integrals |
| `binQ2(q2)` | Maps a Q² value to the bin index `{1.50, 1.65, ..., 5.00}` |

---

### `BinnedHistStore.h`

**Generic N-dimensional histogram booking and filling engine** used by `Main_Figs_Binned.cpp`.
Separates the physics event calculation from histogram bookkeeping. Understands how to produce
both the differential and integrated flat tables that `BuildGraphs.cpp` then converts to graphs.

#### Supporting types

**`Axis`** — describes one selector axis (e.g. pMiss, Q²):
- `name` — string identifier
- `edges` — vector of bin edges
- `nbins()`, `findBin(x)` — bin count and bin lookup; returns -1 for out-of-range
- `center(i)` — bin centre

**`Selection`** — enum: `EP` or `EPP`

**`FillTask<EventT>`** — describes one family of histograms to fill:

| Field | Description |
|---|---|
| `name` | Task name (used in graph naming) |
| `selection` | Which events to fill (`EP` or `EPP`) |
| `selectorAxes` | List of `Axis` objects that define the "selector" dimensions (e.g. pMiss bin, Q² bin) |
| `selectorValues` | Functions `EventT → double` to evaluate each selector axis per event |
| `valueAxis` | The axis whose variable is actually histogrammed (x-axis of the output graph) |
| `valueFunc` | Function `EventT → double` to evaluate the histogrammed variable |
| `fillWeight` | Function `EventT → double` for the fill weight (often just `wep` or `wepp`) |

**`DiffRow` / `IntegratedRow`** — output structs for the flat tables written to `TTree`s.

#### `HistStore<EventT>`

Owns all booked histograms for one store label (e.g. `"nominal"`, `"toy_007"`).

| Method | Description |
|---|---|
| `HistStore(label, tasks)` | Books all histograms for the given list of `FillTask`s |
| `fill(event, mode)` | For each task, evaluates selector axes and fills the correct histogram(s) |
| `write(outDir)` | Writes all histograms into a subdirectory named `<label>` under `outDir` |
| `rawHists(taskName, selection)` | Returns the map of selector-bin → `TH1D*` for a given task |

#### Free functions

| Function | Description |
|---|---|
| `buildDiffRows(nominalStore, toyStores, task, sel, cvm)` | Builds a vector of `DiffRow`s. For each selector bin + value bin, the nominal histogram gives the central value; toy spread gives systematic uncertainty. |
| `buildIntegratedRows(nominalStore, toyStores, task, sel, collapseAxis, cvm)` | Same but collapses one or all selector axes — used for the `integratedTable`. |
| `buildRatioDiffRows(nominalStore, toyStores, numTask, denTask, cvm)` | Builds ratio `DiffRow`s by dividing cloned histograms bin-by-bin; toy percentiles (p16/p84) provide asymmetric systematics. |

---

## Analysis Programs

### `Andrew_skim.cpp`

Event skimmer — keeps only events with exactly one electron passing standard cuts and exactly one
SRC lead proton.

**Usage:**
```
./Andrew_skim outputfile.hipo inputfiles.hipo...
```

**What it does:** Reads HIPO via `HipoChain`, applies `clas12ana` electron cuts plus
`xB > 1`, `Q² > 1`, calls `getLeadRecoilSRC`, and requires exactly one lead proton candidate.
Passes the full event to a `HipoChainWriter`.

**Note:** Uses deuteron target mass for SRC reconstruction; QADB is disabled.

---

### `iso_p_LUND.cpp`

Toy generator that writes uncorrelated electron + proton two-body events to a LUND text file and
a ROOT tree. Useful for detector acceptance studies.

**Usage:**
```
./iso_p_LUND nEvents output.lund output.root
```

**Output ROOT tree** (`genT`): branches `pe[3]` (electron 3-momentum), `pp[3]` (proton 3-momentum).

**Key functions:**

| Function | Description |
|---|---|
| `randomVertex(target)` | Samples beam-spot (x, y) and z from a target map (`4-foil`, `1-foil`, `Ar`, `Ca`, `liquid`). Global z offset = −3 cm. |
| `addParticle(...)` | Formats one LUND particle line |
| `generate_event(ve, vp)` | Draws electron momentum 1–6 GeV forward-biased; proton 0.2–3.5 GeV isotropically |

---

### `p_LUND.cpp`

Physics-motivated elastic electron-proton generator with Rosenbluth cross-section weighting and
soft-photon radiative corrections.

**Usage:**
```
./p_LUND nEvents output.lund output.root
```

**Output ROOT tree**: branches `pe[3]`, `pp[3]`, invariant mass `W`, and event `weight`.
Beam energy: 5.98636 GeV.

**Key functions:**

| Function | Description |
|---|---|
| `generate_event(weight, ve, vp)` | Samples Q², computes elastic kinematics, applies radiative corrections, evaluates Rosenbluth weight. Sets weight = 0 for kinematically invalid events. |
| `sigma_onShell_by_Etheta` | Evaluates Rosenbluth cross section using Kelly form factors |
| `radiateElectron` | Applies soft-photon energy loss |
| `radiationFactor(Q2)` | Multiplicative radiative correction from `deltaHard(Q²)` |
| `deltaHard(Q2)` | Hard-radiative term |

---

### `ep_Kinematics.cpp`

Basic (e,e'p) kinematics diagnostic — fills and saves distributions of standard SRC observables.

**Usage:**
```
./ep_Kinematics output.root output.pdf inputfiles.hipo...
```

**Cuts:** Standard `clas12ana` electron cuts (sampling fraction, PCAL, edge, PID, vertex) +
exactly one SRC lead proton.

**Histograms produced:** `Emiss`, `xB`, `Q2`, `pmiss`, `Emiss_Mmiss` (2D), `omega_Emiss` (2D),
`omega_Q2` (2D).

Missing energy uses ⁴He mass to compute binding-energy-like `Emiss` via `miss_Am1`.

---

### `epp_Kinematics.cpp`

(e,e'p) and (e,e'pp) kinematics study using `many_plots` to track inclusive and per-Q²
distributions together.

**Usage:**
```
./epp_Kinematics output.root output.pdf inputfiles.hipo...
```

**Fills two `many_plots` bundles:**

- `hist_list_ep`: xB, Q², ω, θ_e, φ_e, p_lead, θ_lead, φ_lead, p_miss, m_miss, E_miss,
  θ_pq, θ_miss,q, p/q
- `hist_list_epp`: p_recoil, p_rel, θ_miss,recoil, θ_cm,rel, p_CM, p_CM components, E₂_miss

Also separates FD/CD lead angle histograms.

---

### `Electron_Cuts.cpp`

QA program comparing electron observables before and after the `clas12ana` selection chain.

**Usage:**
```
./Electron_Cuts output.root output.pdf inputfiles.hipo...
```

Runs two analyzers (`clasAna`, `clasAna2`) to compare all reconstructed electrons vs accepted
electrons. Fills: Q², xB, φ-θ coverage, HTCC nphe, PCAL energy deposition, sector-by-sector
momentum vs sampling-fraction, PCAL-vs-ECIN SF correlation, vertex z, PCAL Lv/Lw vs SF, DC
edge distances, χ²/DoF.

---

### `Proton_Cuts.cpp`

QA for FD and CD proton candidates, showing distributions before and after fiducial, vertex, and
PID cuts.

**Usage:**
```
./Proton_Cuts output.root output.pdf inputfiles.hipo...
```

**Key function:** `CD_fiducial(phi, theta, momT)` — rejects three azimuthal dead bands plus
`θ < 40°` or `θ > 125°`.

**FD cuts applied:** all three DC edge distances > 10 cm, |vz_e − vz_p| < 2 cm, |χ²_PID| < 4.  
**CD cuts applied:** `CD_fiducial`, same vertex and χ²_PID cuts.

Fills separate FD/CD QA histograms for ΔToF, β, χ²_PID, momentum, and angles, before and after
each cut stage.

---

### `Central_Detector_Proton.cpp`

Dedicated CD proton timing/fiducial study used to derive and validate the CD proton PID
parameterization.

**Usage:**
```
./Central_Detector_Proton isMC output.root output.pdf inputfiles.hipo...
```

**Key:** `Cut_Params[8]` — empirical timing band parameters for the CD. `pass_cut(mom, DT, w)`
applies `μ ± w·σ` using these parameters.

Produces a detailed QA report: vertex cuts, CVT edge cuts, ΔToF vs momentum in θ/φ slices,
β vs momentum, χ²_PID distributions, and three `TGraphErrors` fiducial-boundary graphs.

---

### `SRC_Cuts.cpp`

Full SRC cut-flow validation. Reproduces the staged (e,e'p) and (e,e'pp) selection using
corrected kinematics and shows where FD/CD protons pass or fail each step.

**Usage:**
```
./SRC_Cuts isMC A output.root output.pdf inputfiles.hipo...
```
(`A` = target mass number)

**Output:** ROOT `TTree T` of per-proton kinematics, plus histograms for all/FD/CD protons at
each cut stage (`_nocuts`, `_all`, `_fd`, `_cd`, `_final`, `_pp`).

**Variables filled:** m_miss, p_miss, k_miss_ZQ, θ_pq, θ_pmiss,q, lead/recoil flags.

---

### `KMiss.cpp`

Study of light-cone missing momentum definitions and their reconstruction resolution.

**Usage:**
```
./KMiss isMC output.root output.pdf inputfiles.hipo...
```

Computes two k_miss conventions:
- `kMissZQ` — k_miss along q (active)
- `kMissZB` — k_miss along beam z (computed but not plotted)

For MC, also computes truth-level k_miss from `mcparts()` and fills resolution histograms
`DkMissZQ_ep` and `DpMiss_ep`.

---

### `KMiss_comparison_plot.cpp`

Overlay plotter comparing `kMissZQ_ep` histograms from data and simulation ROOT files.

**Usage:**
```
./KMiss_comparison_plot output.pdf data.root sim.root
```

Scales simulation so bin 45 matches the data height, then draws both on one canvas. Output is
PDF only — no new ROOT file is written.

---

### `AvgPMiss.cpp`

Computes the per-event mean p_miss in each of the four broad pMiss bins. Used to calibrate the
quasi-free threshold lines drawn in `plot_Emiss_pmiss.py`.

**Usage:**
```
./AvgPMiss isMC A output.root inputfiles.hipo...
```

**Output ROOT tree** (`avgPMiss`): branches `selection` (ep/epp), `bin`, `lo`, `hi`,
`avg_pmiss`, `n_events`. Also prints results to stdout.

pMiss bins: `{0.4, 0.55, 0.7, 0.85, 1.0}`. Uses MC reweighter weights if `isMC`.

---

### `Q2_dependence.cpp`

Studies how Q², E_miss, and CM observables change with p_miss, separately for (e,e'p) and
(e,e'pp).

**Usage:**
```
./Q2_dependence isMC output.root output.pdf inputfiles.hipo...
```

Fills: p_miss spectrum, Q² vs p_miss, E_miss vs p_miss, vertex correlations between electron/lead/recoil,
recoil CM components per p_miss bin, recoil angular distributions, final p_miss,Rec / p_miss,SRC
ratio page.

Q² bins: `{1.5, 1.65, 1.80, 1.95, 2.10, 2.25, 2.40, 2.70, 3.00, 3.50, 5.0}` (11 fine bins).

---

### `Radiation_Offset.cpp`

Offline fitter for electron ΔE′ histograms, extracting peak shifts and widths vs electron angle
and sector to quantify radiative / energy-scale offsets.

**Usage:**
```
./Radiation_Offset isMC output.pdf input.root
```

Reads histograms named `Delta_Eprime_<thetaMin>` and `Delta_Eprime_<thetaMin>_sector_<sector>`
from the input ROOT file. Fits each angle bin (9 bins: 10–13°, 13–16°, …) with a two-Gaussian
+ exponential model (data) or a single Gaussian (MC), then plots fitted μ and σ vs angle.

**Output:** PDF only (no new ROOT file).

---

### `Electron_DeltaE.cpp`

Large calibration study of elastic (e,e'p) energy-balance residuals, used to derive
angle/phi-dependent corrections for electrons and protons in FD and CD.

**Usage:**
```
./Electron_DeltaE isMC output.root output.pdf inputfiles.hipo...
```

**Key functions:**

| Function | Description |
|---|---|
| `getGraph(TH2D*, TGraphErrors*)` | Projects each x bin of a 2D histogram, fits a Gaussian to the y projection, extracts the mean correction, fills the TGraphErrors |
| `getFunction(TGraphErrors*)` | Fits a four-parameter trigonometric function to the graph |
| `GetThetaCorrected(p)` | Evaluates pre-tabulated theta/phi correction functions |

**Key binning:**
- Electron FD θ bins: 10–37° in 3° slices
- CD θ bins: `{35, 40, 45, 50, 55, 60, 70}`
- φ bins: `{−35, −15, −5, 0, 5, 10, 15, 25, 35}`

Beam energy: 5.984792 GeV with σ = 0.00299 GeV.

---

## Analysis-Note Programs

### `Analysis_Note_Q2.cpp`

Main legacy analysis-note production program. Applies corrected kinematics and standard SRC cuts,
then writes both raw histograms and `many_plots` summary objects.

**Usage:**
```
./Analysis_Note_Q2 isMC A output.root output.pdf inputfiles.hipo...
```

**Cuts:** `xB > 1.2`, `1.5 < Q² < 5`, `p_lead > 1`, `k_miss_ZQ > 0.3`,
`0.65 < m_miss < 1.1`, `θ_lead < 37°`, FD lead only.

**`many_plots` bundles written:**
- ep: xB, Q², ω, electron/lead angles, p_miss, m_miss, E_miss, θ_pq, θ_miss,q, p/q, k_miss_ZQ/ZB
- epp: p_recoil, p_rel, θ_miss,recoil, θ_cm,rel, p_CM, p_CM components, E₂_miss

---

### `Analysis_Note_Q2_comp.cpp`

Plot-only comparison utility that overlays the `many_plots` objects from a data ROOT file and a
simulation ROOT file.

**Usage:**
```
./Analysis_Note_Q2_comp output.pdf data.root sim.root
```

Loads `many_plots` objects from both files by base name and draws comparison pages.
`many_plots::Draw_*_same` normalizes simulation to the data integral in each panel.

---

### `Analysis_Note_Q2_manyTheory.cpp`

Comparison driver overlaying data with simulation reweighted under multiple theory choices
(Kelly AV18, dipole AV18, Kelly N2LO).

**Usage:**
```
./Analysis_Note_Q2_manyTheory A output.root output.pdf data.hipo sim.hipo
```

Creates three `reweighter` objects for the three theory choices, loops over data and simulation,
fills ep/epp observables, and overlays theory curves against data. Currently only Kelly overlays
are actively plotted (dipole/N2LO hooks are partly present but commented out).

---

## Main Figures Programs

### `Main_Figs.cpp`

Earlier hand-written version of the main figures analysis, later superseded by
`Main_Figs_Binned.cpp`.

**Usage:**
```
./Main_Figs isMC A output.root output.pdf inputfiles.hipo...
```

Fills hard-coded histograms: Q² yields, p_miss/k_miss distributions and ratios, epp CM
components vs Q², E₀/E₁/E₂_miss vs Q² in p_miss and k_miss bins, angular QA histograms.

Also uses `getGraph(TFile*, TCanvas*, ...)` to project each x bin of a 2D histogram, fit a
Gaussian, and write `g_mu_*` / `g_sigma_*` graphs.

---

### `Main_Figs_Binned.cpp`

**Modern refactored main figures analysis.** Computes event kinematics once, defines all
histograms declaratively via `BinnedHistStore`, fills a nominal histogram store plus 100 toy
systematic stores, and writes flat ROOT tables.

**Usage:**
```
./Main_Figs_Binned isMC A output.root [--mode legacy|modern] \
    [--q2-reweight weights.root] inputfiles.hipo...
```

**Output ROOT file contains:**
- `diffTable` — `TTree` with one row per (task, selection, selector_bins, value_bin)
- `integratedTable` — `TTree` with integrated yields
- `hists/nominal/` — nominal histograms
- `hists/toy_000/` through `hists/toy_099/` — systematically varied histograms

**Key structs:**

| Struct | Description |
|---|---|
| `EventKinematics` | Stores all per-event derived variables (corrected 4-vectors, p_miss, m_miss, k_miss, CM quantities, etc.) |
| `CutVariation` | Stores nominal thresholds and randomized toy variations; `apply(ek)` sets `passep`/`passepp` |

**`computeEventKinematics(...)` — central function** — performs corrected kinematic
reconstruction and returns an `EventKinematics` struct.

**`buildFillTasks(legacyCompatMode)`** — defines all histogram families:
- Plain yields: p_miss, k_miss, Q², p_lead, p_rel, angular distributions
- Analysis-note style: 50-bin uniform p_miss/k_miss distributions
- E₀/E₁/E₂_miss in p_miss/k_miss/Q² bins
- p_CM-x/y/z in Q² bins
- Q²-vs-p_miss/k_miss 2D yields

**Core bin edges:**
```
bE_Q2   = {1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 5.0}
bE_pmiss = {0.4, 0.55, 0.7, 0.85, 1.0}
bE_kmiss = {0.3, 0.45, 0.6, 0.75, 0.9}
```

`--mode legacy` reproduces older analysis behavior including legacy cut randomization.

---

### `Main_Figs_Binned_StatOnly.cpp`

Fast nominal-only companion to `Main_Figs_Binned.cpp`. Skips toys, runs `BuildGraphs` internally,
and adds fit-derived statistical-only graphs.

**Usage:** Mirrors `Main_Figs_Binned` (same flags, same output layout).

In addition to the standard output, also writes `TGraphErrors` for fit-derived quantities using
**`statOnlyBuildFitSpecs`**:
- `sigma_pcmx`, `sigma_pcmy` — Gaussian-fit width of p_CM-x/y vs Q² bin
- `stddev_` and `mean_` graphs for E_miss histograms in p_miss/k_miss bins

Uses Gaussian fits for p_CM widths and plain ROOT `GetMean()`/`GetStdDev()` for E_miss moments.

---

### `Main_Figs_Sys_Err.cpp`

**Reference implementation** of the main figures workflow with 100 toy variations for systematic
uncertainty estimation. Superseded by `Main_Figs_Binned.cpp` but kept as the physics reference.

**Usage:**
```
./Main_Figs_Sys_Err isMC A output.root output.pdf inputfiles.hipo...
```

**Key functions:**

| Function | Description |
|---|---|
| `runEvent(...)` | Computes corrected ep/epp kinematics |
| `CutRandom(...)` | Applies nominal or randomized cuts |
| `setUpHistGroup` / `fillUpHistGroup` | Books / fills one complete set of histograms (one call per toy) |
| `getGraphWithError` | Builds a `TGraphErrors` whose y-errors are `sqrt(sys² + stat²)` from the 100 toy histograms |
| `getMeanStddev` | Computes mean and std dev of 100 toy histogram values for a given bin |

**Known legacy issue:** `TRandom3(0)` for cut randomization is re-seeded inside the event loop,
which means each toy run produces the same random sequence. Fixed in `Main_Figs_Binned.cpp`.

---

### `BuildGraphs.cpp`

Post-processing utility that converts the flat `diffTable` and `integratedTable` `TTree`s
produced by `Main_Figs_Binned.cpp` into `TGraphAsymmErrors` objects.

**Usage:**
```
./BuildGraphs file.root
```

Opens the file in UPDATE mode. Reads both tables and writes graphs into the `graphs/` directory.

**Graph naming:** matches `graph_names.py` exactly — `"<task>|<selection>|<axis_bins>"` for
differential and `"<task>|<selection>|<pattern>"` for integrated.

**Error convention:** total y-error = `sqrt(stat² + sys_up² )` / `sqrt(stat² + sys_down²)`.
Falls back to symmetric `sys_error` if asymmetric branches are absent.

**Key functions:**

| Function | Description |
|---|---|
| `buildDiffGraphs(diffTree, outDir)` | Groups rows by (task_name, selection, axis_bin), sorts by value_bin, writes one graph per group |
| `buildIntegratedGraphs(intTree, outDir)` | Groups rows by (task_name, selection, pattern), writes 1D graphs |
| `joinBins(vector<int>)` | Creates axis-bin suffixes like `"2_1"` |

---

## σ_CM Extraction Pipeline

### `Main_sigmaCM.cpp`

Stage-1 producer: fills reco-level C.M.-momentum histograms for data and a bank of 100 MC
templates, each reweighted to a different assumed input σ_CM value.

**Usage:**
```
./Main_sigmaCM A output.root data.hipo sim.hipo
```

**Output ROOT file:** integrated and Q²-binned histograms for p_CM-x, y, z, T. For each of the
100 σ_CM values, a separate set of reweighted simulation histograms is stored.
σ_CM scan range: 0.08 to 0.25 GeV (100 steps).

**Key functions:**

| Function | Description |
|---|---|
| `runEvent(...)` | Computes epp CM coordinates and selection flag |
| `getChi2(data, sim, normFactor, fitRange)` | Scans normalization factor to minimize χ² between data and one template |
| `getG` / `getGraph` | Fit/convert Q²-binned histograms to mean/sigma graphs |

---

### `Main_sigmaCM_Hists.cpp`

Stage-2 analysis: reads the `Main_sigmaCM.cpp` output, builds χ²-vs-σ curves for each CM
component, and estimates best-fit σ and ±1 χ² intervals.

**Usage:**
```
./Main_sigmaCM_Hists output.root output.pdf CutRangeUpper input.root
```

**Key functions:**

| Function | Description |
|---|---|
| `getChi2(data, sim, normFactor)` | Scans norm and returns χ² for one data/template pair |
| `getValue(graph, min_sigma, max_sigma, chi2NDF, center, lower, upper)` | Finds graph minimum and points where χ² rises by 1 (±1σ confidence interval) |

Parameters: `linbin = 100`, `min_sigma = 0.08`, `max_sigma = 0.25`.

---

### `Main_sigmaCM_WidthMatch.C`

**Cleaner ROOT macro** for the σ_CM extraction. Instead of fitting a χ² parabola, it measures a
width proxy in data and in each template, inverts the monotonic reco-width vs input-σ map, and
reports the matched σ_CM.

**Usage:**
```
root -l -b -q 'Main_sigmaCM_WidthMatch.C("stage1.root","out.root")'
```

**Key structs/functions:**

| Item | Description |
|---|---|
| `Component` | Describes one CM momentum component: histogram naming, fit range, display label |
| `FitWidth` / `MatchResult` | Store width extraction and matched-σ results |
| `fitWidth(hist)` | Measures RMS over a restricted axis range as the width proxy |
| `isMonotone(widths)` | Checks that template widths vary monotonically with input σ |
| `invertLinear(widths, sigmas, dataWidth)` | Inverts the width map by linear interpolation |
| `matchOne(component, q2_bin)` | Full pipeline for one component/Q² bin: builds width map → extracts data width → inverts → derives asymmetric errors |
| `drawOverlay(component, q2_bin, result)` | Writes overlay canvas comparing data to the nearest matched template |

σ_CM scan range: 0.08 to 0.35 GeV; `linbin = 100`.

---

### `Q2_Reweight.cpp`

Derives a Q² reweight map by matching accepted data and MC Q² spectra under the same nominal cuts.

**Usage:**
```
./Q2_Reweight output.root simA simType [--selection ep|epp] [--max N] \
    data.hipo ... --sim sim.hipo ...
```

**Output ROOT file** (`q2_reweighting/` directory):
- `h_Q2_data_raw`, `h_Q2_sim_raw` — raw Q² histograms
- `h_Q2_data_norm`, `h_Q2_sim_norm` — normalized versions
- `h_Q2_reweight` — bin-by-bin `data/sim` weight histogram
- `h_Q2_sim_reweighted` — reweighted simulation histogram for QA

Uses `computeEventKinematics` and `CutVariation::Nominal` from `Main_Figs_Binned.cpp` so the
selection is identical. Q² bins are exactly the `bE_Q2` edges, compatible with `Q2Reweight.h`.
