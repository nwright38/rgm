# myPlots — Script Reference

This directory contains Python plotting scripts and ROOT macros that read the output of the
`Ana/Q2_Ana` analysis pipeline and produce publication-quality figures. The workflow is:

```
Main_Figs_Binned  →  output.root
BuildGraphs       →  output.root (adds graphs/ directory)
[ExtractFitQuantities.C] →  output.root (adds fit-derived graphs)
Python driver scripts  →  PDF figures
```

---

## Utility / Helper Modules

These modules contain no executable code of their own; they are imported by the driver scripts.

---

### `graph_io.py`

Handles all ROOT file reading. Every other module works with plain Python lists — `graph_io` is the
only place a ROOT file is ever opened.

Uses **uproot** (pure Python), not PyROOT, so these scripts can run locally without a ROOT
installation. The C++ side (`Main_Figs_Binned`, `BuildGraphs`) still requires ROOT on the cluster.

**Functions:**

| Function | Arguments | Returns |
|---|---|---|
| `open_file(path)` | file path string | uproot file handle |
| `read_graph_asymm(root_file, name)` | uproot file handle, graph name string | `(x, y, yerr_low, yerr_high)` as plain lists |
| `read_graph(root_file, name)` | uproot file handle, graph name string | `(x, y, yerr)` — symmetric, `yerr = max(low, high)` per point |

`read_graph_asymm` reads graphs from the `graphs/` subdirectory of the file (written by
`BuildGraphs.cpp`). If no `graphs/` directory is found, it raises a `KeyError` explaining that
`BuildGraphs` needs to be run first. Both `TGraphErrors` (symmetric) and `TGraphAsymmErrors`
(asymmetric) are handled; symmetric graphs are promoted to asymmetric form automatically.

---

### `graph_names.py`

The single place that encodes the naming convention `BuildGraphs.cpp` uses when writing graphs.
**Always use these functions instead of hand-building name strings** — if the convention ever
changes, only this file and `BuildGraphs.cpp` need updating.

**Convention:**
```
differential graph  : "<task_name>|<selection>|<b0>_<b1>_..."
                      (or "<task_name>|<selection>|none" for tasks with no selector axes)
integrated graph    : "<task_name>|<selection>|<pattern>"
```

**Constants:**
- `PATTERN_NONE` = `"collapse:none"` — one point per selector bin (no collapsing)
- `PATTERN_ALL` = `"collapse:all"` — all selector axes collapsed to a single total

**Functions:**

| Function | Arguments | Returns |
|---|---|---|
| `pattern_collapse_axis(axis_name)` | axis name, e.g. `'Q2'` | pattern string for collapsing that axis |
| `diff_graph_name(task_name, selection, axis_bin)` | task name; `'ep'` or `'epp'`; list of bin indices | graph name string |
| `integrated_graph_name(task_name, selection, pattern)` | task name; selection; pattern string | graph name string |
| `ratio_diff_graph_name(numerator_task, denominator_task, axis_bin)` | two task names; list of bin indices | graph name string for the ratio graph |

`axis_bin` order matches the order axes were listed in the `FillTask` definition inside
`Main_Figs_Binned.cpp` for that task.

---

### `labels.py`

A single constant:

```python
PMISS_LABELS = [
    r'$0.4<p_{Miss}<0.55$',
    r'$0.55<p_{Miss}<0.7$',
    r'$0.7<p_{Miss}<0.85$',
    r'$0.85<p_{Miss}<1.0$',
]
```

Must match the `bE_pmiss` bin edges in `Main_Figs_Binned.cpp`. Import and use in any script that
labels pMiss-bin panels.

---

### `normalization.py`

Computes the scale factors used to put different files on the same display scale in overlay plots
(i.e., nucleus / He or simulation / data).

Reads the integrated Q² yield straight from the `integratedTable` graphs that `BuildGraphs.cpp`
wrote, so the normalization is always the correct total event count.

**Functions:**

| Function | Arguments | Returns |
|---|---|---|
| `scale_factor(target_file, reference_file, selection)` | two uproot file handles; `'ep'` or `'epp'` | `(factor, error)` — multiply target by `factor` to match reference scale |

Example:
```python
import normalization
factor, err = normalization.scale_factor(f_sim, f_data, 'ep')
```

---

### `nuclei_config.py`

Plain-data configuration: which ROOT files exist and how each nucleus / simulation model should be
styled. No files are opened here — driver scripts call `graph_io.open_file(entry['file'])` when
they actually need one.

**Key variables:**

| Variable | Description |
|---|---|
| `DATA_DIR` | Path to multi-nucleus data ROOT files on the cluster |
| `SIM_DIR` | Path to simulation ROOT files on the cluster |
| `N_DATA_DIR` | Path to the He reference data and AV18 sim ROOT files |
| `REFERENCE_KEY` | `'He'` — the nucleus everything else is normalized to |
| `NUCLEI` | List of dicts, one per nucleus: `key`, `label` (LaTeX), `file`, `color`, `marker`, `offset` |
| `SIM_MODELS` | List of dicts, one per sim model: `key`, `label`, `file`, `color` |

**Functions:**

| Function | Arguments | Returns |
|---|---|---|
| `find(entries, key)` | list of dicts; key string | matching dict, raises `KeyError` if not found |
| `reference_nucleus()` | — | the He entry from `NUCLEI` |

To add a nucleus or model, add one entry to `NUCLEI` or `SIM_MODELS` — no driver script changes needed.

---

### `ratio.py`

One pure function for computing ratios with propagated errors. Used for data/sim ratio subpanels
and for normalization factors. Data and simulation are treated as independent samples (standard
Gaussian error propagation).

**Functions:**

| Function | Arguments | Returns |
|---|---|---|
| `ratio_with_error(y1, e1, y2, e2)` | four scalars | `(y1/y2, error)` — returns `(0, 0)` if either value is 0 |
| `ratio_series_with_error(y1_list, e1_list, y2_list, e2_list)` | four equal-length lists | `(ratios, errors)` — applied bin-by-bin |

---

### `plot_primitives.py`

Low-level matplotlib drawing functions. These only take plain `(x, y, yerr)` lists and never touch
a ROOT object. Scaling and x-shifting are applied by the caller before passing to these.

**Functions:**

| Function | Arguments | Description |
|---|---|---|
| `step_with_error(ax, x, y, yerr_low, yerr_high, color)` | matplotlib axis + data | Data-style: step line + plain error bars, no marker |
| `errorbar_marker(ax, x, y, yerr_low, yerr_high, color, marker)` | matplotlib axis + data | Data-style with a marker symbol, for cross-nucleus overlays |
| `line_with_band(ax, x, y, yerr_low, yerr_high, color, alpha=0.3)` | matplotlib axis + data | Simulation-style: solid line + shaded error band |
| `step_with_error_q2(ax, x, y, yerr_low, yerr_high, color)` | matplotlib axis + data | Same as `step_with_error`, but pads the first/last point out to x = 0 / 10 so the step fills the full Q² panel |
| `line_with_band_q2(ax, x, y, yerr_low, yerr_high, color, alpha=0.3)` | matplotlib axis + data | Same as `line_with_band` with the same Q² padding |

---

### `plot_helpers.py`

Higher-level building blocks layered on top of `graph_io`, `graph_names`, `plot_primitives`, and
`ratio`. This is the main import for all driver scripts.

#### `Series` class

Represents one thing to draw on a plot.

```python
Series(label, root_file, color, kind,
       scale_ep=1.0, scale_epp=1.0, offset=0.0, marker=None)
```

| Attribute | Description |
|---|---|
| `label` | Legend label (can be LaTeX) |
| `file` | uproot file handle |
| `color` | matplotlib color string |
| `kind` | `'data'` → step+errorbars; `'sim'` → line+band; `'data_ex'` → marker with x-offset (other nuclei) |
| `scale_ep` / `scale_epp` | normalization factors for ep / epp plots; ignored for `selection='ratio'` |
| `offset` | base x-shift in data units (scaled by `offset_scale` in draw calls) |
| `marker` | matplotlib marker string (only used for `'data_ex'` kind) |

#### Key functions

| Function | Description |
|---|---|
| `get_xy_err(series, task_name, selection, axis_bin=(), ...)` | Reads, scales, and x-shifts a graph for one Series. Returns `(x, y, yerr)`. Pass `integrated=True` + `pattern=...` to read from the integrated table. |
| `get_xy_err_asymm(...)` | Same but returns `(x, y, yerr_low, yerr_high)`. |
| `draw(ax, series, task_name, selection, axis_bin=(), ...)` | Draws one Series onto a matplotlib axis using the right primitive for its `kind`. `q2_panel=True` uses the wide-step/band padding. |
| `annotate(ax, items)` | Places text labels from a list of dicts with keys `x`, `y`, `text`, optional `color`, `fontsize`, `bbox`. |

#### High-level figure functions

| Function | Description |
|---|---|
| `plot_overlay(pdf, task_name, series, xlabel, ylabel, xlim, ylim, ...)` | Single-axis figure with N overlaid Series and an optional data/sim ratio subpanel underneath. Set `with_ratio=False` for ratio-quantity or A-dependence plots. `unit_scale=True` skips normalization (use for ratio histograms). |
| `plot_emiss_4x2(pdf, series, var_label, task_ep, task_epp, ...)` | 4-row × 2-col E_miss panel: (e,e'p) on the left, (e,e'pp) on the right, one pMiss bin per row. |
| `plot_emiss_4x2_ratio_only(pdf, ref, sims, var_label, task_ep, task_epp, ...)` | Same layout but shows only data/sim ratios. |
| `plot_emiss_4x1(pdf, series, var_label, task_epp, ...)` | 4-row × 1-col (e,e'pp) E_miss panel. |
| `plot_emiss_4x1_ratio_only(pdf, ref, sims, var_label, task_epp, ...)` | Same but data/sim only. |
| `plot_emiss_4x2_note(pdf, series, var_label, task_ep, task_epp, ...)` | Analysis-note style: boxed pMiss-bin label centered above each row, optional threshold line at quasi-free E_miss. |
| `plot_emiss_4x1_note(pdf, series, var_label, task_epp, ...)` | Same for 4×1. |
| `plot_emiss_waterfall(pdf, pmiss_bins, data, sims, task_name, ...)` | Waterfall overlay: all four pMiss bins on the same axis, offset vertically for clarity. |
| `plot_q2_single(pdf, data, sims, task_name, ylabel, ...)` | Single Q²-vs-quantity panel with data step and simulation lines. |
| `plot_q2_2x2(pdf, data, sims, task_name, ylabel, pmiss_labels, ...)` | 2×2 Q²-vs-quantity panel, one pMiss bin per panel. |

All figure functions return the `fig` object and always save to the multi-page `pdf` book. Pass
`save_as='path/to/file.pdf'` to also write the figure as an individual PDF.

---

### `plotTOOL.py`

**Older PyROOT-based utilities**, still used by `Sys_Err.py`. All functions take a PyROOT
`TFile` and a histogram/graph name string as inputs and call PyROOT's `TGraph::GetX()` / `GetY()`
/ `GetErrorY()` directly.

The newer scripts (`makePlots.py` and all `plot_*.py`) use `graph_io` + `plot_helpers` instead.
`plotTOOL.py` is kept as-is for backward compatibility.

**Key functions:**

*TGraph plotting (step/line):*

| Function | Description |
|---|---|
| `plotTGEStep(ax, file, name, color)` | Reads a TGraph from a ROOT file via PyROOT and draws it as a step histogram with error bars |
| `plotTGEStepEx(ax, file, name, color, factor, shift, marker)` | Same but scales y by `factor` and shifts x by `shift`, drawing markers instead of a step line (used for cross-nucleus overlays) |
| `plotTGEStepQ2(ax, file, name, color)` | Same as `plotTGEStep` but pads x to 0 and 10 for Q² panel style |
| `plotTGELine(ax, file, name, color, factor)` | Reads a TGraph and draws it as a solid line with a shaded error band scaled by `factor` (simulation style) |
| `plotTGELineQ2(ax, file, name, color)` | Same but with Q² panel x-padding |
| `plotTGraphAsymmErrors(file, name, lineColor, fillColor)` | Draws a `TGraphAsymmErrors` as a filled band + line |

*Data extraction helpers:*

| Function | Description |
|---|---|
| `getXforGraph*(file, name)` | Extracts x points from a named graph |
| `getYforGraph*(file, name)` | Extracts y points from a named graph |
| `getZforGraph*(file, name)` | Extracts z points from a `TGraph2D` |
| `getXforHist(file, name)` | Extracts bin centers from a named histogram |
| `getYforHist(file, name)` | Extracts bin contents from a named histogram |
| `getYErrforHist(file, name)` | Extracts bin errors from a named histogram |
| `printHist(file, name)` | Prints histogram bin contents to stdout |

*Normalization helpers:*

| Function | Description |
|---|---|
| `getFactor_ep(file, reference)` | Returns the ep normalization factor by scanning bins 1–7 of `h_Q2_ep_SRC_Q2_i` — **note: known bug** where the loop overwrites rather than accumulates, so result reflects only the last Q² bin |
| `getFactor_epp(file, reference)` | Same for epp, same bug |
| `getFactorTGraph(file, name, reference, ref_name)` | Returns the normalization factor for a single named histogram/graph, used for per-Q²-bin normalization (no bug) |
| `getFactorHist(file, name, reference, ref_name)` | Same for `TH1` histograms |
| `getFactor(file, name, reference, ref_name)` | Generic version |

*Special plot layouts:*

| Function | Description |
|---|---|
| `plotSig(pdf, Q2bins, dataFile, simFile, component, xmin, xmax, ylabel, xlabel, factor_array)` | Legacy pseudo-3D Q² waterfall: plots the C.M. momentum distributions per Q² bin with Gaussian fits, plus a summary σ-vs-Q² page |
| `FindMin(graph)` | Scans a graph for its minimum and prints approximate ±1 confidence-interval bounds |

*Contour/2D helpers:* `plotContourTest`, `plotContourSpacing`, `plotNC`, `plotVecPoint`, `fixZ`

---

## Driver Scripts — Python

These scripts are run directly and produce PDF output.

---

### `Sys_Err.py`

**Older combined-plots script** that uses `plotTOOL.py` + PyROOT directly. Produces a single
multi-page PDF (`plots/output.pdf`) containing all the main figures for a He data sample compared
against multiple simulation models and other nuclei.

**Inputs:** Hardcoded ROOT file paths at the top of the script.
**Output:** `plots/output.pdf`

**Configuration:**
```python
Q2bins = [1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0]
pMbins = [0.4, 0.55, 0.7, 0.85, 1.0]
kMbins = [0.3, 0.45, 0.6, 0.75, 0.9]
```

**What it produces (in order):**

1. p_CM·x and p_CM·y distributions (He data + GCF/N2LO simulation)
2. epp/ep ratio vs p_miss (He data + 5 simulation models: AV18, AV4', N2LO×2, Norfolk)
3. E₁_miss 4×2 panel in pMiss bins (ep left, epp right)
4. E₂_miss 4×1 panel in pMiss bins (epp only)
5. σ_CM vs Q² for x and y components (via `plotTOOL.plotSig`)
6. σ_CM-x and σ_CM-y vs Q² (step + line panels)
7. epp/ep vs Q² in 4 pMiss bins (2×2 panel)
8. Std dev E₁_miss vs Q² in pMiss bins (ep and epp separately)
9. Std dev E₂_miss vs Q² in pMiss bins (epp)
10. Mean E₁_miss and E₂_miss vs Q² in pMiss bins
11. p_CM·x and p_CM·y cross-nucleus overlays (He + C + Ca40 + Sn)
12. epp/ep vs p_miss cross-nucleus (He + C + Ca40 + Ca48 + Sn)
13. E₁_miss 4×2 cross-nucleus
14. E₂_miss 4×1 cross-nucleus

This script is **not command-line driven** — edit the file path variables at the top to change inputs.

---

### `makePlots.py`

**Newer combined-plots script** that replaces `Sys_Err.py`'s manual per-figure boilerplate with
`plot_helpers` function calls. Produces the same figures but is shorter and more maintainable.

**Inputs:** Hardcoded ROOT file paths near the top of the script (mix of cluster paths and local
file names).  
**Output:** `plots/output.pdf`

**Same configuration as `Sys_Err.py`:**
```python
Q2bins = [1.5, 1.80, 2.10, 2.40, 2.70, 3.00, 3.50, 5.0]
pMbins = [0.4, 0.55, 0.7, 0.85, 1.0]
```

Uses `plotTOOL.py` for normalization factors (`getFactor_ep`, `getFactor_epp`,
`getFactorTGraph`, `plotSig`) and `plot_helpers` for all figure generation.

Defines `Series` objects at the top for He data, all five simulation models, and the four non-He
nuclei, then calls `ph.plot_overlay`, `ph.plot_emiss_4x2`, `ph.plot_emiss_4x1`, and `ph.plot_q2_*`
to produce the same page sequence as `Sys_Err.py`.

---

### `plot_1D.py`

Writes individual single-figure PDFs for 14 standard 1D distributions under `pdf/1D/`, one file
per variable. Intended for dropping figures directly into a note without running the full
combined script.

**Usage:**
```
python plot_1D.py Data_He.root
python plot_1D.py Data_He.root Sim_He_AV18.root --sim-label AV18
python plot_1D.py Data_He.root Sim_AV18.root Sim_N2LO.root \
    --sim-label AV18 --sim-label N2LO --sim-color red --sim-color blue
```

**Arguments:**
- `data_file` — Main_Figs_Binned + BuildGraphs output for data
- `sim_files` (0–4) — simulation outputs (same pipeline)
- `--sim-label`, `--sim-color` — one per sim file, matched in order
- `--data-label`, `--data-color`
- `--out-dir` — default `pdf/1D`

**Variables produced:**

| Output file | Selection | Variable |
|---|---|---|
| `pMiss_ep.pdf` | ep | p_miss |
| `pMiss_epp.pdf` | epp | p_miss |
| `kMiss_ep.pdf` | ep | k_miss |
| `kMiss_epp.pdf` | epp | k_miss |
| `pRel_epp.pdf` | epp | p_rel |
| `pcmz_epp.pdf` | epp | p_CM·z |
| `pcmx_epp.pdf` | epp | p_CM·x |
| `pcmy_epp.pdf` | epp | p_CM·y |
| `theta_pmiss_ep.pdf` | ep | θ(p_miss, q) |
| `theta_pmiss_epp.pdf` | epp | θ(p_miss, q) |
| `q_ep.pdf` | ep | \|q\| |
| `q_epp.pdf` | epp | \|q\| |
| `theta_pLeadq_ep.pdf` | ep | θ(p_lead, q) |
| `theta_pLeadq_epp.pdf` | epp | θ(p_lead, q) |

Each plot has an optional data/sim ratio subpanel. y-limit is set automatically to 1.2× the
histogram peak. Simulation is normalized to data via `normalization.scale_factor` on the ep yield;
the same factor is applied to epp histograms.

---

### `plot_A_dependence.py`

Cross-nucleus comparisons (He + C + ⁴⁰Ca + ⁴⁸Ca + Sn) written to `plots/A_dependence.pdf`.
There is no simulation reference here, so `with_ratio=False` throughout.

**No command-line arguments** — file paths come from `nuclei_config.py`. Edit that file to change
which nuclei / files are used.

**Figures produced:**
1. p_CM·x overlay (He + C + ⁴⁰Ca + Sn, ⁴⁸Ca omitted)
2. p_CM·y overlay (same)
3. epp/ep vs p_miss, all 5 nuclei (already a ratio — no extra scaling)
4. E₁_miss 4×2 panel (ep + epp), all 5 nuclei
5. E₂_miss 4×1 panel (epp), all 5 nuclei

Normalization is computed via `normalization.scale_factor` against He for each non-He nucleus.

---

### `plot_data_vs_sim.py`

Comprehensive data vs simulation comparison, command-line driven. Accepts one data file and any
number of simulation files, all pre-processed by Main_Figs_Binned + BuildGraphs. Writes a combined
multi-page PDF to `--out` (default `plots/data_vs_sim.pdf`).

**Usage:**
```
python plot_data_vs_sim.py Data_He.root
python plot_data_vs_sim.py Data_He.root Sim_He_AV18.root --sim-label 'AV18' --sim-color red
python plot_data_vs_sim.py Data_He.root Sim_He_AV18.root Sim_He_N2LO.root \
    --sim-label 'AV18' --sim-label 'N2LO' --sim-color red --sim-color blue \
    --out plots/he_multi.pdf
```

**Arguments:**
- `data_file` — data ROOT file
- `sim_files` (0 or more) — simulation ROOT files
- `--sim-label`, `--sim-color` — per sim file
- `--data-label`, `--data-color`
- `--out` — output PDF path

**Figures produced** (each has an optional data/sim ratio subpanel):

- p_CM·x and p_CM·y distributions (epp)
- epp/ep vs p_miss and vs k_miss (ratio task — no subpanel)
- Q² yields (ep and epp, vs Q²)
- p_miss distributions per Q² bin (2×4 panel)
- k_miss distributions per Q² bin (2×4 panel)
- E₁_miss 4×2 (ep + epp, pMiss bins)
- E₁_miss 4×2 data/sim ratio panels
- E₂_miss 4×1 (epp, pMiss bins)
- E₂_miss 4×1 data/sim ratio panels
- E₁_miss 4×2 (kMiss bins)
- E₂_miss 4×1 (kMiss bins)
- If fit-derived graphs are present (written by `ExtractFitQuantities.C`):
  - σ_CM-x and σ_CM-y vs Q²
  - epp/ep vs Q² in pMiss bins (2×2 panel)
  - Std dev / Mean E₁_miss (ep and epp) vs Q² in pMiss bins
  - Std dev / Mean E₂_miss (epp) vs Q² in pMiss bins
  - Same panels repeated for kMiss bins

Individual figures can be exported as their own PDFs by adding entries to `SAVE_AS` at the top
of the script (a `dict` mapping a plot key to an output path).

---

### `plot_Emiss_pmiss.py`

E_miss distributions in pMiss bins for the analysis note. Writes individual single-figure PDFs
under `pdf/Emiss_pmiss/`.

**Usage:**
```
python plot_Emiss_pmiss.py Data_He.root
python plot_Emiss_pmiss.py Data_He.root Sim_He_AV18.root --sim-label AV18
```

**Arguments:** Same as `plot_1D.py` (data file, optional sim files, labels/colors, `--out-dir`).

**What it produces:**

| Output file | Description |
|---|---|
| `E1miss_4x2.pdf` | E₁_miss in pMiss bins, analysis-note style (boxed bin label, optional threshold line) |
| `E0miss_4x1_ep.pdf` | E₀_miss in pMiss bins for (e,e'p), analysis-note style |
| `E2miss_4x1.pdf` | E₂_miss in pMiss bins for (e,e'pp) |
| `E0miss_waterfall_ep.pdf` | E₀_miss waterfall: all four pMiss bins overlaid on the same axis |

Uses analysis-note figure style (`plot_emiss_4x2_note`, `plot_emiss_4x1_note`,
`plot_emiss_waterfall` from `plot_helpers`) with a single boxed pMiss-bin range label centered
above each row.

**Key constant — mean pMiss per bin:**
```python
PMISS_CENTERS_EP  = [0.499816, 0.63018, 0.769554, 0.91613]
PMISS_CENTERS_EPP = [0.508686, 0.634642, 0.772566, 0.915738]
```
These are measured from data by `AvgPMiss.cpp` and hardcoded here. Rerun `AvgPMiss.cpp` and update
these if the cuts, data sample, or pMiss bin edges change.

---

### `plot_Q2_pmiss_extracted.py`

Writes individual single-figure PDFs for Q²-binned and pMiss-binned extracted quantities
(epp/ep ratio, E_miss σ and mean vs Q²) under `pdf/Q2_pmiss/` (or `--out-dir`).

**Usage:**
```
python plot_Q2_pmiss_extracted.py Data_He.root
python plot_Q2_pmiss_extracted.py Data_He.root Sim_He_AV18.root --sim-label AV18
python plot_Q2_pmiss_extracted.py Data_He.root Sim.root --sim-label GCF \
    --extract-fit   # runs ExtractFitQuantities.C first
```

**Arguments:** Same pattern as other driver scripts plus:
- `--extract-fit` — runs `ExtractFitQuantities.C` on each file via `root -l -b -q` before
  plotting (adds the fit-derived graphs to each file if they aren't there already)

**Figures produced:**
- epp/ep vs Q² in each pMiss bin (2×2 panel)
- σ_CM-x and σ_CM-y vs Q² (one panel each)
- Std dev and mean E₁_miss (ep and epp) vs Q² in pMiss bins (2×2 panels each)
- Std dev and mean E₂_miss (epp) vs Q² in pMiss bins

---

### `plot_ratio.py`

epp/ep ratio vs p_miss and vs k_miss, written as individual PDFs under `pdf/ratio/`.

**Usage:**
```
python plot_ratio.py Data_He.root
python plot_ratio.py Data_He.root Sim_He_AV18.root --sim-label AV18
```

**Arguments:** Data file + optional sim files + labels/colors + `--out-dir`.

**Outputs:**
- `pMiss_ratio.pdf` — epp/ep vs p_miss (ratio task, no normalization scaling needed)
- `kMiss_ratio.pdf` — epp/ep vs k_miss

Uses `selection='ratio'` throughout; `plot_helpers._scale_for` always returns 1.0 for ratio
selections regardless of what scale factors are attached to a `Series`.

---

### `plot_pcmz_q2_scratch.py`

Scratch visualization: overlays p_CM·x (configurable, see `TASK_NAME`) distributions in
each Q² bin, one panel per Q² bin. Simulation is scaled independently per Q² bin to match the
data yield in that panel (emphasizing shape differences, not normalization).

**Usage:**
```
python myPlots/plot_pcmz_q2_scratch.py Data.root Sim.root \
    --sim-label GCF --out myPlots/pdf/pcmz_q2_scratch.pdf
```

**Arguments:**
- `data_file`, `sim_files` (0 or more)
- `--sim-label`, `--sim-color`
- `--data-label`, `--data-color`
- `--out` — output PDF path
- `--xlim lo hi` — x-axis limits (default −0.5 to 0.5)
- `--unit-area` — normalize each panel to unit area instead of matching yields

Writes a single PDF with one page per Q² bin (7 bins total for `bE_Q2`), plus an overlay page
showing all Q² bins on the same axis using different colors.

---

## ROOT Macros (C)

These are run with `root -l -b -q '...'` and read/write ROOT files directly. They are not part of
the Python pipeline but complement it with ROOT-native fit and plot operations.

---

### `ExtractFitQuantities.C`

Post-processing macro that derives fit quantities that `Main_Figs_Binned.cpp` deliberately left
out of its flat tables (because they require fitting, not simple bin sums). Reads the
`hists/nominal` and `hists/toy_NNN` directories written by `Main_Figs_Binned`, and writes results
directly into the `graphs/` directory of the same file, using the **same naming convention** as
`graph_names.py` so the Python plotting scripts can read them without modification.

**Usage:**
```
root -l -b -q 'ExtractFitQuantities.C("Data_He.root")'
```

**Quantities defined in `buildSpecs()`:**
- `sigma_pcmx|epp|<q2_bin>` — Gaussian-fit σ of p_CM·x in each Q² bin
- `sigma_pcmy|epp|<q2_bin>` — Gaussian-fit σ of p_CM·y in each Q² bin
- `stddev_E1miss_ep_SRC_pmiss|ep|<pmiss_bin>_<q2_bin>` and epp variant
- `mean_E1miss_ep_SRC_pmiss|ep|<pmiss_bin>_<q2_bin>` and epp variant
- Same for E₂_miss and for kMiss bins

**Internal spec machinery:**

| Item | Description |
|---|---|
| `enum Method` | `FIT_SIGMA` (Gaussian fit) or `PLAIN_STDDEV` / `PLAIN_MEAN` (ROOT `TH1::GetStdDev` / `GetMean`) |
| `struct FitSpec` | Describes one extraction: histogram name template, selector axes (pMiss/kMiss/none), Q² axis, fit range, extraction method, output graph name prefix |
| `buildSpecs()` | Returns a list of all `FitSpec` entries to process |
| `extract(spec, file, nominalDir, toyDirs)` | For each axis bin × Q² bin, collects toy values, computes toy mean and toy std dev, combines with nominal stat error in quadrature, writes one `TGraphErrors` |

**Key constant:** `kNToys = 100`

---

### `plot_extraction_diagnostics.C`

Sanity-check macro that produces a multi-page diagnostic PDF from a `Main_Figs_Binned` output file.
Does not modify any ROOT files.

**Usage:**
```
root -l -b -q 'plot_extraction_diagnostics.C("Data_He_hists_6GeV.root")'
root -l -b -q 'plot_extraction_diagnostics.C("Data_He_hists_6GeV.root", "pdf/Q2_pmiss/extraction_diagnostics.pdf")'
```

**Diagnostic pages produced:**
1. E_miss distributions in pMiss bins (2×2 pages) — annotated with mean and σ
2. E_miss distributions split by pMiss and Q² (3×3 pages) — annotated with mean and σ
3. σ_CM fit diagnostics — per-Q²-bin p_CM histograms with Gaussian fit overlay and fit
   metrics (μ, σ, χ²/ndf) shown on each panel

---

### `plot_sigmaCM.C`

Reads the output of the σ_CM fitting pipeline and produces a multi-page PDF showing σ_CM for
each C.M. momentum component.

**Usage:**
```
root -l -q 'plot_sigmaCM.C("fitOutput.root","sigmaCM_plots.pdf")'
```

**Expected inputs in the ROOT file:**
- `sigmacmx_int`, `sigmacmy_int`, `sigmacmz_int`, `sigmacmT_int` — `TGraphAsymmErrors` for the
  integrated (Q²-summed) σ_CM, one point each
- `sigmacmx_Q2`, `sigmacmy_Q2`, `sigmacmz_Q2`, `sigmacmT_Q2` — `TGraphAsymmErrors` for σ_CM vs
  Q² bins

Produces one PDF page per component (x, y, z, T), each showing both the Q²-differential and
the integrated value.

---

### `plot_sigmaCM_overlays.C`

Overlays data histograms against an ensemble of simulated histograms from the σ_CM template bank
produced by `Main_sigmaCM.cpp`. Makes one integrated page and one page per Q² bin, useful for
visually verifying which template width best matches the data shape.

**Usage:**
```
root -l -q 'plot_sigmaCM_overlays.C("stage1.root","sigmaCM_overlays.pdf")'
```

Input file must contain:
- `pcmx_epp`, `pcmy_epp`, etc. — integrated data histograms
- `h_pcmx_epp_simSCM_<j>` (j = 0…99) — sim template histograms (integrated)
- `h_pcmx_epp_SRC_simSCM_Q2_<j>_<iq>` — sim templates per Q² bin

**Key functions:**

| Function | Description |
|---|---|
| `getHistClone(file, name, cloneName)` | Fetches and detaches a histogram from the ROOT file |
| `drawSimOverlayWithData(...)` | Overlays one data histogram with a family of sim histograms; optionally normalizes sim integrals to match data |
| `plot_sigmaCM_overlays(fileName, pdfName)` | Page 1: 2×2 integrated panels for pcm-x/y/z/T; then one 2×2 page per Q² bin |

Defaults: `nSigma = 100`, `nQ2 = 7`, `normalizeSimToData = true`. Uses ROOT palette colors to
color the template family from low (blue) to high (red) σ_CM.

---

### `plot_toy_diagnostics.C`

Inspects the toy-ensemble distributions underlying yields and the epp/ep ratio for a selected
pMiss bin and Q² bin. Useful for checking that the 100 toy variations are behaving as expected.

**Usage:**
```
root -l -b -q 'plot_toy_diagnostics.C("Data_He.root")'
root -l -b -q 'plot_toy_diagnostics.C("Data_He.root","pdf/Q2_pmiss/toy_diagnostics.pdf",3,6)'
```

Arguments: filename, optional output PDF path, optional `pMissBin` (default 3), optional
`q2ValueBin` (default 6). Bin indices are 0-based in the macro arguments but shifted +1 for ROOT
histogram access.

**Key functions:**

| Function | Description |
|---|---|
| `toyLabels(file)` | Scans the ROOT file for `hists/toy_*` directories and returns their labels |
| `drawToySummary(ax, toyValues, nominal)` | Overlays toy distributions with mean, median, ±1σ percentile band, and optional nominal value marker |
| `drawRatioToyDistributionPage(file, pMissBin, q2ValueBin, pdf)` | Computes toy-by-toy `Q2_epp_SRC_pmiss / Q2_ep_SRC_pmiss` ratio in the selected bin and histograms the 100 values |
| `drawYieldToyDistributionPage(file, taskName, pMissBin, q2ValueBin, pdf)` | Same for a single task's yield |

**Pages produced:** ep yield distribution across toys, epp yield distribution, epp/ep ratio
distribution — each for the selected pMiss + Q² bin combination.
