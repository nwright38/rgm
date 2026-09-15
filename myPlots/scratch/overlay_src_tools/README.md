# overlay_src_tools

Refactored ROOT macros for SRC overlays with shared logic in one header.

Files:
- `OverlaySrcCommon.h`: shared variables, detector definitions, styles, cut/weight helpers.
- `overlay_default_multi.C`: default detector-panel overlays for any number of files.
- `overlay_default_multi_cut_columns.C`: multi-file overlays with nominal cuts in the left column and nominal-plus-extra-cut panels in the right column.
- `overlay_data_by_detector.C`: one-file, one-variable-per-page detector overlays.
- `overlay_q2_by_detector.C`: one-file Q2-sliced overlays in detector panels.
- `overlay_with_without_cut.C`: one-file with-cut vs without-cut overlays.
- `overlay_he4_c12_data_sim_lead_ratios.C`: simple pMiss-only He4/C12 comparison macro with FD/CD lead panels, normalized He4/C12 ratios, and data-vs-sim ratio overlays.
- `convert_events2N_to_srcTree.C`: converts `events_2N.root` (`events` tree) into a `srcTree` compatible with the overlay macros.

All macros are ROOT-interpreted and require no compilation.

## 1) Default overlay for >2 files

Use comma-separated lists for files, weights, and labels.

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_default_multi.C("~/data/RGM_DATA/c12_src_skim.root,~/data/RGM_DATA/c12_sim_skim.root,~/data/RGM_DATA/c12_sim_skim_100MeV_allD.root","srcTree","overlay_default_3files.pdf",true,"pCM > 0","pCM > 0 && pMiss < 1. && recP < 1.","(weight_epp),(weight_epp),(weight_epp)","Data,Sim,Sim+FSI","",true)'
```

To draw only a subset of entries, pass `maxEvents` and `firstEvent` at the end:

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_default_multi.C("~/data/RGM_DATA/events_2N_srcTree.root,~/data/RGM_DATA/events_2N_srcTree.root","srcTree","overlay_default_limited.pdf",true,"pCM > 0","pCM > 0 && pMiss < 1. && recP < 1.","(weight_epp),(weight_epp)","2N-A,2N-B","",false,200000,0)'
```

`overlay_default_multi.C` now does a split for variables that do not require recoil quantities:
- First page: e'p selection (lead-detector panels)
- Following page: e'pp selection (lead/recoil detector-combination panels)

You can control the e'p page selection with trailing args `epCut` and `epBaseCut`:

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_default_multi.C("~/data/RGM_DATA/events_2N_srcTree.root,~/data/RGM_DATA/events_2N_srcTree.root","srcTree","overlay_default_ep_epp_split.pdf",true,"pCM > 0","pMiss < 1. && recP < 1.","(weight_epp),(weight_epp)","2N-A,2N-B","",false,-1,0,"1","pMiss < 1.","(weight_ep),(weight_ep)")'
```

Notes:
- e'pp pages use `weightsCsv` (typically `weight_epp`-based).
- e'p pages use `epWeightsCsv`, which defaults to `weight_ep`-based weights.

## 1b) Multi-file cut-column overlay

Use this when you want the usual detector-combination pages, but with the nominal selection in the left column and the same selection plus one extra cut in the right column.

Default extra cut: `xB > 1.3`

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_default_multi_cut_columns.C("~/data/RGM_DATA/events_2N_srcTree.root,~/data/RGM_DATA/events_2N_srcTree.root","srcTree","pdf/overlay_default_multi_cut_columns.pdf",true,"pCM > 0","pMiss < 1. && recP < 1.","(weight_epp),(weight_epp)","2N-A,2N-B","",false,3000,0,"xB > 1.3","xB > 1.3")'
```

Notes:
- Default layout is 2 plots per column: `Lead FD Rec CD` and `Lead CD Rec CD`.
- By default, `Lead CD Rec FD` is omitted via `omitLeadCdRecFd=true`.
- To include `Lead CD Rec FD`, pass one more trailing argument: `false`.
- Left column: nominal `e'pp` selection.
- Right column: nominal `e'pp` selection plus `extraCut`.
- The macro keeps the same multi-file overlay behavior as `overlay_default_multi.C`.

## 2) Data detector overlays

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_data_by_detector.C()'
```

## 3) Q2 overlays by detector

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_q2_by_detector.C()'
```

## 4) With/without extra cut

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_with_without_cut.C("~/data/RGM_DATA/c12_src_skim.root","srcTree","overlay_cut_compare.pdf","pCM > 0","pCM > 0 && pMiss < 1. && recP < 1.","pCMz > .4","pCMz > .4","(weight_epp)","",false,false)'
```

Notes:
- Set `includeFdFd=false` to omit FD+FD.
- All macros support `weight_epp`-style weighting through weight-expression arguments.
- If fewer labels/weights are given than files in default mode, missing values are auto-filled.

## 4b) He4/C12 pMiss lead-detector ratio pages

Use this when you have four `srcTree`-style inputs: He4 data, C12 data, He4 sim, and C12 sim.

The macro writes five pMiss pages:
- one 4-panel page with He4/C12 overlays for data FD, data CD, sim FD, sim CD
- one 4-panel page with normalized He4/C12 ratios for the same panels
- one full-page overlay of data vs sim for the He4/C12 ratio in FD leads
- one full-page overlay of data vs sim for the He4/C12 ratio in CD leads
- one single-panel overlay of double ratios: (data He4/C12) / (sim He4/C12) with FD and CD on the same axes

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_he4_c12_data_sim_lead_ratios.C("~/data/RGM_DATA/he4_src_skim.root","~/data/RGM_DATA/c12_src_skim.root","~/data/RGM_DATA/he4_sim_skim.root","~/data/RGM_DATA/c12_sim_skim.root","srcTree","pdf/he4_c12_pmiss_data_sim_lead_ratios.pdf",true,"goodLead","pMiss < 1.","(weight_ep),(weight_ep)","(weight_ep),(weight_ep)","")'
```

Notes:
- The overlay page respects `normalizeOverlayPage`; the ratio pages always divide He4/C12 after each histogram is normalized to unity.
- `dataWeightsCsv` should contain weights for `He4 data,C12 data`.
- `simWeightsCsv` should contain weights for `He4 sim,C12 sim`.

## 5) Convert `events_2N.root` to `srcTree`

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/convert_events2N_to_srcTree.C("~/data/RGM_DATA/events_2N.root","~/data/RGM_DATA/events_2N_srcTree.root","events",true,12,-1.0,true,true,true)'
```

Arguments:
- `inputFileName`: source file, usually `~/data/RGM_DATA/events_2N.root`
- `outputFileName`: output file containing `srcTree`
- `inputTreeName`: input tree name (`events`)
- `useFSIAwareMomenta`: if `true`, use post-FSI momenta when `doFSI!=0` and pre-FSI otherwise
- `targetA`: target mass number used for missing-energy/light-cone derived quantities (default `12`)
- `eBeamOverride`: optional beam-energy override in GeV; use `<0` to keep values from the input tree
- `requireEpp`: if `true` (default), require both `lead_type` and `rec_type` to be protons (`2212`)
- `applyMCSmearing`: if `true` (default), apply MC momentum smearing to electron/lead/recoil using the same FD/CD resolution functions used in `simpleSRCSkim_archive.cpp`

Basic SRC cuts are now applied by default in the converter (matching the skim-style baseline):
- `xB >= 1.2`
- `Q2 >= 1.5`
- `leadP >= 1.0 GeV/c`
- `recP >= 0.3 GeV/c`
- `0.65 <= mMiss <= 1.1 GeV`
- `0.3 <= kMiss <= 1.0 GeV/c`

You can disable or tune them with additional arguments:

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/convert_events2N_to_srcTree.C("~/data/RGM_DATA/events_2N.root","~/data/RGM_DATA/events_2N_srcTree_nocuts.root","events",true,12,-1.0,true,true,false)'
```

If you want to keep basic SRC cuts but allow non-e'pp PID combinations, set `requireEpp=false` while keeping `applyBasicSrcCuts=true`:

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/convert_events2N_to_srcTree.C("~/data/RGM_DATA/events_2N.root","~/data/RGM_DATA/events_2N_srcTree_noepp.root","events",true,12,-1.0,false,true,true)'
```
