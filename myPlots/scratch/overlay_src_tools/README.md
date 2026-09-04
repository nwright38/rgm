# overlay_src_tools

Refactored ROOT macros for SRC overlays with shared logic in one header.

Files:
- `OverlaySrcCommon.h`: shared variables, detector definitions, styles, cut/weight helpers.
- `overlay_default_multi.C`: default detector-panel overlays for any number of files.
- `overlay_data_by_detector.C`: one-file, one-variable-per-page detector overlays.
- `overlay_q2_by_detector.C`: one-file Q2-sliced overlays in detector panels.
- `overlay_with_without_cut.C`: one-file with-cut vs without-cut overlays.
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
