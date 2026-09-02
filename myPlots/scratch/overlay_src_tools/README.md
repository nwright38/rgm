# overlay_src_tools

Refactored ROOT macros for SRC overlays with shared logic in one header.

Files:
- `OverlaySrcCommon.h`: shared variables, detector definitions, styles, cut/weight helpers.
- `overlay_default_multi.C`: default detector-panel overlays for any number of files.
- `overlay_data_by_detector.C`: one-file, one-variable-per-page detector overlays.
- `overlay_q2_by_detector.C`: one-file Q2-sliced overlays in detector panels.
- `overlay_with_without_cut.C`: one-file with-cut vs without-cut overlays.

All macros are ROOT-interpreted and require no compilation.

## 1) Default overlay for >2 files

Use comma-separated lists for files, weights, and labels.

```bash
root -l -b -q 'myPlots/scratch/overlay_src_tools/overlay_default_multi.C("~/data/RGM_DATA/c12_src_skim.root,~/data/RGM_DATA/c12_sim_skim.root,~/data/RGM_DATA/c12_sim_skim_100MeV_allD.root","srcTree","overlay_default_3files.pdf",true,"pCM > 0","pCM > 0 && pMiss < 1. && recP < 1.","(weight_epp),(weight_epp),(weight_epp)","Data,Sim,Sim+FSI","",true)'
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
