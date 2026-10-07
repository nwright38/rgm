#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="$HOME/data/RGM_DATA"

inputs=(
  "$data_dir/he4_src_stat_sys_plots.root"
  "$data_dir/c12_src_stat_sys_plots.root"
  "$data_dir/ca40_src_stat_sys_plots.root"
  "$data_dir/ca48_src_stat_sys_plots.root"
  "$data_dir/sn_src_stat_sys_plots.root"
)

labels=(
  "He4 (Data)"
  "C (Data)"
  "Ca40 (Data)"
  "Ca48 (Data)"
  "Sn (Data)"
)


label_args=()
for label in "${labels[@]}"; do
  label_args+=(--label "$label")
done

python "$script_dir/plot_1D_overlay.py" \
  "${inputs[@]}" \
  "${label_args[@]}" \
  --out-dir "$script_dir/pdf/1D_overlay/scratch"


  # python "$script_dir/plot_1D_overlay.py" \
  # "${inputs[@]}" \
  # "${label_args[@]}" \
  # --with-ratio \
  # --ratio-reference-index 0 \
  # --ratio-ylim 0 2 \
  # --out-dir "$script_dir/pdf/1D_overlay/scratch"

# python plot_1D_overlay.py ~/data/RGM_DATA/Data_He_StatOnly_FD_pcmCut.root ~/data/RGM_DATA/Data_C_StatOnly_FD_pcmCut.root ~/data/RGM_DATA/Data_40Ca_StatOnly_FD_pcmCut.root ~/data/RGM_DATA/Sim_C_AV18_StatOnly_FD.root --label He --label C --label 40Ca --label "C (Sim)" --with-ratio --ratio-reference-index 1 --ratio-ylim 0 2 --out-dir pdf/1D_overlay/scratch

# python "$script_dir/plot_Q2_pmiss_extracted.py" \
#   "$data_dir/c12_src_stat_only_plots.root" \
#   "$data_dir/c12_sim_stat_only_plots.root" \
#   --out-dir "$script_dir/pdf/1D_overlay/scratch"