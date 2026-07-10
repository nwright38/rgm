#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="$HOME/data/RGM_DATA"
input_suffix="StatOnly_CD.root"

nuclei=(
#  D
  He
  C
  40Ca
  # 48Ca
  # Sn
)

sim_nucleus=(
  C_Sim
)

inputs=()
labels=()
for nucleus in "${nuclei[@]}"; do
  inputs+=("$data_dir/Data_${nucleus}_${input_suffix}")
  labels+=("$nucleus")
done

inputs+=("$data_dir/Sim_C_AV18_StatOnly_CD_new.root")
labels+=("C (Sim)")


label_args=()
for label in "${labels[@]}"; do
  label_args+=(--label "$label")
done

python "$script_dir/plot_1D_overlay.py" \
  "${inputs[@]}" \
  "${label_args[@]}" \
  --with-ratio \
  --ratio-reference-index 1 \
  --ratio-ylim 0 2 \
  --out-dir "$script_dir/pdf/1D_overlay/scratch"

# python plot_1D_overlay.py ~/data/RGM_DATA/Data_He_StatOnly_FD_pcmCut.root ~/data/RGM_DATA/Data_C_StatOnly_FD_pcmCut.root ~/data/RGM_DATA/Data_40Ca_StatOnly_FD_pcmCut.root ~/data/RGM_DATA/Sim_C_AV18_StatOnly_FD.root --label He --label C --label 40Ca --label "C (Sim)" --with-ratio --ratio-reference-index 1 --ratio-ylim 0 2 --out-dir pdf/1D_overlay/scratch