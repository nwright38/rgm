#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

inputs=(
  "$HOME/data/RGM_DATA/Data_He_StatOnly.root"
  "$HOME/data/RGM_DATA/Data_C_StatOnly.root"
  "$HOME/data/RGM_DATA/Data_40Ca_StatOnly.root"
  "$HOME/data/RGM_DATA/Data_48Ca_StatOnly.root"
  "$HOME/data/RGM_DATA/Data_Sn_StatOnly.root"
)

labels=(
  He
  C
  40Ca
  48Ca
  Sn
)

label_args=()
for label in "${labels[@]}"; do
  label_args+=(--label "$label")
done

python "$script_dir/plot_1D_overlay.py" \
  "${inputs[@]}" \
  "${label_args[@]}" \
  --with-ratio \
  --ratio-reference-index 0 \
  --ratio-ylim 0 2.5 \
  --out-dir "$script_dir/pdf/1D_overlay/All"
