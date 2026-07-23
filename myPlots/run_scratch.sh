python plot_1D_overlay.py \
~/data/RGM_DATA/Data_C_StatOnly_CD.root \
~/data/RGM_DATA/Data_C_StatOnly_FD.root \
~/data/RGM_DATA/Sim_C_AV18_StatOnly_CD.root \
~/data/RGM_DATA/Sim_C_AV18_StatOnly_FD_new.root \
--label CD --label FD --label CD_sim_AV18 --label FD_sim_AV18 \
--out-dir pdf/1D_overlay/scratch \
  --ratio-ylim 0 2 \
  --with-ratio \
  --ratio-mode pairwise \
  --ratio-pair 0 2 \
  --ratio-pair 1 3 \

# --ratio-mode inputs-over-reference \
# --with-ratio --ratio-reference-index 0 


#--ratio-pair 0 1 means input 1 / input 2