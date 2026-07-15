python plot_1D_overlay.py \
~/data/RGM_DATA/Data_C_StatOnly_CD.root \
~/data/RGM_DATA/Sim_C_AV18_StatOnly_CD.root \
~/data/RGM_DATA/Sim_C_N2LO10_StatOnly_CD.root \
~/data/RGM_DATA/Sim_C_N2LO12_StatOnly_CD.root \
--label CD --label CD_sim_AV18 --label CD_sim_N2LO10 --label CD_sim_N2LO12 \
--ratio-mode inputs-over-reference \
--with-ratio --ratio-reference-index 0 --out-dir pdf/1D_overlay/scratch \
--ratio-ylim 0 2 \