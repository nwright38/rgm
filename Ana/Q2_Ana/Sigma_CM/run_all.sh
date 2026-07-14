

DATA_FILE='/work/clas12/users/nwright/RGM/Q2/hipo/Data_C_skim.hipo'
SIM_FILE='/work/clas12/users/nwright/RGM/Q2/hipo/Sim_C_AV18_allD.hipo'
OUT_PREFIX='/work/clas12/users/nwright/RGM/Q2/SigCM/C'


A=12

./sigmacm_make_skim $A data ${OUT_PREFIX}/data_skim.root $DATA_FILE
./sigmacm_make_skim $A mc ${OUT_PREFIX}/mc_skim.root $SIM_FILE

./sigmacm_extract --from-skim ${OUT_PREFIX}/nominal.root ${OUT_PREFIX}/data_skim.root ${OUT_PREFIX}/mc_skim.root
./sigmacm_run_cut_toys ${OUT_PREFIX}/data_skim.root ${OUT_PREFIX}/mc_skim.root ${OUT_PREFIX}/cut_toys.root --n-cut-toys=200 --n-bootstrap=200
./sigmacm_run_combined_toys ${OUT_PREFIX}/data_skim.root ${OUT_PREFIX}/mc_skim.root ${OUT_PREFIX}/combined_toys.root --n-toys=200
./sigmacm_run_fit_range_scan ${OUT_PREFIX}/data_skim.root ${OUT_PREFIX}/mc_skim.root ${OUT_PREFIX}/fit_ranges.root
./sigmacm_run_closure ${OUT_PREFIX}/mc_skim.root ${OUT_PREFIX}/closure.root