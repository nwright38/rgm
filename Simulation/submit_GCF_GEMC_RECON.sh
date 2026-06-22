#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000
#SBATCH --account=clas12
#SBATCH --job-name=mc_rgm_gcf
#SBATCH --partition=production
#SBATCH --time=10:00:00
#SBATCH --output=/farm_out/nwright/%x-%j-%N.out
#SBATCH --error=/farm_out/nwright/%x-%j-%N.err
#SBATCH --array=1-401 #Number of files 1-N

#ADD YOUR PATHS FROM INPUT TO OUTPUT HERE
#exacutables
GCF_EX=/work/clas12/users/nwright/GCF_Generator_Suite/build/programs/genQE/genQE
LUND_EX=/work/clas12/users/nwright/rgm_andrew/Simulation/GCF_to_LUND.C

#output
#change this to your specific output path
OUTPATH=/volatile/clas12/users/nwright/RGM_SIM
FILE_PREFIX=src_qe_c12_6gev

#input
#Here is where I put the configurations files I want for this simulation
INPATH=/work/clas12/users/nwright/rgm_andrew/Simulation
PHASE_SPACE=${INPATH}/phase.txt
GCARD=${INPATH}/submit/gcards/rgm_fall2021_Cx4.gcard
YAML=${INPATH}/submit/rgm_fall2021-ai_6Gev.yaml

Z=6
N=6
BEAM_E=5.98636  #5.98636, 4.02962, 2.07052
NEVENTS=1000  #10 times the next number
NEVENTS_BKG=100 #Max 10000 because this is how many background merging events there are per file
TORUS=-1.0 #-1.0 for inbending(6,4 GeV) 0.5 for outbending (2 Gev)
TARGET=4-foil #Targets: liquid, 4-foil, 1-foil, Ar, Ca
SIGMACM=0.200 #GeV/c


SLURM_ARRAY_TASK_ID=$1
PADDED_ID=$(printf "%05d" ${SLURM_ARRAY_TASK_ID})

BKGDFILE=/cache/clas12/rg-m/production/bkgfiles/tor-1.00_sol-1.00/Cx4_5986MeV/c_${PADDED_ID}.hipo

#DON'T NEED TO TOUCH BELOW HERE UNLESS YOU NEED TO
ROOTOUT=${OUTPATH}/rootfiles
LUNDOUT=${OUTPATH}/lundfiles
MCOUT=${OUTPATH}/mchipo
RECONOUT=${OUTPATH}/reconhipo

# SLURM_ARRAY_TASK_ID=1 #for testing, comment out for array job submission

#GCF Generator
$GCF_EX $Z $N $BEAM_E $ROOTOUT/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.root $NEVENTS -s $SIGMACM -P $PHASE_SPACE -v -C -O

# #TO LUND File
root -b -q "${LUND_EX}(\"${ROOTOUT}/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.root\",\"${LUNDOUT}/lund_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.txt\",\"${TARGET}\")"

# #SUBMIT GEMC MC
# #gemc -USE_GUI=0  -SCALE_FIELD="binary_torus, $TORUS" -SCALE_FIELD="binary_solenoid, -1.0" -N=$NEVENTS -INPUT_GEN_FILE="lund, /volatile/clas12/rg-m/awild/mc/lundfiles/lund_6gev_4He_StandardGCF_00000.txt" -OUTPUT="hipo, ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo" $GCARD
gemc -USE_GUI=0  -SCALE_FIELD="binary_torus, $TORUS" -SCALE_FIELD="binary_solenoid, -1.0" -N=$NEVENTS -INPUT_GEN_FILE="lund, ${LUNDOUT}/lund_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.txt" -OUTPUT="hipo, ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo" $GCARD

bg-merger -d "ALL" -n $NEVENTS_BKG -b $BKGDFILE -i ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}.hipo -o ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}_bkg.hipo

# #RECONSTRUCTION
recon-util -y $YAML -n $NEVENTS -i ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}_bkg.hipo -o ${RECONOUT}/recon_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}_bkg.hipo
