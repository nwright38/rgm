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
NEVENTS=${2:-100000}  #now overridable as 2nd script arg; 10x NEVENTS_BKG (10% bkg-merged)
NEVENTS_BKG=$(( NEVENTS )) #10% of events are bg-merged; each bg file has max 10000 events available
TORUS=-1.0 #-1.0 for inbending(6,4 GeV) 0.5 for outbending (2 Gev)
TARGET=4-foil #Targets: liquid, 4-foil, 1-foil, Ar, Ca
SIGMACM=0.200 #GeV/c

# sanity check: bg-merger can't hand you more bg events than exist per file
if [ "$NEVENTS_BKG" -gt 10000 ]; then
    echo "ERROR: NEVENTS_BKG ($NEVENTS_BKG) exceeds max available per bg file (10000)."
    echo "Reduce NEVENTS or keep the 10:1 ratio at or below 100000:10000."
    exit 1
fi

SLURM_ARRAY_TASK_ID=$1

# Only 100 bg-merge files exist, so cycle through them regardless of how many
# array jobs / how many events per job you're running.
BKG_ID=$(( (SLURM_ARRAY_TASK_ID - 1) % 100 + 1 ))
PADDED_ID=$(printf "%05d" ${BKG_ID})

BKGDFILE=/cache/clas12/rg-m/production/bkgfiles/tor-1.00_sol-1.00/Cx4_5986MeV/c_${PADDED_ID}.hipo

#DON'T NEED TO TOUCH BELOW HERE UNLESS YOU NEED TO
ROOTOUT=${OUTPATH}/rootfiles/Jul10
LUNDOUT=${OUTPATH}/lundfiles/Jul10
MCOUT=${OUTPATH}/mchipo/Jul10
RECONOUT=${OUTPATH}/reconhipo/Jul10

# SLURM_ARRAY_TASK_ID=1 #for testing, comment out for array job submission

#GCF Generator
$GCF_EX $Z $N $BEAM_E $ROOTOUT/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.root $NEVENTS -s $SIGMACM -P $PHASE_SPACE -v -C -O

# #TO LUND File
root -b -q "${LUND_EX}(\"${ROOTOUT}/${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.root\",\"${LUNDOUT}/lund_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.txt\",\"${TARGET}\")"

# #SUBMIT GEMC MC
gemc -USE_GUI=0  -SCALE_FIELD="binary_torus, $TORUS" -SCALE_FIELD="binary_solenoid, -1.0" -N=$NEVENTS -INPUT_GEN_FILE="lund, ${LUNDOUT}/lund_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.txt" -OUTPUT="hipo, ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus$TORUS.hipo" $GCARD

bg-merger -d "ALL" -n $NEVENTS_BKG -b $BKGDFILE -i ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}.hipo -o ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}_bkg.hipo

# #RECONSTRUCTION
recon-util -y $YAML -n $NEVENTS -i ${MCOUT}/mc_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}_bkg.hipo -o ${RECONOUT}/recon_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}_torus${TORUS}_bkg.hipo