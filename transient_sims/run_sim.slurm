#!/bin/bash
#SBATCH --qos=premium
#SBATCH --time=00:30:00
#SBATCH --nodes=8
#SBATCH --job-name=CMBS4_transient_sim
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107a

# Python environment
ulimit -c unlimited
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export HDF5_USE_FILE_LOCKING=FALSE

# TOAST variables
export TOAST_FUNCTIME=1

# Parallelization
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 258 hardware threads
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"

telescope=chlat
site=chile
TELESCOPE=LAT0_CHLAT
bands=(f150)

fnames=split_schedule/*txt
echo "Found ${#fnames[@]} schedule files"

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

echo "Looking for schedule at" `date`

for band in ${bands[*]}; do
    for schedule in ${fnames[*]}; do
        rootname=`basename $schedule .txt`
        outdir=outputs/${TELESCOPE}/${band}/${rootname}
        mkdir -p $outdir
        logfile=$outdir/log
        if [[ ! -e $logfile ]]; then
            echo "Writing $logfile at" `date`
            date > ${logfile}

            srun -N $nnode -n $ntask -c $ncore --cpu_bind=cores \
                 toast_sim_ground.py \
                 --config config.toml \
                 --focalplane focalplane_${TELESCOPE}_${band}.h5 \
                 --telescope $TELESCOPE \
                 --schedule $schedule \
                 --out $outdir \
                 --job_group_size ${ntask} \
                 --sim_catalog.enable \
                 --sim_catalog.catalog_file "catalog.txt" \
                 >> ${logfile} 2>&1
            exit
        else
            echo "$logfile exists"
        fi
    done
done
