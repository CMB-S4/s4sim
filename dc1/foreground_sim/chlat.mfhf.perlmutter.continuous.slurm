#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --nodes=16
#SBATCH --job-name=DC1_mfhf_foregrounds
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107

# Perlmutter-specific fixes
export FI_CXI_OPTIMIZED_MRS="false"
export MPI4PY_RC_RECV_MPROBE="False"

# Python environment
ulimit -c unlimited
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export HDF5_USE_FILE_LOCKING=FALSE

# TOAST variables
export TOAST_FUNCTIME=1

# Parallelization
export OMP_NUM_THREADS=2
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 258 hardware threads
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS
# Make sure nnode is divisible by nnode_group
let nnode_group=4
let ntask_group=$nnode_group*$ntask_node
let groupsize=nnode_group*ntask_node
let ngroup=$nnode/$nnode_group

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "      nnode_group = ${nnode_group}"
echo "      ntask_group = ${ntask_group}"
echo "        groupsize = ${groupsize}"
echo "           ngroup = ${ngroup}"

telescope=chlat
site=chile
TELESCOPE=LAT0_CHLAT
# cannot have more than 384 MPI tasks in a group to process 30 and 40GHz
#bands=( f030 )
#bands=(f030 f040)
#bands=(f090 f150)
bands=(f090 f150 f220 f280)

echo "Listing schedules at" `date`

#fnames=split_schedules_1/${telescope}/split_schedule_????.txt
#fnames=(`python3 get_fnames.py split_schedules_1_upto2mm/${telescope} logs/${TELESCOPE} ${TELESCOPE} ${bands[*]}`)
fnames=(`python3 get_fnames_serial.py split_schedules_1_upto2mm/${telescope} logs/${TELESCOPE} ${TELESCOPE} ${bands[*]}`)
echo "Found ${#fnames[@]} schedule files"

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

echo "Looking for schedule at" `date`

let running=0
for band in ${bands[*]}; do
    logdir=logs/${TELESCOPE}/${band}
    mkdir -p $logdir
done

for schedule in ${fnames[*]}; do
    rootname=`basename $schedule .txt`
    #rootname=${TELESCOPE}_${rootname}
    for band in ${bands[*]}; do
        logdir=logs/${TELESCOPE}/${band}
        #logfile=$logdir/${rootname}_${band}.log
        logfile=$logdir/${rootname}.log
        logfile2=cleared_$logdir/${rootname}.log
        if [[ ! -e $logfile && ! -e $logfile2 ]]; then
            echo "Writing $logfile at" `date`
            date > ${logfile}

            outdir=outputs/${TELESCOPE}/${band}/${rootname}
            mkdir -p --mode=g+rXs-w $outdir
            chgrp -R cmbs4 $outdir

            srun -N $nnode_group -n $ntask_group -c $ncore --cpu_bind=cores toast_sim_ground.py \
                 --config common.toml scanning_${telescope}.toml \
                 atmosphere_${site}.toml reduce_${telescope}.toml \
                 --focalplane ../focalplanes/focalplane_${TELESCOPE}_${band}.h5 \
                 --telescope $TELESCOPE \
                 --schedule $schedule \
                 --out $outdir \
                 --job_group_size ${groupsize} \
                 --obsmaps \
                 --mem_count.enable \
                 --scan_map.enable \
                 --scan_map.file input_maps/foreground.${telescope}.${band}.h5 \
                 --sim_noise.disable \
                 --sim_atmosphere.disable \
                 --raw_statistics.disable \
                 --filtered_statistics.disable \
                 --mapmaker.report_memory \
                 --save_hdf5.disable \
                 >> ${logfile} 2>&1 &
            sleep 1
            let running=`ps | grep srun | wc -l`/2
            echo "There are $running running processes at" `date`
            while [[ $running -ge $ngroup ]]; do
                echo "Waiting for 10 seconds."
                sleep 10
                let running=`ps | grep srun | wc -l`/2
                echo "There are $running running processes at" `date`
            done
            echo "Finding more work at" `date`
        else
            echo "$logfile exists"
        fi
    done
done

echo "Waiting for $running jobs to complete at" `date`

wait

echo "Jobs completed at" `date`
