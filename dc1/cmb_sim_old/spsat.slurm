#!/bin/bash
#SBATCH --partition=regular
#SBATCH --time=00:30:00
#SBATCH --nodes=32
#SBATCH --job-name=DC1_cmb
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=knl
#SBATCH --core-spec=4
#SBATCH --account=mp107

ulimit -c unlimited
export MALLOC_MMAP_THRESHOLD_=131072
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export TOAST_FUNCTIME=1
export HDF5_USE_FILE_LOCKING=FALSE

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS
# Make sure nnode is divisible by nnode_group
let nnode_group=1
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

telescope=spsat
site=pole
TELESCOPE=SAT0_SAT
#bands=(f085 f095 f145 f155 f220 f280)
bands=(f095 f145 f155 f220 f280)
#bands=(f085)

for band in ${bands[*]}; do
    logdir=logs/${telescope}/${band}
    mkdir -p $logdir
done

echo "Listing schedules at" `date`

#fnames=split_schedules/${telescope}/split_schedule_????.txt
#fnames=(`python3 get_fnames.py split_schedules/${telescope} logs/${telescope} $telescope ${bands[*]}`)
fnames=(`python3 get_fnames_serial.py split_schedules/${telescope} logs/${telescope} $telescope ${bands[*]}`)

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

echo "Looking for schedule at" `date`

let running=0

for schedule in ${fnames[*]}; do
    rootname=`basename $schedule .txt`
    rootname=${telescope}_${rootname}
    for band in ${bands[*]}; do
        logdir=logs/${telescope}/${band}
        logfile=$logdir/${rootname}_${band}.log
        if [[ ! -e $logfile ]]; then
            echo "Writing $logfile at" `date`
            date > ${logfile}

            outdir=outputs/${telescope}/${band}/${rootname}
            mkdir -p $outdir
            chgrp -R cmbs4 $outdir

            srun -N $nnode_group -n $ntask_group -c $ncore --cpu_bind=cores toast_sim_ground.py \
                 --config common.toml scanning_${telescope}.toml atmosphere_${site}.toml reduce_${telescope}.toml \
                 --focalplane ../focalplanes/focalplane_${TELESCOPE}_${band}.h5 \
                 --schedule $schedule \
                 --out $outdir \
                 --job_group_size ${groupsize} \
                 --obsmaps \
                 --mem_count.enable \
                 --scan_map.file input_maps/cmb.${telescope}.${band}.h5 \
                 --sim_noise.disable \
                 --sim_atmosphere.disable \
                 --raw_statistics.disable \
                 --filtered_statistics.disable \
                 --mapmaker.report_memory \
                 --save_hdf5.disable \
                 >> ${logfile} 2>&1 &
            let running++
            [[ $running -eq $ngroup ]] && break
            #date >> ${logfile}
            #echo "Done with $logfile at" `date`
        else
            echo "$logfile exists"
        fi
    done
    [[ $running -eq $ngroup ]] && break
done

echo "Waiting for $running jobs to complete at" `date`

wait

echo "Jobs completed at" `date`
