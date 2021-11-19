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
export OMP_NUM_THREADS=2
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export TOAST_FUNCTIME=1

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS
# Make sure nnode is divisible by nnode_group
let nnode_group=32
let groupsize=nnode_group*ntask_node

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "      nnode_group = ${nnode_group}"
echo "        groupsize = ${groupsize}"

telescope=chlat
site=chile
TELESCOPE=LAT0_CHLAT
# cannot have more than 384 MPI tasks in a group to process 30 and 40GHz
#bands=(f030 f040)
bands=(f090 f150 f220 f280)

for band in ${bands[*]}; do
    logdir=logs/${telescope}/${band}
    mkdir -p $logdir
done

fnames=split_schedules/${telescope}/split_schedule_????.txt

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

echo "Looking for schedule at" `date`

for schedule in $fnames; do
    rootname=`basename $schedule .txt`
    rootname=${telescope}_${rootname}
    for band in ${bands[*]}; do
        logdir=logs/${telescope}/${band}
        outdir=outputs/${telescope}/${band}/${rootname}
        mkdir -p $outdir
        logfile=$logdir/${rootname}_${band}.log
        if [[ ! -e $logfile ]]; then
            echo "Writing $logfile at" `date`
            date > ${logfile}

            srun -n $ntask -c $ncore --cpu_bind=cores toast_sim_ground.py \
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
                 >> ${logfile} 2>&1
            date >> ${logfile}
            echo "Done with $logfile at" `date`
            exit

        else
            echo "$logfile exists"
        fi
    done
done
