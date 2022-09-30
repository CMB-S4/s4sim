#!/bin/bash
#SBATCH --partition=regular
#SBATCH --qos=premium
#SBATCH --time=02:00:00
#SBATCH --nodes=64
#SBATCH --job-name=AoAsims
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=knl
#SBATCH --core-spec=4
#SBATCH --account=mp107

ulimit -c unlimited
export MALLOC_MMAP_THRESHOLD_=131072
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export TOAST_FUNCTIME=1
export HDF5_USE_FILE_LOCKING=FALSE
export MKL_FAST_MEMORY_LIMIT=0

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

for flavor in spsat chsat_s4 chsat_so; do
#for flavor in chsat_so; do
    case $flavor in
	spsat)
	    params=params.spsat.toml
	    ;;
	*)
	    params=params.chsat.toml
	    ;;
    esac
    for band in f030 f040 f085 f095 f145 f155 f220 f280; do
    #for band in f280; do
	case $band in
	    f030|f040)
		fpfile=../focalplanes/focalplane_SAT5_SAT_${band}.h5
		pwv=3.0
		;;
	    f085|f095|f145|f155)
		fpfile=../focalplanes/focalplane_SAT0_SAT_${band}.h5
		pwv=3.0
		;;
	    f220|f280)
		fpfile=../focalplanes/focalplane_SAT0_SAT_${band}.h5
		pwv=2.0
		;;
	esac

        if [[ $flavor == chsat* ]]; then
            fpfile=${fpfile/_SAT_/_CHSAT_}
            fpfile=${fpfile/_SAT/_CHSAT}
        fi

	outdir=outputs_full/sat/$flavor/$band
	mkdir -p $outdir

	if [[ -e $outdir/mapmaker_cov.fits || -e $outdir/madam_wcov.fits ]]; then
	    echo "$flavor $band already complete! Skipping..."
	    continue
	fi

	logfile=log.full.${flavor}.${band}

	echo "Writing to $logfile"
        srun -N $nnode -n $ntask -c $ncore --cpu_bind=cores toast_sim_ground.py \
	     --config $params \
	     --focalplane $fpfile \
	     --pwv_limit $pwv \
	     --telescope SAT \
	     --schedule ../scan_strategies/${flavor}/schedules/schedule_sat.txt \
	     --sample_rate 0.1 \
	     --pixels_radec.nside 128 \
	     --out $outdir \
	     --job_group_size $groupsize \
	     >& $logfile
        # remove the useless CFITSIO warnings
        grep -v CFITSIO $logfile | egrep -v '$^' > temp.log
        mv temp.log $logfile
        exit
    done
done
