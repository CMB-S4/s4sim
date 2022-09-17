#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --job-name=AoA
##SBATCH --licenses=SCRATCH
#SBATCH --constraint=cpu
#SBATCH --account=mp107

ulimit -c unlimited

# This script is for the Perlmutter CPU nodes

# Perlmutter-specific fixes
export FI_CXI_OPTIMIZED_MRS="false"
export MPI4PY_RC_RECV_MPROBE="False"

# Python environment
export MALLOC_MMAP_THRESHOLD_=131072
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
#export HOME=$SCRATCH

# TOAST variables
export TOAST_FUNCTIME=1

# Parallelization
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 258 hardware threads
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS
#let groupsize=1*$ntask_node
let groupsize=1

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "        groupsize = ${groupsize}"

for flavor in spsat chsat_s4 chsat_so; do
#for flavor in chsat_so; do
    case $flavor in
	spsat)
	    #params=params.spsat.madam.toml
	    params=params.spsat.toml
	    ;;
	*)
	    #params=params.chsat.madam.toml
	    params=params.chsat.toml
	    ;;
    esac
    for band in f030 f040 f085 f095 f145 f155 f220 f280; do
    #for band in f030 f040 f085 f095 f145 f155; do
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
        exit
        # remove the useless CFITSIO warnings
        # grep -v CFITSIO $logfile | egrep -v '$^' > temp.log
        # mv temp.log $logfile
    done
done
