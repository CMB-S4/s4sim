#!/bin/bash
#SBATCH --partition=regular
#SBATCH --qos=premium
#SBATCH --time=02:00:00
#SBATCH --nodes=8
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
export OMP_NUM_THREADS=4
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
#let nnode_group=1
#let ntask_group=$nnode_group*$ntask_node
#let groupsize=nnode_group*ntask_node
#let ngroup=$nnode/$nnode_group

echo "Start: " `date`
echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
#echo "      nnode_group = ${nnode_group}"
#echo "      ntask_group = ${ntask_group}"
#echo "        groupsize = ${groupsize}"
#echo "           ngroup = ${ngroup}"

#for flavor in spsat s4_like so_like; do
#for flavor in splat chlat_so chlat_s4 splat_from_chile splat_flat; do
#for flavor in splat_for_spsat splat_for_spsat chlat_for_chsat_so chlat_for_chsat_s4 chlat_for_spsat chlat_wide chlat_for_spsat_wo_wide; do
#for flavor in splat_for_spsat chlat_for_chsat_so chlat_for_chsat_s4 chlat_for_spsat chlat_wide; do
for flavor in chlat_for_chsat_s4 chlat_for_chsat_so chlat_for_spsat chlat_for_spsat_w_wide chlat_wide splat_for_spsat; do
    case $flavor in
	splat_for_spsat)
	    bands=( f020 f030 f040 f090 f150 f220 f280 )
	    fproot=../focalplanes/focalplane_LAT2_SPLAT
	    params=params.splat.toml
	    ;;
	*)
	    bands=( f030 f040 f090 f150 f220 f280 )
	    fproot=../focalplanes/focalplane_LAT0_CHLAT
	    params=params.chlat.toml
	    ;;
    esac
    bands=( f150 )
    for band in ${bands[*]}; do
	fpfile=${fproot}_${band}.h5
	case $band in
	    f020|f030|f040)
		pwv=3.0
		thinfp=32
		;;
	    f220|f280)
		pwv=2.0
		thinfp=128
		;;
	    *)
		pwv=3.0
		thinfp=128
		;;
	esac

	outdir=outputs/lat/$flavor/$band
	mkdir -p $outdir

	if [[ -e $outdir/mapmaker_cov.fits ]]; then
	    echo "$flavor $band already complete! Skipping..."
	    continue
	fi

	logfile=log.${flavor}.${band}

	echo "Writing to $logfile"
        srun -N $nnode -n $ntask -c $ncore --cpu_bind=cores toast_sim_ground.py \
	     --config $params \
	     --focalplane $fpfile \
	     --pwv_limit $pwv \
	     --telescope LAT \
	     --schedule ../scan_strategies/${flavor}/schedules/schedule_lat.txt \
	     --sample_rate 1.00 \
	     --pixels_radec.nside 128 \
	     --thinfp $thinfp \
             --fknee 1e-6 \
	     --out $outdir \
	     --job_group_size 1 \
	     >& $logfile
    done
done
