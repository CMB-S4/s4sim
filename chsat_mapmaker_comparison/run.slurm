#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=2
#SBATCH --job-name=CMBS4_CHSAT_sim
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
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 258 hardware threads
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS
# Make sure nnode is divisible by nnode_group
let nnode_group=1
let ntask_group=$nnode_group*$ntask_node
#let groupsize=nnode_group*ntask_node
let groupsize=16
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

schedule="sample.txt"
# schedule=schedule.txt
# head -n 4 sample.txt > $schedule

# for flavor in noise hwpss ps atmosphere ground; do
for flavor in noise; do
    logfile=${flavor}.log
    if [[ -e $logfile ]]; then
	echo "$logfile exists. Skipping..."
	#continue
    fi

    echo "Writing to $logfile"

    outdir=outputs/$flavor
    mkdir -p $outdir

    case $flavor in
        noise)
            args="--sim_noise.enable"
            ;;
	hwpss)
	    args="--sim_hwpss.enable"
	    args+=" --sim_ground.hwp_angle 'hwp_angle'"
	    args+=" --sim_ground.hwp_rpm 120"
	    #args+=" --demodulate.enable"
	    args+=" --hwpfilter.enable"
	    args+=" --hwpfilter.filter_order 10"
	    args+=" --hwpfilter.trend_order 5"
	    ;;
	ps)
	    args="--scan_healpix_map.enable"
	    args+=" --scan_healpix_map.file sim_ps_map_150_24arcmin_celestial.fits"
	    args+=" --pixels_healpix_radec.nside 2048"
	    args+=" --pixels_healpix_radec.no_nest"
	    ;;
	atmosphere)
	    args="--sim_atmosphere.enable"
	    args+=" --sim_atmosphere.cache_dir 'atm_cache'"
	    args+=" --sim_atmosphere.field_of_view Quantity('35deg')"
	    args+=" --sim_atmosphere.lmin_center Quantity('0.001m')"
	    args+=" --sim_atmosphere.lmin_sigma Quantity('0.0001m')"
	    args+=" --sim_atmosphere.lmax_center Quantity('1.0m')"
	    args+=" --sim_atmosphere.lmax_sigma Quantity('0.1m')"
	    args+=" --sim_atmosphere.xstep Quantity('4m')"
	    args+=" --sim_atmosphere.ystep Quantity('4m')"
	    args+=" --sim_atmosphere.zstep Quantity('4m')"
	    args+=" --sim_atmosphere.zmax Quantity('200m')"
	    args+=" --sim_atmosphere.gain 4e-5"
	    args+=" --sim_atmosphere.wind_dist Quantity('1000m')"
	    ;;
	ground)
	    args="--sim_sss.enable"
	    args+=" --sim_sss.path ground_pickup_nside512_gauss45deg_plus_1overl2_cosEL2_150GHz_sm3deg.fits"
	    args+=" --sim_sss.pol"
	    ;;
	*)
	    echo "Unknown flavor: $flavor"
	    continue
	    ;;
    esac

    srun -N $nnode -n $ntask -c $ncore --cpu_bind=cores \
        toast_sim_ground.py \
        --job_group_size ${groupsize} \
        --mem_count.enable \
	--focalplane focalplane_SAT1_SAT_f155_ST0.h5 \
	--telescope SAT \
	--schedule $schedule \
	--sample_rate 100.0 \
	--pixels_healpix_radec_final.enable \
	--pixels_healpix_radec_final.nside 512 \
	--out $outdir \
	--sim_atmosphere.disable \
	--sim_noise.disable \
	--mapmaker.no_write_map \
	--mapmaker.no_write_rcond \
	--baselines.disable \
	--polyfilter1D.disable \
	${args} \
	>& $logfile

    exit

done
