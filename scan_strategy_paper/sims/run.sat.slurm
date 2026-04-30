#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=8
#SBATCH --job-name=S4_SAT_scans
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
let groupsize=8

echo "$(date) : Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"
echo "        groupsize = ${groupsize}"

flavors=(sat_deep sat_deep_bk)
# freqs=(30 40 90 150 230 290)
freqs=(90)

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

for flavor in ${flavors[*]}; do
    case $flavor in
        "sat_deep") schedule=../sat_deep/schedule.sat.txt ;;
        "sat_deep_bk") schedule=../sat_deep/schedule.sat.bk.txt ;;
        *)
            echo "Unknown flavor: $flavor"
            exit
    esac
    # DEBUG begin
    # head -n 4 $schedule > single.txt
    # schedule=single.txt
    # DEBUG end
    for freq in ${freqs[*]}; do
        band=`printf "SAT_f%03i" $freq`
        case $freq in
            30|40)
	        fpfile=../focalplanes/focalplane_SAT3_${band}_ST8.h5
                pwv_limit=3
                ;;
            90|150)
	        fpfile=../focalplanes/focalplane_SAT1_${band}_ST0.h5
                pwv_limit=3
                ;;
            230|290)
	        fpfile=../focalplanes/focalplane_SAT1_${band}_ST2.h5
                pwv_limit=2
                ;;
            *)
                echo "$(date) : Unknown freq: $freq"
                exit -1
        esac

        outdir=outputs/$flavor/$band
        mkdir -p $outdir
        logfile=$outdir/log

        rootname=${flavor}_${band}
        
        if [[ -e $outdir/mapmaker_cov.fits ]]; then
	    echo "$(date) : $flavor $band already complete! Skipping..."
	    continue
        fi

        if [[ -e $logfile ]]; then
	    echo "$(date) : $logfile exists! Skipping..."
	    continue
        fi

        echo "$(date) : Writing to $logfile"
        srun -N $nnode -n $ntask -c $ncore --cpu_bind=cores \
             toast_sim_ground.py \
             --config params.sat.toml \
             --schedule $schedule \
             --focalplane $fpfile \
             --pwv_limit $pwv_limit \
             --sample_rate 1 \
             --out $outdir \
             --job_group_size $groupsize \
             >& $logfile
        
        echo "$(date) : Done!"

        exit

    done
done

echo "$(date) : All done!"
