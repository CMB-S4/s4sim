#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=8
#SBATCH --job-name=S4_LAT_scans
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

flavors=(lat_deep lat_deep_bk lat_wide)
# freqs=(20 30 40 90 150 220 280)
freqs=(90)

# Random wait time to reduce clashes
sleep $((RANDOM % 15))

for flavor in ${flavors[*]}; do
    case $flavor in
        "lat_deep")
            schedule=../lat_deep/schedule_lat_deep.txt
            params=params.lat.toml
            ;;
        "lat_deep_bk")
            schedule=../lat_deep/schedule_lat_deep.bk.txt
            params=params.lat.toml
            ;;
        "lat_wide")
            schedule=../lat_wide/schedule_lat_wide.txt
            params=params.lat.wide.toml
            ;;
        *)
            echo "Unknown flavor: $flavor"
            exit
    esac
    # DEBUG begin
    # head -n 4 $schedule > single.txt
    # schedule=single.txt
    # DEBUG end
    for freq in ${freqs[*]}; do
        band=`printf "f%03i" $freq`
        fpfile=../focalplanes/focalplane_LAT0_CHLAT_${band}.h5
        case $freq in
            30|40)
                pwv_limit=3
                thinfp=""
                ;;
            90|150)
                pwv_limit=3
                thinfp="--thinfp 16"
                ;;
            220|280)
                pwv_limit=2
                thinfp="--thinfp 16"
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
             --config $params \
             --schedule $schedule \
             --focalplane $fpfile \
             --pwv_limit $pwv_limit \
             --sample_rate 1 \
             --out $outdir \
             --job_group_size $groupsize \
             $thinfp \
             >& $logfile
        
        echo "$(date) : Done!"

        exit

    done
done

echo "$(date) : All done!"
