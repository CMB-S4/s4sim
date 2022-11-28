#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=04:00:00
#SBATCH --nodes=16
#SBATCH --job-name=DC1_foregrounds
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
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 258 hardware threads
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim
indir1=${indir_root}/outputs_rk
indir2=${indir_root}/outputs_DUMMY1
indir3=${indir_root}/outputs_DUMMY2
let nsplit=2

for telescope in LAT0_CHLAT; do
    outdir=outputs/coadd/${telescope}
    for band in f030; do
    #for band in f030 f040 f090 f150 f220 f280; do
        for isplit in `seq 1 $nsplit`; do
            echo ${telescope} ${band}GHz $isplit
            if [[ $nsplit -gt 1 ]]; then
                splitroot="_${isplit}of${nsplit}"
            else
                splitroot=""
            fi
            logfile=coadd_${telescope}_${band}${splitroot}.log
            input_maps=""
            let ntotal=0
            let nfail=0
            let nfound=0
            fname_maps="coadd_maps_${band}${splitroot}.txt"
            rm -f $fname_maps
            for schedule in split_schedules_1_upto2mm/chlat/*txt; do
                if [[ $nsplit -gt 1 ]]; then
                    start=`awk '{if (NR == 4) print $1}' $schedule`
                    let testsplit=`echo $start | sed 's/-/ /g' | awk "{print \\$3 % $nsplit}"`+1
                    [[ ! $testsplit -eq $isplit ]] && continue
                fi

                obs=`basename --suffix=.txt $schedule`
                fname1="${indir1}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
                fname2="${indir2}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
                fname3="${indir3}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
                for fname in $fname1 $fname2 $fname3 FAILED; do
                    [[ -e $fname ]] && break
                done
                if [[ $fname == FAILED ]]; then
                    echo "No input map for $band $obs" >> $logfile
                    let nfail++
                else
                    echo $fname >> $fname_maps
                    let nfound++
                fi
                let ntotal++
                #[[ $ntotal == 100 ]] && break
            done
            echo "Found a total of ${nfound} / ${ntotal} maps. ${nfail} maps were missing." >> $logfile
            mkdir -p $outdir
            outroot=$outdir/coadd_${telescope}_${band}${splitroot}
            echo "Writing $logfile"
            date >> $logfile
            srun -n $ntask -c $ncore --cpu_bind=cores toast_healpix_coadd \
                 --outmap ${outroot}_map.fits \
                 --rcond ${outroot}_rcond.fits \
                 --rcond_limit 1e-3 \
                 --invcov ${outroot}_invcov.fits \
                 --cov ${outroot}_cov.fits \
                 $fname_maps \
                 >> $logfile 2>&1
            date >> $logfile
        done
    done
done
