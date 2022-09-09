#!/bin/bash
#SBATCH --partition=regular
#SBATCH --qos=premium
#SBATCH --time=02:00:00
#SBATCH --nodes=32
#SBATCH --job-name=coadd_maps
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=haswell
#SBATCH --account=mp107

ulimit -c unlimited
export MALLOC_MMAP_THRESHOLD_=131072
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export OMP_NUM_THREADS=16
export HDF5_USE_FILE_LOCKING=FALSE

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=32/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=2*$OMP_NUM_THREADS

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim
indir1=${indir_root}/outputs_rk
indir2=${indir_root}/outputs_float32
indir3=${indir_root}/outputs
let nsplit=2

for telescope in LAT0_CHLAT; do
    outdir=outputs/coadd/${telescope}
    for band in f030 f040; do
    #for band in f030 f044 f090 f150 f220 f280; do
    #for band in f090 f150 f220 f280; do
        for isplit in `seq 0 $nsplit`; do
            echo ${telescope} ${band}GHz $isplit
            if [[ $isplit -eq 0 ]]; then
                splitroot=""
            else
                splitroot="_${isplit}of${nsplit}"
            fi
            logfile=coadd_${telescope}_${band}${splitroot}.log
            input_maps=""
            let ntotal=0
            let nfail=0
            let nfound=0
            fname_maps="coadd_maps_${band}${splitroot}.txt"
            rm -f $fname_maps
            for schedule in split_schedules_1_upto2mm/chlat/*txt; do
                if [[ $isplit -gt 0 ]]; then
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
