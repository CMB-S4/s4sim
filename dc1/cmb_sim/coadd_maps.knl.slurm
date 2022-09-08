#!/bin/bash
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=128
#SBATCH --job-name=coadd_maps
#SBATCH --licenses=SCRATCH
#SBATCH --constraint=knl
#SBATCH --core-spec=4
#SBATCH --account=mp107

ulimit -c unlimited
export MALLOC_MMAP_THRESHOLD_=131072
export PYTHONSTARTUP=""
export PYTHONNOUSERSITE=1
export HOME=$SCRATCH
export OMP_NUM_THREADS=32
export HDF5_USE_FILE_LOCKING=FALSE

let nnode=$SLURM_JOB_NUM_NODES
let ntask_node=64/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=4*$OMP_NUM_THREADS

echo "Running with"
echo "            nnode = ${nnode}"
echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
echo "       ntask_node = ${ntask_node}"
echo "            ntask = ${ntask}"

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim
indir=${indir_root}/outputs

for telescope in chlat; do
    outdir=outputs/coadd/${telescope}
    mkdir -p $outdir
    #for band in f090; do
    for band in f280; do
        echo ${telescope} ${band}GHz
        logfile=coadd_${telescope}_${band}.log
        let ntotal=0
        let nfail=0
        let nfound=0
        fname_maps="coadd_maps_${band}.txt"
        rm -f $fname_maps
        for schedule in ../noise_sim/split_schedules_1_upto2mm/chlat/*txt; do
            obs=`basename --suffix=.txt $schedule`
            fname1="${indir}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
            for fname in $fname1 FAILED; do
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
        outroot=$outdir/coadd_${telescope}_${band}
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
