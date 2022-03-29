#!/bin/bash

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim
indir1=${indir_root}/outputs
indir2=${indir_root}/outputs_float32
indir3=${indir_root}/outputs_rk

for telescope in LAT0_CHLAT; do
    outdir=outputs/coadd/${telescope}
    band=f090
    let i=0
    for schedule in split_schedules_1_upto2mm/chlat/*txt; do
        obs=`basename --suffix=.txt $schedule`
        echo ${telescope} ${band}GHz
        fname1=`find $indir1/$telescope/$band/$obs -name ${obs}_*.h5 2> /dev/null`
        fname2=`find $indir2/$telescope/$band/$obs -name ${obs}_*.h5 2> /dev/null`
        fname3=`find $indir3/$telescope/$band/$obs -name ${obs}_*.h5 2> /dev/null`
        echo "fname1=$fname1"
        echo "fname2=$fname2"
        echo "fname3=$fname3"
        let i++
        [[ $i == 10 ]] && break
        #if [[ -e $indir1/$telescope/$band/$obs/${obs}_*.h5 ]]
        #indir=${indir_root}/${telescope}/${band}
        #mkdir -p $outdir
        #outroot=$outdir/coadd_${telescope}_${band}
        #toast_healpix_coadd \
        #    --outmap ${outroot}_map.fits \
        #    --rcond ${outroot}_rcond.fits \
        #    --rcond_limit 1e-3 \
        #    --invcov ${outroot}_invcov.fits \
        #    --cov ${outroot}_cov.fits \
        #    ${indir}/${telescope}_split_schedule_000?/mapmaker_*_noiseweighted_map.h5
    done
done
