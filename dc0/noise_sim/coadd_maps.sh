#!/bin/bash

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim
indir1=${indir_root}/outputs_rk
indir2=${indir_root}/outputs_float32
indir3=${indir_root}/outputs

for telescope in LAT0_CHLAT; do
    outdir=outputs/coadd/${telescope}
    for band in f090; do
        input_maps=""
        let i=0
        for schedule in split_schedules_1_upto2mm/chlat/*txt; do
            obs=`basename --suffix=.txt $schedule`
            echo ${telescope} ${band}GHz
            fname1="${indir1}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
            fname2="${indir2}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
            fname3="${indir3}/$telescope/${band}/${obs}/mapmaker_${obs}_noiseweighted_map.h5"
            for fname in $fname1 $fname2 $fname3 FAILED; do
                [[ -e $fname ]] && break
            done
            if [[ $fname == FAILED ]]; then
                echo "No input map for $band $obs"
            else
                input_maps+=" $fname"
                echo "$band $obs : $fname"
            fi
            let i++
            [[ $i == 10 ]] && break
        done
        echo $input_maps
        mkdir -p $outdir
        outroot=$outdir/coadd_${telescope}_${band}
        toast_healpix_coadd \
            --outmap ${outroot}_map.fits \
            --rcond ${outroot}_rcond.fits \
            --rcond_limit 1e-3 \
            --invcov ${outroot}_invcov.fits \
            --cov ${outroot}_cov.fits \
            $input_maps
    done
done
