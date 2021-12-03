#!/bin/bash

indir_root=/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs

for telescope in LAT0_CHLAT; do
    outdir=outputs/coadd/${telescope}
    for band_dir in ${indir_root}/${telescope}/f???; do
        band=`basename $band_dir`
        echo ${telescope} ${band}GHz
        indir=${indir_root}/${telescope}/${band}
        mkdir -p $outdir
        outroot=$outdir/coadd_${telescope}_${band}
        toast_healpix_coadd \
            --outmap ${outroot}_map.fits \
            --rcond ${outroot}_rcond.fits \
            --rcond_limit 1e-3 \
            --invcov ${outroot}_invcov.fits \
            --cov ${outroot}_cov.fits \
            ${indir}/${telescope}_split_schedule_000?/mapmaker_*_noiseweighted_map.h5
    done
done
