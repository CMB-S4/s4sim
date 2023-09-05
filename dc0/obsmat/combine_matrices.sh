#!/bin/bash

band=f030
indirs=outputs/*/f030/*
obsmats=""
for indir in $indirs; do
    obs=$(basename $indir)
    obsmat=$indir/filterbin_${obs}_noiseweighted_obs_matrix.npz
    obsmats+=" $obsmat"
    if [[ -e $obsmat ]]; then
        echo "Found: $obsmat"
        continue
    fi
    logfile=${obsmat/.npz/.log}
    echo "Writing ${logfile}"
    toast_obsmatrix_combine ${obsmat/.npz/} >& ${logfile}
done

echo "$obsmats"
outdir="obsmats/full/${band}"
mkdir -p $outdir
obsmat_full="$outdir/obsmat_${band}.npz"
invcov_full="$outdir/invcov_${band}.fits"
cov_full="$outdir/cov_${band}.fits"

logfile=${obsmat_full/.npz/.log}
echo "Writing ${logfile}"
toast_obsmatrix_coadd \
    --outmatrix $obsmat_full \
    --invcov $invcov_full \
    --cov $cov_full \
    --rcond_limit 1e-3 \
    --double_precision \
    $obsmats \
    >& ${logfile}

#usage: toast_obsmatrix_coadd [-h] [--outmatrix OUTMATRIX] [--invcov INVCOV] [--cov COV] [--nside_submap NSIDE_SUBMAP] [--rcond_limit RCOND_LIMIT]
#                             [--double_precision]
#                             inmatrix [inmatrix ...]
