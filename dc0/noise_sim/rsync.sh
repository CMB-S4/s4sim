#!/bin/bash

rsync -avrP /pscratch/sd/k/keskital/s4sim/dc0/noise_sim/outputs/LAT0_CHLAT /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/
exit
rsync \
    -avrP \
    --exclude='*_invcov.fits' \
    --exclude='*_cov.fits' \
    --exclude='*_rcond.fits' \
    /pscratch/sd/k/keskital/s4sim/dc0/noise_sim/outputs/coadd \
    /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/
