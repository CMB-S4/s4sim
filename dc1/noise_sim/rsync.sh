#!/bin/bash

outroot=/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim

for cleared_logfile in `find cleared_logs -name '*.log'`; do
    # Use the following test to rsync a particular frequency
    [[ ! $cleared_logfile == *f150*log ]] && continue
    indir=${cleared_logfile/cleared_logs/outputs}
    indir=${indir/.log/}
    inroot=`dirname ${indir}`
    inroot=${inroot/outputs/outputs_rk}
    outdir=${outroot}/${inroot}
    echo "Staging ${indir} to ${outdir}"
    mkdir -p $outdir
    time rsync -avrP $indir $outdir/
    exit_value=$?
    if [ $exit_value -ne 0 ]; then
        echo "Failed to sync $indir!"
        exit
    fi
    echo "Success! Removing $indir"
    rm -rf $indir
    staged_logfile=${cleared_logfile/cleared_/staged_}
    staged_logdir=`dirname ${staged_logfile}`
    echo "Moving $cleared_logfile $staged_logfile"
    mkdir -p $staged_logdir
    mv $cleared_logfile $staged_logdir
done

exit

rsync -avrP /pscratch/sd/k/keskital/s4sim/dc1/cmb_sim_new/outputs/LAT0_CHLAT /global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs_rk/
rsync -avrP /pscratch/sd/k/keskital/s4sim/dc1/cmb_sim_new/outputs/coadd /global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs_rk/
#rsync -avrP /pscratch/sd/k/keskital/s4sim/dc1/foreground_sim/outputs_lowcomplexity/LAT0_CHLAT /global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim/outputs_lowcomplexity_rk/
#rsync -avrP /pscratch/sd/k/keskital/s4sim/dc1/foreground_sim/outputs_highcomplexity/LAT0_CHLAT /global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim/outputs_highcomplexity_rk/
