#!/bin/bash

#schedule="schedules/schedule.trimmed.txt"
schedule="test_schedule.txt"
head -n 12 schedules/pole_schedule_lat.txt > $schedule
outdir="hitmaps_track_az"

export OMP_NUM_THREADS=1
ntask=4

TELESCOPE="LAT2_SPLAT"
band=f150

# Simulate a hit map using a decimated focalplane,
# low sampling rate and slow scan rate

fsample="1.0"
thinfp="64"

mkdir -p $outdir

mpirun -np ${ntask} toast_sim_ground.py \
       `# --config ../../params/scanning_splat.toml` \
       --config scanning_splat.toml \
       --focalplane ../../focalplanes/focalplane_${TELESCOPE}_${band}.h5 \
       --telescope $TELESCOPE \
       --schedule $schedule \
       --sample_rate ${fsample} \
       --thinfp ${thinfp} \
       --out_dir ${outdir} \
       --job_group_size 1 \
       --sim_noise.disable \
       --baselines.disable \
       --sim_atmosphere.disable \
       --pixels_healpix_radec.nside 64 \
       --sim_ground.track_azimuth \
    | tee "${outdir}/log"

