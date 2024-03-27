#!/bin/bash

mkdir -p maps

# SPLAT

if [[ ! -e maps/single_obs_splat_f150.h5 ]]; then
    toast_healpix_coadd \
        --outmap maps/single_obs_splat_f150.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT2_SPLAT/f150/POLE2-100-0/*noiseweighted_map.h5
fi

if [[ ! -e maps/single_pole1_splat_f150.h5 ]]; then
    toast_healpix_coadd \
        --outmap maps/single_pole1_splat_f150.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT2_SPLAT/f150/POLE1-100-*/*noiseweighted_map.h5
fi

if [[ ! -e maps/single_pole2_splat_f150.h5 ]]; then
    toast_healpix_coadd \
        --outmap maps/single_pole2_splat_f150.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT2_SPLAT/f150/POLE2-100-*/*noiseweighted_map.h5
fi

if [[ ! -e maps/single_pole3_splat_f150.h5 ]]; then
    toast_healpix_coadd \
        --outmap maps/single_pole3_splat_f150.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT2_SPLAT/f150/POLE3-100-*/*noiseweighted_map.h5
fi

# CHLAT

if [[ ! -e maps/single_obs_chlat_f150.h5 ]]; then
    toast_healpix_coadd \
        --outmap maps/single_obs_chlat_f150.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-12/*noiseweighted_map.h5
fi

if [[ ! -e maps/single_day_chlat_f150.h5 ]]; then
    toast_healpix_coadd \
        --outmap maps/single_day_chlat_f150.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-10/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-11/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-12/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-13/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-14/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-15/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-0/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-1/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-2/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-3/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-4/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-5/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-6/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-7/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-8/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-9/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-10/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-11/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-12/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-13/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-14/*noiseweighted_map.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-157-15/*noiseweighted_map.h5
fi
    
# SPSAT

if [[ ! -e maps/single_obs_spsat_f155.h5 ]]; then
    toast_healpix_coadd \
        --outmap maps/single_obs_spsat_f155.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/SAT1_SAT/f155/POLE_DEEP-100-0/*noiseweighted_filtered_map.h5
fi

if [[ ! -e maps/single_pass_spsat_f155.h5 ]]; then
    toast_healpix_coadd \
        --outmap maps/single_pass_spsat_f155.h5 \
        /global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/SAT1_SAT/f155/POLE_DEEP-100-*/*noiseweighted_filtered_map.h5
fi
