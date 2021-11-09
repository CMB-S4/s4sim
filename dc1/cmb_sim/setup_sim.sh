#!/bin/bash

# split the observing schedules

nline=1

#for telescope in chlat splat spsat; do
for telescope in splat spsat; do
    case $telescope in
        chlat)
            schedule_in=../scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned.txt
            ;;
        splat)
            schedule_in=../scan_strategy/pole_lat/schedules/pole_schedule_lat.pruned.txt
            ;;
        spsat)
            schedule_in=../scan_strategy/pole_sat/schedules/pole_schedule_sat.pruned.txt
            ;;
    esac
        
    outdir=split_schedules/${telescope}
    mkdir -p ${outdir}
    rm -rf ${outdir}/*

    schedule_out="${outdir}/split_schedule_"
    let nces=`wc -l ${schedule_in} | awk '{print $1 - 3}'`
    echo "NCES = ${nces}"
    head -n 3 $schedule_in > header.txt
    awk '{if (NR > 3) print}' $schedule_in | \
        split --lines=${nline} --numeric-suffixes \
              --suffix-length=4 --additional-suffix=.txt - \
          $schedule_out
    for fname in ${schedule_out}*; do
        mv ${fname} temp.txt
        cat header.txt temp.txt > ${fname}
        rm temp.txt
    done
    
    rm header.txt
done

exit


# get a default parameter file. Don't worry, this produces an error message
# about missing input files.

toast_sim_ground.py \
    --defaults_toml defaults.toml \
    --focalplane XXX \
    --schedule XXX \
    >& /dev/null

# create focalplanes

toast_fake_focalplane --fov_deg 7.1 --sample_rate 10 --min_pix 20 --out dummy_LAT_focalplane.h5
toast_fake_focalplane --out dummy_focalplane.220Hz.h5 --min_pix 200 --sample_rate 220
toast_fake_focalplane --out dummy_focalplane.20Hz.h5 --min_pix 200 --sample_rate 20

# 400Hz for LAT and 100Hz for SAT
# 200Hz for LAT and 20Hz for SAT

s4_hardware_to_toast3.py --telescope LAT0 --fsample 220
s4_hardware_to_toast3.py --telescope LAT1 --fsample 220
s4_hardware_to_toast3.py --telescope LAT2 --fsample 220
s4_hardware_to_toast3.py --telescope SAT0 --fsample 20
s4_hardware_to_toast3.py --telescope SAT1 --fsample 20
s4_hardware_to_toast3.py --telescope SAT2 --fsample 20
s4_hardware_to_toast3.py --telescope SAT3 --fsample 20
s4_hardware_to_toast3.py --telescope SAT4 --fsample 20
s4_hardware_to_toast3.py --telescope SAT5 --fsample 20

exit

# Extract one day of scheduled observations

fname_in=../scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned.txt
fname_out=schedule.chlat.txt
head -n 3 $fname_in > $fname_out
grep "2027-07-01" $fname_in >> $fname_out

fname_in=../scan_strategy/pole_lat/schedules/pole_schedule_lat.pruned.txt
fname_out=schedule.splat.txt
head -n 3 $fname_in > $fname_out
grep "2027-07-01" $fname_in >> $fname_out

fname_in=../scan_strategy/pole_sat/schedules/pole_schedule_sat.pruned.txt
fname_out=schedule.spsat.txt
head -n 3 $fname_in > $fname_out
grep "2027-07-01" $fname_in >> $fname_out

# Reorder and scale the input CMB map

python -c '
import numpy as np
from healpy import read_map;
from toast.pixels_io import write_healpix
m = read_map(
    "/global/cfs/cdirs/cmbs4/dm/dstool_202102/input_pysm"
    "/4096/cmb_unlensed_solardipole/0000/"
    "cmbs4_cmb_unlensed_solardipole_uKCMB_LAT-MFL1_nside4096_0000.fits",
    None,
    nest=True,
    verbose=False,
    dtype=np.float32,
);
write_healpix(
    "cmb.chlat.f090.h5",
    m * 1e-6,
    coord="G",
    nest=True,
)
'

python -c '
import numpy as np
from healpy import read_map;
from toast.pixels_io import write_healpix
m = read_map(
    "/global/cfs/cdirs/cmbs4/dm/dstool_202102/input_pysm"
    "/4096/cmb_unlensed_solardipole/0000/"
    "cmbs4_cmb_unlensed_solardipole_uKCMB_LAT-MFPL1_nside4096_0000.fits",
    None,
    nest=True,
    verbose=False,
    dtype=np.float32,
);
write_healpix(
    "cmb.splat.f090.h5",
    m * 1e-6,
    coord="G",
    nest=True,
)
'

python -c '
import numpy as np
from healpy import read_map;
from toast.pixels_io import write_healpix
m = read_map(
    "/global/cfs/cdirs/cmbs4/dm/dstool_202102/input_pysm"
    "/512/cmb_unlensed_solardipole/0000/"
    "cmbs4_cmb_unlensed_solardipole_uKCMB_SAT-MFLS1_nside512_0000.fits",
    None,
    nest=True,
    verbose=False,
    dtype=np.float32,
);
write_healpix(
    "cmb.spsat.f085.h5",
    m * 1e-6,
    coord="G",
    nest=True,
)
'
