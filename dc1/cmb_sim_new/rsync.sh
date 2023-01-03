#!/bin/bash

rsync -avrP /pscratch/sd/k/keskital/s4sim/dc1/cmb_sim_new/outputs/LAT0_CHLAT /global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs_rk/
rsync -avrP /pscratch/sd/k/keskital/s4sim/dc1/cmb_sim_new/outputs/coadd /global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs_rk/
#rsync -avrP /pscratch/sd/k/keskital/s4sim/dc1/foreground_sim/outputs_lowcomplexity/LAT0_CHLAT /global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim/outputs_lowcomplexity_rk/
#rsync -avrP /pscratch/sd/k/keskital/s4sim/dc1/foreground_sim/outputs_highcomplexity/LAT0_CHLAT /global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim/outputs_highcomplexity_rk/
