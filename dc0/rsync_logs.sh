#!/bin/bash

for flavor in noise_sim multimap_sim; do
    rsync -avrP /pscratch/sd/k/keskital/s4sim/dc0/$flavor/cleared_logs /global/cfs/cdirs/cmbs4/dc/dc0/staging/$flavor/outputs_rk/
    rsync -avrP /pscratch/sd/k/keskital/s4sim/dc0/$flavor/failed_logs /global/cfs/cdirs/cmbs4/dc/dc0/staging/$flavor/outputs_rk/
done
