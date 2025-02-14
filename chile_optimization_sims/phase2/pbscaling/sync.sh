#!/bin/bash


#for flavor in cmb with_pbscaling no_pbscaling; do
for flavor in no_pbscaling_no_artifact; do
    rsync -avrP $flavor /global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2/
    chmod g+rX-w -R /global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2/$flavor
    chgrp cmbs4 -R /global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2/$flavor
done
