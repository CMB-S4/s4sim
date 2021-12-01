#!/bin/bash

# Remove local files that have been synced to CFS

fnames=(`find /global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs -name "*.h5"`)

for fname_copy in ${fnames[*]}; do
    echo $fname_copy
    fname_orig=${fname_copy/\/global\/cfs\/cdirs\/cmbs4\/dc\/dc1\/staging\/cmb_sim\/}
    if [[ -e $fname_orig ]]; then
        echo "Deleting $fname_orig"
        rm -f $fname_orig
    fi
done
