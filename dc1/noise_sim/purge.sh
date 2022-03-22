#!/bin/bash

# Remove local files that have been synced to CFS

fnames=(`find outputs -name "*.h5"`)

for fname_orig in ${fnames[*]}; do
    fname_copy=/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/$fname_orig
    if [[ -e $fname_copy ]]; then
        size_orig=`wc -c $fname_orig | awk '{print $1}'`
        size_copy=`wc -c $fname_copy | awk '{print $1}'`
        if [[ $size_orig == $size_copy ]]; then
            echo "Purging $fname_orig == $fname_copy"
            rm -f $fname_orig
        else
            echo "MISMATCH : $fname_orig != $fname_copy"
            echo "MISMATCH : ($size_orig != $size_copy)"
        fi
    fi
    #fname_copy=/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/${fname_orig/outputs/outputs_float32}
    fname_copy=/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/${fname_orig/outputs/outputs_rk}
    if [[ -e $fname_copy ]]; then
        size_orig=`wc -c $fname_orig | awk '{print $1}'`
        size_copy=`wc -c $fname_copy | awk '{print $1}'`
        if [[ $size_orig == $size_copy ]]; then
            echo "Purging $fname_orig == $fname_copy"
            rm -f $fname_orig
        else
            echo "MISMATCH : $fname_orig != $fname_copy"
            echo "MISMATCH : ($size_orig != $size_copy)"
        fi
    fi
done
