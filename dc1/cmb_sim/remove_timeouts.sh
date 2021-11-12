#!/bin/bash

fnames=`grep -lR --include '*.log' CANCELLED logs`
for fname in ${fnames[*]}; do
    echo "Removing $fname"
    rm $fname
done
