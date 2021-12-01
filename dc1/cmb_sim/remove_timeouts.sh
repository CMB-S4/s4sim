#!/bin/bash

mkdir -p failed_logs

fnames=(`grep -lR --include '*.log' CANCELLED logs`)
echo "Removing ${#fnames[@]} cancelled logs"
for fname in ${fnames[*]}; do
    echo "Removing $fname"
    mv $fname failed_logs/
done

fnames=(`grep -lR --include '*.log' Terminated logs`)
echo "Removing ${#fnames[@]} terminated logs"
for fname in ${fnames[*]}; do
    echo "Removing $fname"
    mv $fname failed_logs/
done
