#!/bin/bash

mkdir -p failed_logs

fnames=(`grep -lR --include '*.log' -e CANCELLED -e Terminated -e myquota -e client_io_handler_start logs`)
echo "Removing ${#fnames[@]} cancelled logs"
for fname in ${fnames[*]}; do
    echo "Removing $fname"
    mv $fname failed_logs/
done
