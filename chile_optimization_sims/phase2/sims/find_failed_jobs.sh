#!/bin/bash

let nmove=0
for fname_in in outputs/*/*/*/log; do
    if [[ x$(grep "Workflow completed" $fname_in) = "x" ]]; then
        for i in {1..10}; do
            fname_out=${fname_in}.$i
            if [[ ! -e $fname_out ]]; then
                echo "$fname_in -> $fname_out"
                mv $fname_in $fname_out
                let nmove++
                break
            fi
        done
    fi
done
echo "Moved $nmove logfiles"

echo "Completed jobs:"
grep -l "Workflow completed" outputs/*/*/*/log
