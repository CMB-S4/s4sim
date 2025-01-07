#!/bin/bash

if [[ -e focalplane_LAT0_CHLAT_f020.h5 ]]; then
    echo "Focalplane exists"
else
    s4_hardware_to_toast3.py --telescope LAT0
fi
