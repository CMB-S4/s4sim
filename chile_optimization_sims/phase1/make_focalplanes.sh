#!/bin/bash

# Generate TOAST focalplane files

mkdir -p focalplanes
cd focalplanes
s4_hardware_to_toast3.py --telescope LAT0
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
# s4_hardware_to_toast3.py --telescope SAT2 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube
