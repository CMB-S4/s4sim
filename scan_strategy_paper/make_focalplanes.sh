#!/bin/bash

# Generate TOAST focalplane files

# Used the Phase2_concepts for this on 2026-04-24

mkdir -p focalplanes
cd focalplanes
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
# s4_hardware_to_toast3.py --telescope SAT2 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube
s4_hardware_to_toast3.py --telescope LAT0
