#!/bin/bash

# Generate TOAST focalplane files

# split3 = split4 = split5

cd ../.. && git checkout Phase2_concepts && pip install . --prefix=$PREFIX
mkdir -p temp/focalplanes_split3
cd temp/focalplanes_split3
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube

# split6

cd ../.. && git checkout split6 && pip install . --prefix=$PREFIX
mkdir -p temp/focalplanes_split6
cd temp/focalplanes_split6
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube

# split7

cd ../.. && git checkout split7 && pip install . --prefix=$PREFIX
mkdir -p temp/focalplanes_split7
cd temp/focalplanes_split7
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube

# split8 = split9

cd ../.. && git checkout split89 && pip install . --prefix=$PREFIX
mkdir -p temp/focalplanes_split8
cd temp/focalplanes_split8
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube

# split0 = split1 = split2
# Run this last so we end up back in the master branch

cd ../.. && git checkout master && pip install . --prefix=$PREFIX
mkdir -p temp/focalplanes_split0
cd temp/focalplanes_split0
s4_hardware_to_toast3.py --telescope LAT0
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube
cd ../../chile_optimization_sims/phase2_prep
rsync -avrP ../../temp/focalplanes_split? .
