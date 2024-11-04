#!/bin/bash

# Generate TOAST focalplane files

# split3 = split4 = split5

cd ../.. && git checkout Phase2_concepts && pip install . --prefix=$PREFIX && cd chile_optimization_sims/phase2_prep
mkdir -p focalplanes_split3
cd focalplanes_split3
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube
cd ..

# split6

cd ../.. && git checkout split6 && pip install . --prefix=$PREFIX && cd chile_optimization_sims/phase2_prep
mkdir -p focalplanes_split6
cd focalplanes_split6
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube
cd ..

# split7

cd ../.. && git checkout split7 && pip install . --prefix=$PREFIX && cd chile_optimization_sims/phase2_prep
mkdir -p focalplanes_split7
cd focalplanes_split7
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube
cd ..

# split8 = split9

cd ../.. && git checkout split89 && pip install . --prefix=$PREFIX && cd chile_optimization_sims/phase2_prep
mkdir -p focalplanes_split8
cd focalplanes_split8
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube
cd ..

# split0 = split1 = split2
# Run this last so we end up back in the master branch

cd ../.. && git checkout master && pip install . --prefix=$PREFIX && cd chile_optimization_sims/phase2_prep
mkdir -p focalplanes_split0
cd focalplanes_split0
s4_hardware_to_toast3.py --telescope LAT0
s4_hardware_to_toast3.py --telescope SAT1 --by-tube
s4_hardware_to_toast3.py --telescope SAT3 --by-tube
cd ..
