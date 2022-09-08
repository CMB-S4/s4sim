#!/bin/bash

fknee="--fknee 1e-6 --fmin 1e-10"

s4_hardware_to_toast3.py $fknee --telescope CHSAT0
s4_hardware_to_toast3.py $fknee --telescope CHSAT1
s4_hardware_to_toast3.py $fknee --telescope CHSAT2
s4_hardware_to_toast3.py $fknee --telescope CHSAT3
s4_hardware_to_toast3.py $fknee --telescope CHSAT4
s4_hardware_to_toast3.py $fknee --telescope CHSAT5

s4_hardware_to_toast3.py $fknee --telescope SAT0
s4_hardware_to_toast3.py $fknee --telescope SAT1
s4_hardware_to_toast3.py $fknee --telescope SAT2
s4_hardware_to_toast3.py $fknee --telescope SAT3
s4_hardware_to_toast3.py $fknee --telescope SAT4
s4_hardware_to_toast3.py $fknee --telescope SAT5

s4_hardware_to_toast3.py $fknee --telescope LAT0
s4_hardware_to_toast3.py $fknee --telescope LAT2
