#!/bin/bash

rsync -vP ../dc0/scan_strategy/chile_sat/variations/schedules/sample.txt .

s4_hardware_to_toast3.py --telescope SAT1 --by-tube
