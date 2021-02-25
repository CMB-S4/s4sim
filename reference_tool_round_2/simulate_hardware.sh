#!/bin/bash

s4_hardware_sim --out hardware_CMBS4 --overwrite

s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT0 --out hardware_LAT_HFL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT0 --out hardware_LAT_MFL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT0 --out hardware_LAT_LFL --overwrite

pickle_hardware.py hardware_LAT_???_*gz

s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_HFPL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_MFPL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_LFPL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_ULFPL --overwrite

pickle_hardware.py hardware_LAT_*P*_*gz

s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes SAT4 --tubes ST14 --out hardware_SAT_LFS --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes SAT0 --tubes ST0 --out hardware_SAT_MFLS --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes SAT0 --tubes ST1 --out hardware_SAT_MFHS --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes SAT0 --tubes ST2 --out hardware_SAT_HFS --overwrite

pickle_hardware.py hardware_SAT_*gz
