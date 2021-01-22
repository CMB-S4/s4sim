#!/bin/bash

s4_hardware_sim --out hardware_CMBS4 --overwrite

s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT0 --out hardware_LAT_HFL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT0 --out hardware_LAT_MFL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT0 --out hardware_LAT_LFL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_ULFL --overwrite

s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_HFPL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_MFPL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_LFPL --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes LAT2 --out hardware_LAT_ULFPL --overwrite

s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes SAT0 --tubes ST0 --out hardware_SAT_MFLS --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes SAT2 --tubes ST6 --out hardware_SAT_MFHS --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes SAT4 --tubes ST12 --out hardware_SAT_HFS --overwrite
s4_hardware_trim --hardware hardware_CMBS4.toml.gz --telescopes SAT5 --tubes ST16 --out hardware_SAT_LFS --overwrite
