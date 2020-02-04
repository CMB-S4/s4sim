""" This script generates SLURM scripts to simulate input maps for the simulation tool.
"""

import os
import sys

"""
  --bands BANDS         Comma-separated list of bands: ULFL1 (20 GHz, LAT),
                        LFL1 (27 GHz LAT), LFL2 (39 GHz, LAT), LFS1 (30 GHz,
                        SAT), LFS2 (40 GHz, SAT), MFL1 (93 GHz, LAT), MFL2
                        (145 GHz, LAT), MFLS1 (85 GHz, SAT), MFLS2 (145.1 GHz,
                        SAT), MFHS1 (95 GHz, SAT), MFHS2 (155.1 GHz, SAT),
                        HFL1(225 GHz, LAT), HFL2 (278 GHz, LAT), HFS1 (220
                        GHz, SAT), HFS2 (270 GHz, SAT).Length of list must
                        equal --tubes
  --tubes TUBES         Comma-separated list of optics tubes: LT0 (HFL), LT1
                        (HFL), LT2 (HFL), LT3 (HFL), LT4 (HFL), LT5 (MFL), LT6
                        (MFL), LT7 (MFL), LT8 (MFL), LT9 (MFL), LT10 (MFL),
                        LT11 (MFL), LT12 (MFL), LT13 (MFL), LT14 (MFL), LT15
                        (MFL), LT16 (MFL), LT17 (LFL), LT18 (LFL), LT19 (HFL),
                        LT20 (HFL), LT21 (HFL), LT22 (HFL), LT23 (HFL), LT24
                        (MFL), LT25 (MFL), LT26 (MFL), LT27 (MFL), LT28 (MFL),
                        LT29 (MFL), LT30 (MFL), LT31 (MFL), LT32 (MFL), LT33
                        (MFL), LT34 (MFL), LT35 (MFL), LT36 (LFL), LT37 (LFL),
                        LT38 (HFL), LT39 (HFL), LT40 (HFL), LT41 (HFL), LT42
                        (MFL), LT43 (MFL), LT44 (MFL), LT45 (MFL), LT46 (MFL),
                        LT47 (MFL), LT48 (MFL), LT49 (MFL), LT50 (MFL), LT51
                        (MFL), LT52 (MFL), LT53 (MFL), LT54 (LFL), LT55 (LFL),
                        LT56 (ULFL), ST0 (MFLS), ST1 (MFLS), ST2 (MFLS), ST3
                        (MFLS), ST4 (MFLS), ST5 (MFLS), ST6 (MFHS), ST7
                        (MFHS), ST8 (MFHS), ST9 (MFHS), ST10 (MFHS), ST11
                        (MFHS), ST12 (HFS),ST13 (HFS), ST14 (HFS), ST15 (HFS),
                        ST16 (LFS), ST17 (LFS).Length of list must equal
                        --bands
"""

flavors = "noise", "atmosphere", "signal"

lat_tubes = [
    "LT38",  # HFL 225 & 278 GHz
    "LT42",  # MFL  93 & 145 GHz
    "LT54",  # LFL  27 &  39 GHz
    "LT56",  # ULFL,      20 GHz
]

sat_tubes = [
    "ST0",   # MFLS 85 & 145.1 GHz
    "ST6",   # MFHS 95 & 155.1 GHz
    "ST12",  # HFS 220 & 270 GHz
    "ST16",  # LFS  30 &  40 GHz
]
