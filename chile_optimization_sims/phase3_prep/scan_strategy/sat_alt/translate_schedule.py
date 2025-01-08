import os
import sys

from astropy.time import Time
import numpy as np

from toast.coordinates import to_UTC, to_JD, to_MJD, DJDtoUNIX

import warnings
from erfa import ErfaWarning
# astropy.time.Time throws ErfaWarnings for future dates
warnings.filterwarnings("ignore", category=ErfaWarning)


if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <schedule_in> <schedule_out>")
    sys.exit()

fname_in = sys.argv[1]
fname_out = sys.argv[2]

rotation = 0.0

counts = {}

with open(fname_out, "w") as f_out:
    f_out.write("#Site           Telescope        Latitude [deg] Longitude [deg]   Elevation [m]\n")
    f_out.write("ATACAMA         SAT                     -22.958         -67.786          5200.0\n")
    f_out.write("#        Start time UTC             Stop time UTC Rotation Patch name            Az min  Az max      El Pass Sub\n")
    with open(fname_in, "r") as f_in:
        jd_offset = None
        for line in f_in:
            if line.startswith("#"):
                jd_offset = float(line.split()[0].split("-")[-1])
                continue
            fields = line.split()
            jd_start = float(fields[0])
            jd_stop = float(fields[1])
            el = float(fields[2])
            az = float(fields[3])
            flag = int(fields[4])
            throw = float(fields[5])
            if flag == 0:
                # Observation violates Sun/Moon avoidance
                continue
            d = 22.5
            if az > 360 - d or az < d:
                name = "N"
            elif az < 3 * d:
                name = "NE"
            elif az < 5 * d:
                name = "E"
            elif az < 7 * d:
                name = "SE"
            elif az < 9 * d:
                name = "S"
            elif az < 11 * d:
                name = "SW"
            elif az < 13 * d:
                name = "W"
            else:
                name = "NW"
            if name not in counts:
                counts[name] = 0
            pass_ = counts[name]
            sub = 0
            counts[name] += 1
                
            t_start = Time(jd_start + jd_offset, format="jd")
            t_stop = Time(jd_stop + jd_offset, format="jd")
            az_min = az - throw / 2
            az_max = az + throw / 2
            f_out.write(f"{t_start.iso:25} {t_stop.iso:25} {rotation:6.1f} {name:19} {az_min:8.3f} {az_max:8.3f} {el:6.4f} {pass_:4} {sub:3}\n")
