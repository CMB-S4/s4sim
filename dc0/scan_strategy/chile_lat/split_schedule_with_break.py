import dateutil
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

import toast.weather
import toast.utils

realization = 0
pwv_limit = 6.0

fname_in = "schedules/chile_schedule_lat.pruned.txt"
fname_out_upto2mm = "schedules/chile_schedule_lat.pruned.upto2mm_with_break.txt"
fname_out_upto3mm = "schedules/chile_schedule_lat.pruned.upto3mm_with_break.txt"

first_time = dateutil.parser.parse(f"2027-01-01 00:00:00 +0000")

break_start = dateutil.parser.parse(f"2027-01-01 00:00:00 +0000")
break_stop = dateutil.parser.parse(f"2027-04-01 00:00:00 +0000")

with open(fname_in, "r") as schedule_in:
    header = []
    pwvs = []
    times = []
    lengths = []
    lines = []
    read_header = True
    weather = None
    for line in schedule_in:
        if line.startswith("#"):
            header.append(line)
            continue
        if read_header:
            site_name, telescope_name, site_lat, site_lon, site_alt = line.split()
            site_uid = toast.utils.name_UID(site_name)
            header.append(line)
            read_header = False
            continue


        parts = line.split()
        name = parts[7]
        if name == "CALIBRATION_BREAK":
            continue

        start_date, start_time, stop_date, stop_time = parts[:4]
        start_time = dateutil.parser.parse(f"{start_date} {start_time} +0000")
        stop_time = dateutil.parser.parse(f"{stop_date} {stop_time} +0000")
        if start_time > break_start and stop_time < break_stop:
            continue
        mid_time = start_time + 0.5 * (stop_time - start_time)
        if weather is None:
            weather = toast.weather.SimWeather(
                time=mid_time,
                name="atacama",
                site_uid=site_uid,
                realization=realization,
            )
        else:
            weather.set(time=mid_time, realization=realization, site_uid=site_uid)

        pwvs.append(weather.pwv.to_value("mm"))
        times.append((mid_time - first_time).total_seconds())
        lengths.append((stop_time - start_time).total_seconds())
        lines.append(line)

pwvs = np.array(pwvs)
times = np.array(times)
lengths = np.array(lengths)

upto2 = open(fname_out_upto2mm, "w")
upto3 = open(fname_out_upto3mm, "w")
for fout in upto2, upto3:
    for line in header:
        fout.write(line)

for pwv, time, line in zip(pwvs, times, lines):
    if pwv > pwv_limit:
        continue
    if pwv <= 2:
        upto2.write(line)
    if pwv <= 3:
        upto3.write(line)

upto2.close()
upto3.close()
