import dateutil
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

import toast.weather
import toast.utils


realization = 0

fname_in = "schedules/pole_schedule_sat.txt"

first_time = dateutil.parser.parse(f"2027-01-01 00:00:00 +0000")

# break_start = dateutil.parser.parse(f"2027-01-01 00:00:00 +0000")
# break_stop = dateutil.parser.parse(f"2027-04-01 00:00:00 +0000")

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
        # if start_time > break_start and stop_time < break_stop:
        #     continue
        mid_time = start_time + 0.5 * (stop_time - start_time)
        if weather is None:
            weather = toast.weather.SimWeather(
                time=mid_time,
                name="south_pole",
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

ntotal = 0
n2mm = 0
n3mm = 0
for pwv, time, line in zip(pwvs, times, lines):
    ntotal += 1
    if pwv <= 2:
        n2mm += 1
    if pwv <= 3:
        n3mm += 1

print(f"Total obs : {ntotal}")
print(f"PWV < 2mm : {n2mm} (f_weather = {n2mm / ntotal:.3f})")
print(f"PWV < 3mm : {n3mm} (f_weather = {n3mm / ntotal:.3f})")
