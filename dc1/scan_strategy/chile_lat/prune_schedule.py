import dateutil
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from toast.weather import SimWeather


# Fraction of the year surviving the weather cuts
fyear = 0.6
# Down-time from 3-5 day long breaks
break_length = 4 * 86400
break_length_sigma = 1 * 86400
fdown = 0.8


fname_in = "schedules/chile_schedule_lat.txt"
fname_out = "schedules/chile_schedule_lat.pruned.txt"

first_time = None

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
        if first_time is None:
            first_time = start_time
        mid_time = start_time + 0.5 * (stop_time - start_time)
        if weather is None:
            weather = SimWeather(time=mid_time, name="atacama")
        else:
            weather.set(time=mid_time)

        pwvs.append(weather.pwv.to_value("mm"))
        times.append((mid_time - first_time).total_seconds())
        lengths.append((stop_time - start_time).total_seconds())
        lines.append(line)

pwvs = np.array(pwvs)
times = np.array(times)
lengths = np.array(lengths)

ind = np.argsort(pwvs)
sorted_pwvs = pwvs[ind]
sorted_lengths = lengths[ind]

target_time = np.sum(lengths) * fyear
total_time = 0
n = 0
while total_time < target_time:
    total_time += sorted_lengths[n]
    n += 1
pwv_limit = sorted_pwvs[n]
print(f"PWV limit for fyear = {fyear} is {pwv_limit} mm")

fig = plt.figure()
nrow, ncol = 2, 1

ax = fig.add_subplot(nrow, ncol, 1)
ax.plot(times / 86400, pwvs, '.')

ax = fig.add_subplot(nrow, ncol, 2)
good = pwvs < pwv_limit
ax.plot(times[good] / 86400, pwvs[good], '.')

good_times = times[good]
good_lengths = lengths[good]
t_month = 30.5 * 86400
print("Monthly observing efficiencies after fyear cut.")
for month in range(12):
    tstart = month * t_month
    tstop = tstart + t_month
    ind = np.logical_and(good_times > tstart, good_times < tstop)
    frac = np.sum(good_lengths[ind]) / t_month
    print(f"{month + 1:02} : {frac:.3f}")

# Add downtime
nsecond = 365 * 86400
good_seconds = np.ones(nsecond, dtype=bool)
breaks = []
while np.sum(good_seconds) / nsecond > fdown:
    # Add another break
    start = np.random.rand() * nsecond
    length = break_length + np.random.randn() * break_length_sigma
    stop = start + length
    breaks.append((start, stop))
    good_seconds[int(start) : int(stop)] = False

# write out the observations that survive the cuts

with open(fname_out, "w") as schedule_out:
    for line in header:
        schedule_out.write(line)
    for pwv, time, line in zip(pwvs, times, lines):
        if pwv > pwv_limit:
            continue
        itime = int(time)
        if itime >= nsecond or not good_seconds[itime]:
            continue
        schedule_out.write(line)
