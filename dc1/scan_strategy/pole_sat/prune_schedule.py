import dateutil
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

import toast.weather
import toast.utils


np.random.seed(673653982)
realization = 0

# Additional, random downtime
fdown = 0.864

# Fraction of the year surviving the weather cuts
if False:
    # Measure the pwv_limit based on fyear
    fyear = 0.6
    pwv_limit = None
else:
    # Use a fixed pwv_limit
    fyear = None
    pwv_limit = 6.0

fname_in = "schedules/pole_schedule_sat.txt"
fname_out = "schedules/pole_schedule_sat.pruned.txt"

first_time = dateutil.parser.parse(f"2027-01-01 00:00:00 00:00:00 +0000")

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
        mid_time = start_time + 0.5 * (stop_time - start_time)
        if weather is None:
            weather = toast.weather.SimWeather(
                time=mid_time,
                name="south_pole",
                site_uid=site_uid,
                realization=realization,
            )
        else:
            weather.set(time=mid_time, site_uid=site_uid, realization=realization)

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

if fyear is not None:
    target_time = np.sum(lengths) * fyear
    total_time = 0
    n = 0
    while total_time < target_time:
        total_time += sorted_lengths[n]
        n += 1
    pwv_limit = sorted_pwvs[n]
    print(f"PWV limit for fyear = {fyear} is {pwv_limit} mm")
else:
    n = np.searchsorted(sorted_pwvs, pwv_limit)
    fyear = np.sum(lengths[:n]) / np.sum(lengths)
    print(f"Fyear for PWV limit = {pwv_limit} mm is {fyear}")

# Plot

fig = plt.figure(figsize=[18, 12])
nrow, ncol = 4, 1

ax = fig.add_subplot(nrow, ncol, 1)
ax.plot(times / 86400, pwvs, '.')
ax.axhline(pwv_limit, color='k', linestyle="--")
ax.set_xlabel("DOY")
ax.set_ylabel("PWV [mm]")
ax.set_title(f"All {len(lengths)} observations")
ax.set_xlim([0, 366])

ax = fig.add_subplot(nrow, ncol, 2)
good = pwvs < pwv_limit
ax.plot(times[good] / 86400, pwvs[good], '.')
ax.set_xlabel("DOY")
ax.set_ylabel("PWV [mm]")
ax.set_title(f"{n} observations after PWV cut")
ax.set_xlim([0, 366])

# Compute observing efficiencies

good_times = times[good]
good_lengths = lengths[good]
t_month = 30.5 * 86400
x, y1, y2 = [], [], []  # for plotting
print("Monthly observing efficiencies after fyear cut.")
for month in range(12):
    tstart = month * t_month
    tstop = tstart + t_month
    ind = np.logical_and(good_times > tstart, good_times < tstop)
    frac = np.sum(good_lengths[ind]) / t_month
    print(f"{month + 1:02} : {frac:.3f}")
    # Uncut, for plotting
    ind2 = np.logical_and(times > tstart, times < tstop)
    frac2 = np.sum(lengths[ind2]) / t_month
    x += [tstart / 86400, tstop / 86400]
    y1 += [frac, frac]
    y2 += [frac2, frac2]

ax = fig.add_subplot(nrow, ncol, 3)
ax.plot(x, y1, label="w/ PWV limit")
ax.plot(x, y2, label="w/o PWV limit")
ax.set_ylim([0, 1.1])
ax.axhline(1, color="k")
ax.legend(loc="best")
ax.set_xlabel("DOY")
ax.set_ylabel("Obs. efficiency")
ax.set_xlim([0, 366])

"""
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
"""

# write out the observations that survive the cuts

final_times = []
final_pwvs = []

with open(fname_out, "w") as schedule_out:
    for line in header:
        schedule_out.write(line)
    for pwv, time, line in zip(pwvs, times, lines):
        if pwv > pwv_limit:
            continue
        #itime = int(time)
        #if itime >= nsecond or not good_seconds[itime]:
        #    continue
        if np.random.rand() < fdown:
            schedule_out.write(line)
            final_times.append(time)
            final_pwvs.append(pwv)

ax = fig.add_subplot(nrow, ncol, 4)
ax.plot(np.array(final_times) / 86400, final_pwvs, '.')
ax.set_xlabel("DOY")
ax.set_ylabel("PWV [mm]")
ax.set_title(f"{len(final_times)} observations after breaks")
ax.set_xlim([0, 366])

fig.subplots_adjust(hspace=0.4)
fig.savefig("pwv_limit_SPSAT.png")
