import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

import toast


flavors = "baseline", "no_break", "no_lunar_avoidance", "solar90", "el_min_40", "sidereal", "all"

weather = None
site_uid = toast.utils.name_UID("ATACAMA")
realization = 0

times = {}
for flavor in flavors:
    schedule = toast.schedule.GroundSchedule()
    schedule.read(f"schedules/{flavor}.txt")
    x = []
    y = []
    z = []
    tot = 0
    for scan in schedule.scans:
        # start_time = dateutil.parser.parse(f"{start_date} {start_time} +0000")
        # stop_time = dateutil.parser.parse(f"{stop_date} {stop_time} +0000")
        # mid_time = start_time + 0.5 * (stop_time - start_time)
        mid_time = scan.start + 0.5 * (scan.stop - scan.start)
        if weather is None:
            weather = toast.weather.SimWeather(
                time=mid_time,
                name="atacama",
                site_uid=site_uid,
                realization=realization,
            )
        else:
            weather.set(time=mid_time, realization=realization, site_uid=site_uid)
        pwv = weather.pwv.to_value("mm")
        if pwv > 3.0:
            continue

        x.append(scan.stop)
        td = (scan.stop - scan.start).seconds
        y.append(td)
        tot += td
        z.append(tot)
    times[flavor] = [np.array(x), np.array(y), np.array(z)]

nrow, ncol = 1, 1
fig = plt.figure(figsize=[8 * ncol, 4 * nrow])
ax = fig.add_subplot(nrow, ncol, 1)
for flavor in flavors:
    ax.plot(times[flavor][0], times[flavor][2] / 86400, label=flavor)
ax.set_title("Cumulative observing time")
ax.set_xlabel("Date")
ax.set_ylabel("[days]")
ax.legend(loc="best")

fig.tight_layout()
fname = "observing_time.png"
plt.savefig(fname)

print(f"Plot written to {fname}")
