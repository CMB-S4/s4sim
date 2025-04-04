import os
import sys

import ephem
import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

import toast
from toast import qarray as qa
from toast.coordinates import to_DJD


schedule = toast.schedule.GroundSchedule()
schedule.read("schedule_lat_delensing_sun90max.txt"); name = "delensing_sun90max"

days = []
tobs = []
for doy in range(365):
    print(doy)
    t0 = schedule.scans[0].start.timestamp() + doy * 86400
    tday = 0
    for scan in schedule.scans:
        tstart = scan.start.timestamp()
        tstop = scan.stop.timestamp()
        if tstart < t0 or tstart > t0 + 86400:
            continue
        if (
                scan.az_min.to_value(u.deg) > 90
                and scan.az_max.to_value(u.deg) < 270
        ):
            tday += tstop - tstart
    days.append(doy)
    tobs.append(tday / 3600)

plt.plot(days, tobs)
