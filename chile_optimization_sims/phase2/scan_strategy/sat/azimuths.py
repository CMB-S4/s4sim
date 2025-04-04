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


angles = [90]
for angle in angles:
    az = np.arange(360)
    azimuth_north = np.zeros(360)
    azimuth_south = np.zeros(360)
    schedule = toast.schedule.GroundSchedule()
    fname_schedule = f"schedule_sat.sun{angle}max.txt"
    schedule.read(fname_schedule)
    t0 = schedule.scans[0].start.timestamp()
    for scan in schedule.scans:
        tstart = scan.start.timestamp()
        # tstop = scan.stop.timestamp()
        # tdelta = tstop - tstart
        az_min = scan.az_min.to_value(u.deg)
        az_max = scan.az_max.to_value(u.deg)
        if az_max < 360:
            good = np.logical_and(az >= az_min, az <= az_max)
        else:
            good = np.logical_or(az >= az_min, az <= az_max - 360)
        #in_season = (tstart - t0) / 86400 > 90
        #azimuth_full[good] += 1
        #if in_season:
        if "North" in scan.name:
            azimuth_north[good] += 1
        else:
            azimuth_south[good] += 1

nrow, ncol = 1, 1
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
ax = fig.add_subplot(nrow, ncol, 1)
ax.plot(az + .5, azimuth_south, label="South")
ax.plot(az + .5, azimuth_north, label="North")
ax.set_xlabel("Azimuth [deg]")
ax.set_ylabel("N_obs")
ax.legend(loc="best")
ax.set_title(f"Azimuth in {fname_schedule}")
ax.grid()
fig.savefig("azimuth.png")
