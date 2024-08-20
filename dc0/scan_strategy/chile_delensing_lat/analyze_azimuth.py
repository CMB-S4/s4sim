import os
import sys

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import toast


weather = None
site_uid = toast.utils.name_UID("ATACAMA")
realization = 0

times = {}
schedule = toast.schedule.GroundSchedule()
# schedule.read(f"schedules/chile_schedule_sat.w_break.txt")
schedule.read(sys.argv[1])

tier1 = []
tier2 = []
for scan in schedule.scans:
    if scan.name.startswith("Tier1DEC-"):
        tier1.append(scan)
    elif scan.name.startswith("Tier2DEC-"):
        tier2.append(scan)

nrow, ncol = 1, 2
fig = plt.figure(figsize=[4 * ncol, 4 * nrow])

ax = fig.add_subplot(nrow, ncol, 1)
azmin = 1e10
azmax = -1e10
for scan in tier1:
    az1 = scan.az_min.to_value(u.deg)
    az2 = scan.az_max.to_value(u.deg)
    azmin = min(az1, azmin)
    azmax = max(az1, azmax)
    el = scan.el.to_value(u.deg)
    ax.plot([az1, az2], [el, el], 'k')
ax.set_xlim([0, 360])
ax.set_ylim([45, 65])
ax.set_xlabel("Azimuth [deg]")
ax.set_ylabel("Elevation [deg]")
ax.set_title(f"Tier1 tiles; az = [{azmin:.1f}, {azmax:.1f}]")

ax = fig.add_subplot(nrow, ncol, 2)
azmin = 1e10
azmax = -1e10
for scan in tier2:
    az1 = scan.az_min.to_value(u.deg)
    az2 = scan.az_max.to_value(u.deg)
    azmin = min(az1, azmin)
    azmax = max(az1, azmax)
    el = scan.el.to_value(u.deg)
    ax.plot([az1, az2], [el, el], 'k')
ax.set_xlim([0, 360])
ax.set_ylim([45, 65])
ax.set_xlabel("Azimuth [deg]")
ax.set_ylabel("Elevation [deg]")
ax.set_title(f"Tier2 tiles; az = [{azmin:.1f}, {azmax:.1f}]")
