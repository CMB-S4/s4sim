import glob
import os
import sys

import astropy.units as u
import numpy as np
import toast


include_off_season = False
south_only = True

times = {}
fname_schedule = sys.argv[1]

if not include_off_season:
    print("WARNING: not including January - March")
if south_only:
    print("WARNING: not including Northern observations")

n_tot = 0
n_good = 0
t_good = 0
t_tot = 0
t_scan_tot = 0
t_turn_tot = 0
t_start = None
t_stop = None
el_min = 90
el_max = 0
el_sum = 0

schedule = toast.schedule.GroundSchedule()
schedule.read(fname_schedule)
for scan in schedule.scans:
    n_tot += 1
    if t_start is None:
        t_start = scan.start
    doy = scan.start.timetuple().tm_yday
    if doy < 31 + 28 + 31 and not include_off_season:
        continue
    az_min = scan.az_min.to_value(u.deg)
    az_max = scan.az_max.to_value(u.deg)
    az_mean = 0.5 * (az_min + az_max)
    if az_mean < 90 or az_mean > 270 and south_only:
        continue
    t_stop = scan.stop
    t_delta = (scan.stop - scan.start).total_seconds()
    t_tot += t_delta

    az_delta = scan.az_max - scan.az_min
    el = scan.el.to_value(u.deg)
    el_min = min(el_min, el)
    el_max = max(el_max, el)
    el_sum += el * t_delta


t_delta_tot = (t_stop - t_start).total_seconds()
if not include_off_season:
    t_delta_tot -= (31 + 28 + 31) * 86400
print(
    f"Scheduled observing time: {t_tot / 86400:.3f} / {t_delta_tot / 86400:.3f} "
    f"days = {t_tot / t_delta_tot:.3f}"
)

el_mean = el_sum / t_delta_tot
print(f"el_min = {el_min:.1f} deg, el_max = {el_max:.1f} deg, el_mean = {el_mean:.1f} deg")


print(f"Number of observations: {n_tot}")
