import glob
import os
import sys

import astropy.units as u
import numpy as np
import toast


weather = None
site_uid = toast.utils.name_UID("ATACAMA")
realization = 0

times = {}
fname_schedules = glob.glob("*/schedule*.txt")

for fname_schedule in fname_schedules:
    print("")
    n_tot = 0
    n_good = 0
    t_good = 0
    t_tot = 0
    t_scan_tot = 0
    t_turn_tot = 0
    t_start = None
    t_stop = None

    if "sat" in fname_schedule:
        scan_rate = 1.0 * u.deg / u.s
        az_accel = 0.97 * u.deg / u.s**2
        modulated = False
    elif "wide" in fname_schedule:
        base_rate = 0.55 * u.deg / u.s
        az_accel = 2.58 * u.deg / u.s**2
        modulated = True
    else:
        scan_rate = 1.0 * u.deg / u.s
        az_accel = 2 * u.deg / u.s**2
        modulated = False

    schedule = toast.schedule.GroundSchedule()
    schedule.read(fname_schedule)
    for scan in schedule.scans:
        n_tot += 1
        if t_start is None:
            t_start = scan.start
        t_stop = scan.stop
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

        az_min = scan.az_min.to_value(u.rad)
        az_max = scan.az_max.to_value(u.rad)
        az_delta = scan.az_max - scan.az_min
        el = scan.el.to_value(u.rad)
        t_delta = (scan.stop - scan.start).total_seconds()
        if pwv < 3.0:
            n_good += 1
            t_good += t_delta
        if modulated:
            t_scan = (np.cos(az_min) - np.cos(az_max)) / base_rate.to_value(u.rad / u.s) * u.s
            # scan rate at the start of a turnaround
            az_rate = base_rate / np.abs(np.sin(az_min))
        else:
            t_scan = az_delta / az_rate
            az_rate = scan_rate / np.cos(el)
        t_turn = 2 * az_rate / az_accel

        t_tot += t_delta
        t_scan_tot += t_delta * t_scan
        t_turn_tot += t_delta * t_turn

    t_delta_tot = (t_stop - t_start).total_seconds()
    t_scan_tot /= t_tot
    t_turn_tot /= t_tot
    t_throw_tot = t_scan_tot + t_turn_tot
    print(
        f"Scheduled observing time: {t_tot / 86400:.3f} / {t_delta_tot / 86400:.3f} "
        f"days = {t_tot / t_delta_tot:.3f}"
    )
    print(
        f"Good observing time: {t_good / 86400:.3f} / {t_tot / 86400:.3f} "
        f"days = {t_good / t_tot:.3f}"
    )
    print(f"Number of good observations: {n_good} / {n_tot} = {n_good / n_tot:.3f}")
    if modulated:
        print(f"Base rate (on-sky): {base_rate}")
        print(f"Max rate (on-sky): {az_rate:.3f}")
    else:
        print(f"Scan rate (on-sky): {scan_rate}")
    print(f"Scan accel (mount): {az_accel}")
    print(
        f"Scan time / throw time: {t_scan_tot:.3f} / {t_throw_tot:.3f} "
        f"= {t_scan_tot / t_throw_tot:.3f}"
    )
