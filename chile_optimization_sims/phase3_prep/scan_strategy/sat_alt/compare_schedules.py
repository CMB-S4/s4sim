import os
import sys

import ephem
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

import toast
from toast import qarray as qa
from toast.coordinates import to_DJD


tol = 1 * u.deg   # Add 1 deg tolerance
sun_avoidance_altitude = 5 * u.deg + tol
avoidance = {
    "Sun" : 90 * u.deg - tol,
    "Moon" : 45 * u.deg - tol,
}

# Load the schedules

schedule1 = toast.schedule.GroundSchedule()
#schedule1.read("../sat_max/schedule_sat.sun90.test.txt"); name1 = "sun90.test"
#schedule1.read("../sat_max/schedule_sat.sun45.v2.txt"); name1 = "sun45.v2"
schedule1.read("../../../phase2/scan_strategy/sat/schedule_sat.sun90max3.txt"); name1 = "phase2.sun90max3"

schedule2 = toast.schedule.GroundSchedule()
#schedule2.read("schedule.txt"); name2 = "Alt"
schedule2.read("../../../phase2/scan_strategy/sat/schedule_sat.sun90.txt"); name2 = "phase2.sun90"
#schedule2.read("../../../phase2/scan_strategy/sat/schedule_sat.sun45.txt"); name2 = "phase2.sun45"
#schedule2.read("../sat_max/schedule_sat.sun90.v3.txt"); name2 = "sun90.v3"


plotdir = f"plots_{name1}_vs_{name2}"
os.makedirs(plotdir, exist_ok=True)

# Create an observer at the site

observer = ephem.Observer()
observer.lon = schedule1.site_lon.to_value(u.radian)
observer.lat = schedule1.site_lat.to_value(u.radian)
observer.elevation = schedule1.site_alt.to_value(u.meter)
observer.epoch = ephem.J2000
observer.compute_pressure()

# Plot Az, El and distance to the Sun and the Moon

nrow, ncol = 2, 5
fig = plt.figure(figsize=[ncol * 4, nrow * 4])

# plot el

ax = fig.add_subplot(nrow, ncol, 1)
ax.set_title(f"{name1} El")
for scan in schedule1.scans:
    tstart = scan.start.timestamp()
    tstop = scan.stop.timestamp()
    azmin = scan.az_min.to_value(u.deg)
    azmax = scan.az_max.to_value(u.deg)
    el = scan.el.to_value(u.deg)
    ax.plot([tstart, tstop], [el, el])

ax = fig.add_subplot(nrow, ncol, 1 + ncol)
ax.set_title(f"{name2} El")
for scan in schedule2.scans:
    tstart = scan.start.timestamp()
    tstop = scan.stop.timestamp()
    azmin = scan.az_min.to_value(u.deg)
    azmax = scan.az_max.to_value(u.deg)
    el = scan.el.to_value(u.deg)
    ax.plot([tstart, tstop], [el, el])

# plot az

ax = fig.add_subplot(nrow, ncol, 2)
az_mount = []
az_sky = []
for scan in schedule1.scans:
    tstart = scan.start.timestamp()
    tstop = scan.stop.timestamp()
    azmin = scan.az_min.to_value(u.deg)
    azmax = scan.az_max.to_value(u.deg)
    el = scan.el.to_value(u.rad)
    az_mount.append(azmax - azmin)
    az_sky.append((azmax - azmin) * np.sin(el))
    ax.fill_between([tstart, tstop], [azmin, azmin], [azmax, azmax])
az_mount = np.mean(az_mount)
az_sky = np.mean(az_sky)
ax.set_title(f"{name1} Az, throw = {az_mount:.3f} {az_sky:.3f}")

ax = fig.add_subplot(nrow, ncol, 2 + ncol)
az_mount = []
az_sky = []
for scan in schedule2.scans:
    tstart = scan.start.timestamp()
    tstop = scan.stop.timestamp()
    azmin = scan.az_min.to_value(u.deg)
    azmax = scan.az_max.to_value(u.deg)
    el = scan.el.to_value(u.rad)
    az_mount.append(azmax - azmin)
    az_sky.append((azmax - azmin) * np.sin(el))
    ax.fill_between([tstart, tstop], [azmin, azmin], [azmax, azmax])
az_mount = np.mean(az_mount)
az_sky = np.mean(az_sky)
ax.set_title(f"{name2} Az, throw = {az_mount:.3f} {az_sky:.3f}")

# plot solar distance

def from_angles(az, el):
    elquat = qa.rotation(YAXIS, np.radians(90 - el))
    azquat = qa.rotation(ZAXIS, np.radians(az))
    return qa.mult(azquat, elquat)

XAXIS, YAXIS, ZAXIS = np.eye(3)

sun = ephem.Sun()
moon = ephem.Moon()

print("Measuring SSO distances")

ax1 = fig.add_subplot(nrow, ncol, 3)
ax2 = fig.add_subplot(nrow, ncol, 4)

bad = {
    "Sun" : {"count" : 0, "time" : 0},
    "Moon" : {"count" : 0, "time" : 0},
}

for scan in schedule1.scans:
    tstart = scan.start.timestamp()
    tstop = scan.stop.timestamp()
    azmin = scan.az_min.to_value(u.deg)
    azmax = scan.az_max.to_value(u.deg)
    el = scan.el.to_value(u.deg)

    azmin, azmax = np.unwrap([azmin, azmax], period=360)

    times = np.linspace(tstart, tstop, 10)

    quats = []
    for az in np.linspace(azmin, azmax, 10):
        quats.append(from_angles(az % 360, el))
    vecs = qa.rotate(quats, ZAXIS)

    min_dists = {"Sun" : [], "Moon" : []}
    in_violation = set()
    for t in times:
        observer.date = to_DJD(t)
        for name, sso in [("Sun", sun), ("Moon", moon)]:
            sso.compute(observer)
            sso_az, sso_el = np.degrees(sso.az), np.degrees(sso.alt)
            if sso_el > sun_avoidance_altitude.to_value(u.deg):
                sso_quat = from_angles(sso_az, sso_el)
                sso_vec = qa.rotate(sso_quat, ZAXIS)
                dpmax = np.amax(np.dot(vecs, sso_vec))
                min_dist = np.degrees(np.arccos(dpmax))
                #if name == "Sun" and min_dist < 80:
                #    print("Solar avoidance violation")
                #    sys.exit()
                if min_dist < avoidance[name].to_value(u.deg):
                    in_violation.add(name)
            else:
                min_dist = np.nan
            min_dists[name].append(min_dist)
    for name in in_violation:
        bad[name]["count"] += 1
        bad[name]["time"] += tstop - tstart

    ax1.plot(times, min_dists["Sun"])
    ax2.plot(times, min_dists["Moon"])

nbad, tbad = bad["Sun"]["count"], bad["Sun"]["time"]
ax1.set_title(f"{name1} Solar dist {nbad} violations, {tbad / 3600:.1f} h")
ax1.axhline(avoidance["Sun"].to_value(u.deg), linestyle="--", color="k")
nbad, tbad = bad["Moon"]["count"], bad["Moon"]["time"]
ax2.set_title(f"{name1} Lunar dist {nbad} violations, {tbad / 3600:.1f} h")
ax2.axhline(avoidance["Moon"].to_value(u.deg), linestyle="--", color="k")

ax1 = fig.add_subplot(nrow, ncol, 3 + ncol)
ax2 = fig.add_subplot(nrow, ncol, 4 + ncol)


bad = {
    "Sun" : {"count" : 0, "time" : 0},
    "Moon" : {"count" : 0, "time" : 0},
}

for scan in schedule2.scans:
    tstart = scan.start.timestamp()
    tstop = scan.stop.timestamp()
    azmin = scan.az_min.to_value(u.deg)
    azmax = scan.az_max.to_value(u.deg)
    el = scan.el.to_value(u.deg)

    azmin, azmax = np.unwrap([azmin, azmax], period=360)

    times = np.linspace(tstart, tstop, 10)

    quats = []
    for az in np.linspace(azmin, azmax, 10):
        quats.append(from_angles(az % 360, el))
    vecs = qa.rotate(quats, ZAXIS)

    min_dists = {"Sun" : [], "Moon" : []}
    in_violation = set()
    for t in times:
        observer.date = to_DJD(t)
        for name, sso in [("Sun", sun), ("Moon", moon)]:
            sso.compute(observer)
            sso_az, sso_el = np.degrees(sso.az), np.degrees(sso.alt)
            if sso_el > sun_avoidance_altitude.to_value(u.deg):
                sso_quat = from_angles(sso_az, sso_el)
                sso_vec = qa.rotate(sso_quat, ZAXIS)
                dpmax = np.amax(np.dot(vecs, sso_vec))
                min_dist = np.degrees(np.arccos(dpmax))
                if min_dist < avoidance[name].to_value(u.deg):
                    in_violation.add(name)
                #if name == "Sun" and min_dist < 80:
                #    print("Solar avoidance violation")
                #    sys.exit()
            else:
                min_dist = np.nan
            min_dists[name].append(min_dist)
    for name in in_violation:
        bad[name]["count"] += 1
        bad[name]["time"] += tstop - tstart

    ax1.plot(times, min_dists["Sun"])
    ax2.plot(times, min_dists["Moon"])

nbad, tbad = bad["Sun"]["count"], bad["Sun"]["time"]
ax1.set_title(f"{name2} Solar dist {nbad} violations, {tbad / 3600:.1f} h")
ax1.axhline(avoidance["Sun"].to_value(u.deg), linestyle="--", color="k")
nbad, tbad = bad["Moon"]["count"], bad["Moon"]["time"]
ax2.set_title(f"{name2} Lunar dist {nbad} violations, {tbad / 3600:.1f} h")
ax2.axhline(avoidance["Moon"].to_value(u.deg), linestyle="--", color="k")

# plot daily and cumulative integration time

ax1 = fig.add_subplot(nrow, ncol, 5)
ax2 = fig.add_subplot(nrow, ncol, 5 + ncol)
ax1.set_title("Cumulative integration time")
ax2.set_title("Daily integration time")
t0 = schedule1.scans[0].start.timestamp()
x = [t0]
y = [0]
daily = np.zeros(366)
for scan in schedule1.scans:
    tstart = scan.start.timestamp()
    tstop = scan.stop.timestamp()
    daily[int((tstart - t0) / 86400)] += tstop - tstart
    x.append(tstop)
    y.append(tstop - tstart)
ax1.plot(np.array(x) - t0, np.cumsum(y) / 86400, label=name1)
ax2.plot(daily / 3600, label=name1)

t0 = schedule2.scans[0].start.timestamp()
x = [t0]
y = [0]
daily = np.zeros(366)
for scan in schedule2.scans:
    tstart = scan.start.timestamp()
    tstop = scan.stop.timestamp()
    daily[int((tstart - t0) / 86400)] += tstop - tstart
    x.append(tstop)
    y.append(tstop - tstart)
ax1.plot(np.array(x) - t0, np.cumsum(y) / 86400, label=name2)
ax2.plot(daily / 3600, label=name2)
ax1.legend(loc="best")
ax2.legend(loc="best")

fig.tight_layout()
fig.savefig(f"{plotdir}/schedule_comparison.{name1}_vs_{name2}.png")
