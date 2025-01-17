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

dust = hp.read_map("/home/reijo/work/npipe6/npipe6v20_353_map.fits", None)
sdust = hp.smoothing(dust, fwhm=np.radians(1), lmax=512)
pdust = np.sqrt(sdust[1]**2 + sdust[2]**2)


# Load the schedules

schedule1 = toast.schedule.GroundSchedule()
#schedule1.read("../sat_max/schedule_sat.sun90.txt"); name1 = "sun90.v1"
#schedule1.read("../sat_max/schedule_sat.sun90.test.txt"); name1 = "sun90.test"
#schedule1.read("../sat_max/schedule_sat.sun90.test2.txt"); name1 = "sun90.test2"
#schedule1.read("../sat_max/schedule_sat.sun90.v4.txt"); name1 = "sun90.v4"
#schedule1.read("../sat_max/schedule_sat.sun90.v2.txt"); name1 = "sun90.v2"
#schedule1.read("../sat_max/schedule_sat.sun45.v2.txt"); name1 = "sun45.v2"
#schedule1.read("../../../phase2/scan_strategy/sat/schedule_sat.sun90.txt"); name1 = "phase2.sun90"
schedule1.read("../../../phase2/scan_strategy/sat/schedule_sat.sun90max.txt"); name1 = "phase2.sun90max"
#schedule1.read("../../../phase2/scan_strategy/sat/schedule_sat.sun90max3.txt"); name1 = "phase2.sun90max3"

schedule2 = toast.schedule.GroundSchedule()
#schedule2.read("schedule.txt"); name2 = "Alt"
#schedule2.read("../../../phase2/scan_strategy/sat/schedule_sat.sun90.txt"); name2 = "phase2.sun90"
#schedule2.read("../../../phase2/scan_strategy/sat/schedule_sat.sun45.txt"); name2 = "phase2.sun45"
#schedule2.read("../sat_max/schedule_sat.sun45.v2.txt"); name2 = "sun45.v2"
#schedule2.read("../sat_max/schedule_sat.sun90.v2.txt"); name2 = "sun90.v2"
#schedule2.read("../sat_max/schedule_sat.sun90.v3.txt"); name2 = "sun90.v3"
schedule2.read("../sat_max/schedule_sat.sun90.test2.txt"); name2 = "sun90.test2"
#schedule2.read("../../../phase2/scan_strategy/sat/schedule_sat.sun90max.txt"); name2 = "phase2.sun90max"
#schedule2.read("../../../phase2/scan_strategy/sat/schedule_sat.sun90max.txt"); name2 = "phase2.sun90max"

plotdir = f"plots_{name1}_vs_{name2}"
os.makedirs(plotdir, exist_ok=True)

# Create an observer at the site

observer = ephem.Observer()
observer.lon = schedule1.site_lon.to_value(u.radian)
observer.lat = schedule1.site_lat.to_value(u.radian)
observer.elevation = schedule1.site_alt.to_value(u.meter)
observer.epoch = ephem.J2000
observer.compute_pressure()

# Visualize what each schedule is observing

nside = 128
npix = 12 * nside**2
fov = 35 * u.deg
nday = 365
tstep = 600
naz = 10

fname1 = f"hitmap_{name1}.npy"
fname2 = f"hitmap_{name2}.npy"

if os.path.isfile(fname1):
    m1 = np.load(fname1)
else:
    m1 = np.zeros([nday, npix])
    for doy in range(nday):
        print(doy)
        t0 = schedule1.scans[0].start.timestamp() + doy * 86400
        for t in np.arange(t0, t0 + 86400, tstep):
            for scan in schedule1.scans:
                tstart = scan.start.timestamp()
                tstop = scan.stop.timestamp()
                if tstop < t or tstart > t:
                    continue
                azmin = scan.az_min.to_value(u.rad)
                azmax = scan.az_max.to_value(u.rad)
                el = scan.el.to_value(u.rad)
                observer.date = to_DJD(t)
                for az in np.linspace(azmin, azmax, naz):
                    ra, dec = observer.radec_of(az, el)
                    vec = hp.dir2vec(np.degrees(ra), np.degrees(dec), lonlat=True)
                    pix = hp.query_disc(
                        nside, vec, radius=fov.to_value(u.rad) / 2
                    )
                    m1[doy, pix] += 1
                break
    np.save(fname1, m1)

if os.path.isfile(fname2):
    m2 = np.load(fname2)
else:
    m2 = np.zeros([nday, npix])
    for doy in range(nday):
        print(doy)
        t0 = schedule1.scans[0].start.timestamp() + doy * 86400
        for t in np.arange(t0, t0 + 86400, tstep):
            for scan in schedule2.scans:
                tstart = scan.start.timestamp()
                tstop = scan.stop.timestamp()
                if tstop < t or tstart > t:
                    continue
                azmin = scan.az_min.to_value(u.rad)
                azmax = scan.az_max.to_value(u.rad)
                el = scan.el.to_value(u.rad)
                observer.date = to_DJD(t)
                for az in np.linspace(azmin, azmax, naz):
                    ra, dec = observer.radec_of(az, el)
                    vec = hp.dir2vec(np.degrees(ra), np.degrees(dec), lonlat=True)
                    pix = hp.query_disc(
                        nside, vec, radius=fov.to_value(u.rad) / 2
                    )
                    m2[doy, pix] += 1
                break
    np.save(fname2, m2)

m1sum = np.sum(m1, 0)
m2sum = np.sum(m2, 0)

ind = int(m1sum.size * 0.97)
limit1 = np.sort(m1sum)[ind]
mask1 = m1sum > limit1
limit2 = np.sort(m2sum)[ind]
mask2 = m2sum > limit2

a1 = np.argmax(m1sum)
a2 = np.argmax(m2sum)

nrow, ncol = 2, 4

fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
ax = fig.add_subplot(nrow, ncol, 1)
ax.set_title("Daily hits to deepest pixel")
ax.plot(m1[:, a1], label=name1)
ax.plot(m2[:, a2], label=name2)
ax.legend(loc="best")

ax = fig.add_subplot(nrow, ncol, 1 + ncol)
ax.set_title("Cumulative hits to deepest pixel")
ax.plot(np.cumsum(m1[:, a1]), label=f"{name1} : {m1sum[a1]}")
ax.plot(np.cumsum(m2[:, a2]), label=f"{name2} : {m2sum[a2]}")
ax.legend(loc="best")

ax = fig.add_subplot(nrow, ncol, 2)
ax.set_title("Daily hits to best 3%")
ax.plot(np.sum(m1[:, mask1], 1), label=name1)
ax.plot(np.sum(m2[:, mask2], 1), label=name2)
ax.legend(loc="best")

ax = fig.add_subplot(nrow, ncol, 2 + ncol)
ax.set_title("Cumulative hits to best 3%")
ax.plot(np.cumsum(np.sum(m1[:, mask1], 1)), label=f"{name1} : {np.sum(m1sum[mask1])}")
ax.plot(np.cumsum(np.sum(m2[:, mask2], 1)), label=f"{name2} : {np.sum(m2sum[mask2])}")
ax.legend(loc="best")


m1sum[m1sum == 0] = hp.UNSEEN
m2sum[m2sum == 0] = hp.UNSEEN
hp.mollview(m1sum, sub=[nrow, ncol, 3], cmap="magma", title=name1)
hp.mollview(m2sum, sub=[nrow, ncol, 3 + ncol], cmap="magma", title=name2)

smask1 = hp.smoothing(mask1, fwhm=np.radians(2), lmax=256) 
smask2 = hp.smoothing(mask2, fwhm=np.radians(2), lmax=256)
for mask in smask1, smask2:
    mask[mask > 0.95] = 0
    mask[mask > 0.05] = 1
    mask[mask < 0.05] = 0

hp.mollview(pdust, sub=[nrow, ncol, 4], max=1e-4, cmap="magma", coord="gc", title="")
hp.mollview(smask1, sub=[nrow, ncol, 4], min=0, max=1, cmap="bwr", cbar=False, title=name1, reuse_axes=True, alpha=0.75*smask1)
hp.mollview(pdust, sub=[nrow, ncol, 4 + ncol], max=1e-4, cmap="magma", coord="gc", title="")
hp.mollview(smask2, sub=[nrow, ncol, 4 + ncol], min=0, max=1, cmap="bwr", cbar=False, title=name2, reuse_axes=True, alpha=0.75*smask2)

plt.tight_layout()

plt.savefig(f"best_pixel.{name1}_vs_{name2}.png")

sys.exit()

sun = ephem.Sun()
moon = ephem.Moon()

def plot_sso(sso, radius, color, n=100, avoid=True):
    ra = np.degrees(sso.ra)
    dec = np.degrees(sso.dec)
    hp.projplot(ra, dec, "o", color=color, ms=radius / 3, lonlat=True)
    if avoid:
        n = 100
        thetas = np.ones(n) * np.radians(radius)
        phis = np.linspace(0, 2 * np.pi, n)
        hp.projplot(thetas, phis, "-", color=color, rot=[ra, dec, 0], zorder=100)
    return

nrow, ncol = 1, 2
for doy in range(nday):
    print(doy)
    t = schedule1.scans[0].start.timestamp() + doy * 86400
    observer.date = to_DJD(t)
    sun.compute(observer)
    moon.compute(observer)
    fig = plt.figure(figsize=[ncol * 4, nrow * 4])
    fig.suptitle(f"DOY {doy}")
    mm = m1[doy]
    mm[mm == 0] = hp.UNSEEN
    h1 = m1[doy, a1]
    hp.mollview(mm, sub=[nrow, ncol, 1], title=f"{name1}, {h1}")
    #if np.degrees(sun.alt) > 5:
    plot_sso(sun, 90, "yellow")
    plot_sso(moon, 45, "white")
    hp.projplot(*hp.pix2ang(nside, a1), "x", color="black", ms=10)
    #hp.projplot(np.degrees(moon.ra), np.degrees(moon.dec), "o", color="white", ms=5, lonlat=True)
    mm = m2[doy]
    mm[mm == 0] = hp.UNSEEN
    h2 = m2[doy, a2]
    hp.mollview(mm, sub=[nrow, ncol, 2], title=f"{name2}, {h2}, {h2 / h1:.3f}")
    plot_sso(sun, 90, "yellow")
    plot_sso(moon, 45, "white")
    hp.projplot(*hp.pix2ang(nside, a2), "x", color="black", ms=10)
    fname = f"{plotdir}/day_{doy:03}.png"
    fig.savefig(fname)
    plt.close()


# Minute-by-minute plot

tstep = 600
nstep = 86400 // tstep
naz = 100

def plot_sky(observer):
    ra, dec = observer.radec_of(0, np.pi / 2)
    vec = hp.dir2vec(np.degrees(ra), np.degrees(dec), lonlat=True)
    margin = fov.to_value(u.deg) / 2
    pix1 = hp.query_disc(nside, vec, radius=np.radians(40 + margin))
    pix2 = hp.query_disc(nside, vec, radius=np.radians(20 - margin))
    mask = np.zeros(npix)
    mask[pix1] = 1
    mask[pix2] = 0
    hp.mollview(mask, cmap="magma", cbar=False, title="", reuse_axes=True, alpha=mask * .25)
    return mask


for doy in range(0, 366, 15):
    m1day = np.zeros([nstep, npix])
    m2day = np.zeros([nstep, npix])

    t0 = schedule1.scans[0].start.timestamp() + doy * 86400
    for istep, t in enumerate(np.arange(t0, t0 + 86400, tstep)):
        for scan in schedule1.scans:
            tstart = scan.start.timestamp()
            tstop = scan.stop.timestamp()
            if tstop < t or tstart > t:
                continue
            azmin = scan.az_min.to_value(u.rad)
            azmax = scan.az_max.to_value(u.rad)
            el = scan.el.to_value(u.rad)
            observer.date = to_DJD(t)
            for az in np.linspace(azmin, azmax, naz):
                ra, dec = observer.radec_of(az, el)
                vec = hp.dir2vec(np.degrees(ra), np.degrees(dec), lonlat=True)
                pix = hp.query_disc(
                    nside, vec, radius=fov.to_value(u.rad) / 2
                )
                m1day[istep, pix] += 1
            break

    for istep, t in enumerate(np.arange(t0, t0 + 86400, tstep)):
        for scan in schedule2.scans:
            tstart = scan.start.timestamp()
            tstop = scan.stop.timestamp()
            if tstop < t or tstart > t:
                continue
            azmin = scan.az_min.to_value(u.rad)
            azmax = scan.az_max.to_value(u.rad)
            el = scan.el.to_value(u.rad)
            observer.date = to_DJD(t)
            for az in np.linspace(azmin, azmax, naz):
                ra, dec = observer.radec_of(az, el)
                vec = hp.dir2vec(np.degrees(ra), np.degrees(dec), lonlat=True)
                pix = hp.query_disc(
                    nside, vec, radius=fov.to_value(u.rad) / 2
                )
                m2day[istep, pix] += 1
            break

    m1daysum = np.sum(m1day, 0)
    m2daysum = np.sum(m2day, 0)

    nrow, ncol = 1, 2
    fig = plt.figure(figsize=[ncol * 4, nrow * 4])
    for istep, t in enumerate(np.arange(t0, t0 + 86400, tstep)):
        print(f"{istep} / {nstep}")
        observer.date = to_DJD(t)
        sun.compute(observer)
        moon.compute(observer)
        tt = istep * tstep
        hour = int(tt / 3600)
        minute = (int((tt - hour * 3600) / 60))
        fig.suptitle(f"DOY {doy}, {hour:02}:{minute:02}")
        mm = m1day[istep]
        mm[mm == 0] = hp.UNSEEN
        hp.mollview(mm, sub=[nrow, ncol, 1], title=name1)
        plot_sky(observer)
        plot_sso(sun, 90, "yellow", avoid=np.degrees(sun.alt) > 5)
        plot_sso(moon, 45, "white", avoid=np.degrees(moon.alt) > 5)
        hp.projplot(*hp.pix2ang(nside, a1), "x", color="black", ms=10)
        mm = m2day[istep]
        mm[mm == 0] = hp.UNSEEN
        hp.mollview(mm, sub=[nrow, ncol, 2], title=name2)
        plot_sky(observer)
        plot_sso(sun, 90, "yellow", avoid=np.degrees(sun.alt) > 5)
        plot_sso(moon, 45, "white", avoid=np.degrees(moon.alt) > 5)
        hp.projplot(*hp.pix2ang(nside, a2), "x", color="black", ms=10)
        fname = f"{plotdir}/day_{doy:03}_{istep:03}.png"
        fig.savefig(fname)
        plt.clf()
