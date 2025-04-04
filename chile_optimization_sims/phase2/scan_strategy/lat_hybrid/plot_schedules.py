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


plot_el = False
cut_altiplanic_winter = True
plot_bk = False
nside = 128
npix = 12 * nside**2
fov = 8 * u.deg
nday = 365
tstep = 600
daz = 3
naz = 10

try:
    dust = hp.read_map("/global/cfs/cdirs/cmb/data/planck2020/npipe/npipe6v20/npipe6v20_353_map.fits", None)
except:
    dust = hp.read_map("/home/reijo/work/npipe6/npipe6v20_353_map.fits", None)
sdust = hp.smoothing(dust, fwhm=np.radians(1), lmax=512)
pdust = np.sqrt(sdust[1]**2 + sdust[2]**2)

if plot_bk:
    rhits = hp.read_map("../bk18_mask_largefield_cel_n0512.fits")
    rhits[np.isnan(rhits)] = 0
    sorted_hits = np.sort(rhits)
    limit = sorted_hits[int(0.97 * rhits.size)]
    bkmask = rhits > limit
    sbkmask = hp.smoothing(bkmask, fwhm=np.radians(2), lmax=256) 
    sbkmask[sbkmask > 0.95] = 0
    sbkmask[sbkmask > 0.05] = 1
    sbkmask[sbkmask < 0.05] = 0

# Load the schedules

schedule1 = toast.schedule.GroundSchedule()
#schedule1.read("schedule_lat_delensing_sun90max+wide.txt"); name1 = "delensing_sun90max+wide"
schedule1.read("schedule_lat_delensing_sun90max.txt"); name1 = "delensing_sun90max"

schedule2 = toast.schedule.GroundSchedule()
#schedule2.read("schedule_lat_delensing_sun90bk+wide.txt"); name2 = "delensing_sun90bk+wide"
schedule2.read("schedule_lat_delensing_sun90bk.txt"); name2 = "delensing_sun90bk"

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

def get_hits(name, schedule):
    fname_radec = f"hitmap_{name}.npy"
    fname_azel = f"hitmap_{name}.azel.npy"
    if os.path.isfile(fname_radec) and os.path.isfile(fname_azel):
        m = np.load(fname_radec)
        m_azel = np.load(fname_azel)
    else:
        m = np.zeros([nday, npix])
        m_azel = np.zeros([nday, npix])
        for doy in range(nday):
            print(doy)
            t0 = schedule.scans[0].start.timestamp() + doy * 86400
            for t in np.arange(t0, t0 + 86400, tstep):
                for scan in schedule.scans:
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
                        m[doy, pix] += 1
                        vec = hp.dir2vec(np.degrees(az), np.degrees(el), lonlat=True)
                        pix = hp.query_disc(
                            nside, vec, radius=fov.to_value(u.rad) / 2
                        )
                        m_azel[doy, pix] += 1
                    break
        np.save(fname_radec, m)
        np.save(fname_azel, m_azel)
    return m, m_azel

m1, m1azel = get_hits(name1, schedule1)
m2, m2azel = get_hits(name2, schedule2)

if cut_altiplanic_winter:
    for m in m1, m2, m1azel, m2azel:
        m[:31 + 28 + 31, :] = 0

m1sum = np.sum(m1, 0)
m2sum = np.sum(m2, 0)

m1azelsum = np.sum(m1azel, 0)
m2azelsum = np.sum(m2azel, 0)

ind = int(m1sum.size * 0.97)
limit1 = np.sort(m1sum)[ind]
mask1 = m1sum > limit1
limit2 = np.sort(m2sum)[ind]
mask2 = m2sum > limit2

a1 = np.argmax(m1sum)
a2 = np.argmax(m2sum)

if plot_el:
    nrow, ncol = 2, 5
else:
    nrow, ncol = 2, 4

fig = plt.figure(figsize=[4 * ncol, 4 * nrow])

ax = fig.add_subplot(nrow, ncol, 1)
ax.set_title("Daily hits to deepest pixel")
ax.plot(m1[:, a1], label=name1)
ax.plot(m2[:, a2], label=name2)
ax.set_xlabel("DOY")
ax.legend(loc="best")

ax = fig.add_subplot(nrow, ncol, 1 + ncol)
ax.set_title("Cumulative hits to deepest pixel")
ax.plot(np.cumsum(m1[:, a1]), label=f"{name1} : {m1sum[a1]}")
ax.plot(np.cumsum(m2[:, a2]), label=f"{name2} : {m2sum[a2]}")
ax.set_xlabel("DOY")
ax.legend(loc="best")

ax = fig.add_subplot(nrow, ncol, 2)
ax.set_title("Daily hits to best 3%")
ax.plot(np.sum(m1[:, mask1], 1), label=name1)
ax.plot(np.sum(m2[:, mask2], 1), label=name2)
ax.set_xlabel("DOY")
ax.legend(loc="best")

ax = fig.add_subplot(nrow, ncol, 2 + ncol)
ax.set_title("Cumulative hits to best 3%")
ax.plot(np.cumsum(np.sum(m1[:, mask1], 1)), label=f"{name1} : {np.sum(m1sum[mask1])}")
ax.plot(np.cumsum(np.sum(m2[:, mask2], 1)), label=f"{name2} : {np.sum(m2sum[mask2])}")
ax.set_xlabel("DOY")
ax.legend(loc="best")

for name, msum, col in [(name1, m1sum, 3), (name2, m2sum, 3 + ncol)]:
    msum[msum == 0] = hp.UNSEEN
    hp.mollview(
        msum,
        sub=[nrow, ncol, col],
        cmap="magma",
        title=name,
        unit="Hits",
    )

smask1 = hp.smoothing(mask1, fwhm=np.radians(2), lmax=256) 
smask2 = hp.smoothing(mask2, fwhm=np.radians(2), lmax=256)
for mask in smask1, smask2:
    mask[mask > 0.95] = 0
    mask[mask > 0.05] = 1
    mask[mask < 0.05] = 0

for name, smask, col in [(name1, smask1, 4), (name2, smask2, 4 + ncol)]:
    hp.mollview(
        pdust, sub=[nrow, ncol, col], max=1e-4, cmap="magma", coord="gc", title=""
    )
    hp.mollview(
        smask,
        sub=[nrow, ncol, col],
        min=0,
        max=1,
        cmap="bwr",
        cbar=False,
        title=name,
        reuse_axes=True,
        alpha=0.75 * smask,
    )
    if plot_bk:
        hp.mollview(
            sbkmask,
            sub=[nrow, ncol, col],
            min=0,
            max=1,
            cmap="plasma",
            cbar=False,
            title="",
            reuse_axes=True,
            alpha=0.75 * sbkmask,
        )

if plot_el:
    ax = fig.add_subplot(1, ncol, 5)
    for name, msum in [(name1, m1azelsum), (name2, m2azelsum)]:
        azs, els = hp.pix2ang(nside, np.arange(msum.size), lonlat=True)
        mean_el = np.sum(els * msum) / np.sum(msum)
        nbin = 30
        el = np.linspace(0, 90, nbin, endpoint=False)
        wbin = el[1] - el[0]
        el_hits = np.zeros(nbin)
        for pix_el, pix_hit in zip(els, msum):
            el_hits[int(pix_el / wbin)] += pix_hit
        el_hits /= np.sum(el_hits)
        good = el_hits != 0
        ax.step(el[good], el_hits[good], label=f"{name}, mean = {mean_el:.1f}°")
    ax.legend(loc="best")
    ax.set_xlabel("Elevation [deg]")
    ax.set_ylabel("Nhit")
    ax.set_title("Relative hits per elevation")

"""
for name, msum, col in [(name1, m1azelsum, 5), (name2, m2azelsum, 5 + ncol)]:
    msum /= np.amax(msum)
    msum[msum == 0] = hp.UNSEEN
    hp.projview(
        msum,
        sub=[nrow, ncol, col],
        cmap="magma",
        title=f"Az/El {name}",
        unit="Relative hits",
        rot=[180, 0],
        # rot_graticule=True,
        #latra=[30, 90],
        min=0,
        max=1,
        projection_type="cart",
        graticule=True,
        graticule_labels=True,
        latitude_grid_spacing=15,
        longitude_grid_spacing=45,
        phi_convention="counterclockwise",
    )
"""

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
