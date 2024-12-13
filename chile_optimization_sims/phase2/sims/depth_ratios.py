import glob
import os
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


# Find SAT field centers

radius = np.radians(10)

cov = hp.read_map("scaled_outputs/sun45_f155_90years_cov.fits")
bad = cov == 0
cov[bad] = 1e10
nside = hp.get_nside(cov)
npix = cov.size
pix = np.arange(npix)
lon, lat = hp.pix2ang(nside, pix, lonlat=True)
centers = {}
for lon_min, lon_max, name in [
        (100, 180, "Field2"),
        (0, 180, "Deep"),
        (280, 330, "Field3"),
        (180, 280, "Field4"),
]:
    temp = cov.copy()
    temp[lon < lon_min] = 1e10
    temp[lon > lon_max] = 1e10
    pix0 = np.argmin(temp)
    lon0, lat0 = lon[pix0], lat[pix0]
    vec0 = hp.ang2vec(lon0, lat0, lonlat=True)
    pixels = hp.query_disc(nside, vec0, radius)
    mask = np.zeros(npix, dtype=bool)
    mask[pixels] = True
    centers[name] = (lon0, lat0, mask)

vmin = np.amin(cov)
vmax = np.sort(cov)[int(npix * 0.03)]
cov[bad] = hp.UNSEEN
nrow, ncol = 1, 2
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
hp.mollview(cov, min=vmin, max=vmax, sub=[nrow, ncol, 1])
for name, (lon, lat, mask) in centers.items():
    # theta, phi = np.radians(90 - lat), np.radians(lon)
    hp.projplot(lon, lat, "ro", lonlat=True)
    hp.projtext(lon, lat + 10, name, lonlat=True)

hp.mollview(cov, min=vmin, max=vmax, sub=[nrow, ncol, 2])
mask_tot = np.zeros(npix, dtype=bool)
for name, (lon, lat, mask) in centers.items():
    mask_tot[mask] = True
    hp.projtext(lon, lat + 10, name, lonlat=True)
hp.mollview(
    mask_tot, reuse_axes=True, cbar=False, alpha=np.ones_like(cov) * .5
)

# Compute ratios
bands = ["f030", "f040", "f085", "f095", "f145", "f155", "f220", "f280"]
indir1 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase1/noise_depth"
indir2 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2/noise_depth"
print(" " * 8, end="")
for band in bands:
    print(f"{band:>8}", end="")
print("Phase2 / Phase1")
for name, (lon, lat, mask) in centers.items():
    print(f"{name:>8}", end="")
    for band in bands:
        m1 = hp.read_map(glob.glob(f"{indir1}/sat_{band}_*years_depth.fits")[0])
        m2 = hp.read_map(glob.glob(f"{indir2}/sun90_{band}_*years_depth.fits")[0])
        ratio = np.median(m2[mask] / m1[mask])
        print(f"{ratio:8.3f}", end="")
    print()

print("Sun45 / Sun90")
for name, (lon, lat, mask) in centers.items():
    print(f"{name:>8}", end="")
    for band in bands:
        m2 = hp.read_map(glob.glob(f"{indir2}/sun90_{band}_*years_depth.fits")[0])
        m3 = hp.read_map(glob.glob(f"{indir2}/sun45_{band}_*years_depth.fits")[0])
        ratio = np.median(m3[mask] / m2[mask])
        print(f"{ratio:8.3f}", end="")
    print()
