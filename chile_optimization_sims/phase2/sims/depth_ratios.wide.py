import glob
import os
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


# Find field centers

radius = np.radians(10)

cov = hp.read_map("scaled_outputs/lat_wide_f150_14years_cov.fits")
bad = cov == 0
cov[bad] = 1e10
cov_sorted = np.sort(cov)
limit = cov_sorted[int(0.60 * cov.size)]
mask = cov < limit
nside = hp.get_nside(cov)
npix = cov.size
pix = np.arange(npix)
lon, lat = hp.pix2ang(nside, pix, lonlat=True)

vmin = np.amin(cov)
vmax = limit
cov[bad] = hp.UNSEEN
nrow, ncol = 1, 2
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
hp.mollview(cov, min=vmin, max=vmax, sub=[nrow, ncol, 1])
#for name, (lon, lat, mask) in centers.items():
#    # theta, phi = np.radians(90 - lat), np.radians(lon)
#    hp.projplot(lon, lat, "ro", lonlat=True)
#    hp.projtext(lon, lat + 10, name, lonlat=True)

hp.mollview(cov, min=vmin, max=vmax, sub=[nrow, ncol, 2])
hp.mollview(
    mask, reuse_axes=True, cbar=False, alpha=np.ones_like(cov) * .5
)

# Compute ratios

bands = ["f030", "f040", "f090", "f150", "f220", "f280"]
indir1 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase1/noise_depth"
indir2 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2/noise_depth"
print(" " * 8, end="")
for band in bands:
    print(f"{band:>8}", end="")
print()
name = "wide"
print(f"{name:>8}", end="")
for band in bands:
    m1 = hp.read_map(f"{indir1}/lat_wide_{band}_14years_depth.fits")
    m2 = hp.read_map(f"{indir2}/lat_wide_{band}_14years_depth.fits")
    ratio = np.median(m2[mask] / m1[mask])
    print(f"{ratio:8.3f}", end="")
print()
