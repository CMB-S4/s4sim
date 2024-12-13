#  Generate a list of observing targets to match SAT depth

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


fname_out = "patches.txt"

cov = hp.read_map(
    "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2/noise_depth/sun90_f155_90years_cov.fits"
)
rhits = np.zeros_like(cov)
good = cov != 0
rhits[good] = 1 / cov[good]
rhits /= np.amax(rhits)

nside = hp.get_nside(cov)
npix = cov.size
pix = np.arange(npix)
lons, lats = hp.pix2ang(nside, pix, lonlat=True)

dlon = 5 # 10
dlat = 10 # 20
mask_sum = np.zeros(cov.size)
lines = []
for lat_min in range(-75, 75, 5):
    lat_max = lat_min + dlat
    for lon_min in range(0, 360, 5):
        lon_max = lon_min + dlon
        mask = np.ones(npix, dtype=bool)
        mask[lats < lat_min] = False
        mask[lats >= lat_max] = False
        if lon_max > 360:
            mask[np.logical_and(
                lons > lon_max - 360, lons < lon_min
            )] = False
        else:
            mask[lons < lon_min] = False
            mask[lons >= lon_max] = False
        hitmean = np.mean(rhits[mask])
        # if hitmean < .1:
        if hitmean < .20:
            continue
        # Lower priority number means more frequent targeting
        priority = 1 / hitmean
        if np.isnan(priority):
            raise RuntimeError()
        if lon_min >= 100 and lon_max <= 190:
            # Lower priority for the Northern tiles
            priority *= 1000
        elif lon_min >= 190 and lon_max <= 340:
            # Lowest the priority for tiles close to the Galaxy
            priority *= 1000000
        mask_sum += mask.astype(float)
        name = f"DEC{lat_min:+04}..{lat_max:+04}_RA{lon_min:+04}..{lon_max:+04}"
        lines.append(
            f"{name},{priority:.3f},{lon_min:.3f},{lat_max:.3f},{lon_max:.3f},{lat_min:.3f}"
        )

with open(fname_out, "w") as f:
    for line in lines:
        f.write("--patch\n")
        f.write(line + "\n")
print(f"Wrote {fname_out}")
